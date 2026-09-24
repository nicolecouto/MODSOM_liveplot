#!/usr/bin/env python3
"""
Fast MOD-SOM parser + realtime plotting (PyQt5 + pyqtgraph) WITH TCP input.

Inputs supported (mutually exclusive):
  --file PATH
  --serial PORT --baud 115200
  --tcp HOST PORT            (connect as TCP client)
  --tcp-listen PORT          (listen as TCP server; accepts ONE client)
  --watch-dir DIR            (tail whichever file is currently newest in DIR)

Speed design:
- Fixed-size circular buffers (NumPy) for all channels (no dynamic list growth).
- Batch append for EFE4/TTV where possible.
- GUI plots only last --window seconds (default 5s) from ring buffers (TTV, VNAV).
- EFE4 PSD uses a fixed sample count (--epsi-scan-length, default 1024), not a time
  window, so every update has the same number of frequency bins. The EFE4 time series
  show 3 scan lengths: the newest one, which the PSD is taken from, after a dashed
  line, and the previous 2 to its left for context.

Serial stop:
- When app exits, sends "som.stop\r\n" only if using --serial.

Notes:
- TCP is a byte stream (like serial): the state machine parses records exactly the same.
- If you want multiple TCP listeners, TCP is not ideal (use UDP or implement fan-out).
- --watch-dir never opens the serial port itself - it watches a directory some other
  process (which owns the port) is writing raw files into, and always tails whichever
  file in there is currently newest, following along as that process rotates to new
  files. This is how you can visualize live while another process is simultaneously
  logging the same stream to disk, without both processes contending for the port.
"""

import argparse
import glob
import os
import threading
import queue
import time
import struct
import socket
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from PyQt5 import QtCore, QtWidgets
import pyqtgraph as pg

try:
    import serial  # type: ignore
except ImportError:
    serial = None

# -----------------------------
# Protocol / tags
# -----------------------------
VALID_TAGS = {
    b"EFE4",
    b"TTV1",
    b"TTV2",
    b"TTV3",
    b"VNAV",
    b"SB49",
    b"SB41",
    b"ECOP",
    b"SOM3",
}

# -----------------------------
# Time conversion
# -----------------------------
MATLAB_EPOCH_DNUM = 719529.0
SECONDS_PER_DAY = 86400.0


def posix_to_matlab_dnum(posix_seconds: float) -> float:
    return posix_seconds / SECONDS_PER_DAY + MATLAB_EPOCH_DNUM


# -----------------------------
# ADC conversion helpers (EFE4)
# -----------------------------
ADC_FULL_SCALE_COUNTS = 2**24
ADC_FULL_SCALE_COUNTS_M1 = 2**23
ADC_VREF_TEMP = 2.5
ADC_VREF_SHEAR = 2.5
ADC_VREF_ACCEL = 1.8
ACC_OFFSET = 1.8 / 2
ACC_FACTOR = 0.4


def counts24_to_volts_unipolar(counts: np.ndarray, full_range: float) -> np.ndarray:
    return (counts.astype(np.float32) / ADC_FULL_SCALE_COUNTS) * full_range


def counts24_to_volts_bipolar(counts: np.ndarray, full_range: float) -> np.ndarray:
    return (counts.astype(np.float32) / ADC_FULL_SCALE_COUNTS_M1) * full_range


def volts_to_g(v: np.ndarray) -> np.ndarray:
    return (v - ACC_OFFSET) / ACC_FACTOR


# -----------------------------
# Circular buffer (fixed length)
# -----------------------------
class RingBuffer:
    """Fixed-size ring buffer for 1D arrays (single writer, multi-reader)."""
    def __init__(self, capacity: int, dtype=np.float64):
        self.capacity = int(capacity)
        self._x = np.empty(self.capacity, dtype=dtype)
        self._n = 0
        self._i = 0
        self._lock = threading.Lock()

    def append_many(self, data: np.ndarray):
        data = np.asarray(data)
        if data.size == 0:
            return
        with self._lock:
            m = int(data.size)
            if m >= self.capacity:
                data = data[-self.capacity:]
                self._x[:] = data
                self._n = self.capacity
                self._i = 0
                return

            end = self._i + m
            if end <= self.capacity:
                self._x[self._i:end] = data
            else:
                k = self.capacity - self._i
                self._x[self._i:] = data[:k]
                self._x[:end - self.capacity] = data[k:]
            self._i = end % self.capacity
            self._n = min(self.capacity, self._n + m)

    def append_one(self, v):
        with self._lock:
            self._x[self._i] = v
            self._i = (self._i + 1) % self.capacity
            self._n = min(self.capacity, self._n + 1)

    def snapshot_last_seconds(self, tbuf: "RingBuffer", seconds: float) -> Tuple[np.ndarray, np.ndarray]:
        seconds = float(seconds)
        with tbuf._lock, self._lock:
            n = min(tbuf._n, self._n)
            if n <= 1:
                return np.empty(0), np.empty(0)
            t = _ordered_copy(tbuf._x, tbuf._i, n)
            y = _ordered_copy(self._x, self._i, n)

        t0 = t[-1] - seconds
        j = np.searchsorted(t, t0, side="left")
        return t[j:], y[j:]

    def snapshot_last_n(self, tbuf: "RingBuffer", n: int) -> Tuple[np.ndarray, np.ndarray]:
        """Like snapshot_last_seconds, but a fixed sample count instead of a time window -
        every call returns the same length (once enough samples exist), so the resulting
        PSD's frequency bins don't drift update to update."""
        n = int(n)
        with tbuf._lock, self._lock:
            m = min(tbuf._n, self._n, n)
            if m <= 1:
                return np.empty(0), np.empty(0)
            t = _ordered_copy(tbuf._x, tbuf._i, m)
            y = _ordered_copy(self._x, self._i, m)
        return t, y

    def size(self) -> int:
        with self._lock:
            return self._n


def _ordered_copy(x: np.ndarray, write_index: int, n: int) -> np.ndarray:
    if n <= 0:
        return np.empty(0, dtype=x.dtype)
    start = (write_index - n) % x.size
    if start < write_index:
        return x[start:write_index].copy()
    else:
        return np.concatenate((x[start:], x[:write_index])).copy()


# -----------------------------
# Record emitted by state machine
# -----------------------------
@dataclass
class Record:
    inst_tag: str
    posix: float
    dnum: float
    payload_size: int
    payload: bytes


# -----------------------------
# Instrument buffers
# -----------------------------
class EFE4Buffers:
    CHANNELS = ("t1", "t2", "s1", "s2", "a1", "a2", "a3")

    def __init__(self, capacity: int):
        # Held by the writer around appending all eight buffers, and by
        # snapshot_last_n, so a snapshot never sees some channels a record ahead.
        self.lock = threading.Lock()
        self.t = RingBuffer(capacity, dtype=np.float64)
        self.t1 = RingBuffer(capacity, dtype=np.float32)
        self.t2 = RingBuffer(capacity, dtype=np.float32)
        self.s1 = RingBuffer(capacity, dtype=np.float32)
        self.s2 = RingBuffer(capacity, dtype=np.float32)
        self.a1 = RingBuffer(capacity, dtype=np.float32)
        self.a2 = RingBuffer(capacity, dtype=np.float32)
        self.a3 = RingBuffer(capacity, dtype=np.float32)

    def snapshot_last_n(self, n: int) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """Last n samples of time and every channel, all the same length and aligned."""
        with self.lock:
            t = np.empty(0)
            ch = {}
            for name in self.CHANNELS:
                t, ch[name] = getattr(self, name).snapshot_last_n(self.t, n)
        return t, ch


class TTVBuffers:
    def __init__(self, capacity: int):
        self.t = RingBuffer(capacity, dtype=np.float64)
        self.tof_up = RingBuffer(capacity, dtype=np.float32)
        self.tof_down = RingBuffer(capacity, dtype=np.float32)
        self.dtof = RingBuffer(capacity, dtype=np.float32)
        self.error = RingBuffer(capacity, dtype=np.float32)
        self.up_peak = RingBuffer(capacity, dtype=np.float32)
        self.down_peak = RingBuffer(capacity, dtype=np.float32)

        # Derived per-tag beams (each tag fills only one; others get NaN)
        self.beam1 = RingBuffer(capacity, dtype=np.float32)
        self.beam2 = RingBuffer(capacity, dtype=np.float32)
        self.beam3 = RingBuffer(capacity, dtype=np.float32)


class VNAVBuffers:
    def __init__(self, capacity: int):
        self.t = RingBuffer(capacity, dtype=np.float64)
        self.mag_x = RingBuffer(capacity, dtype=np.float32)
        self.mag_y = RingBuffer(capacity, dtype=np.float32)
        self.mag_z = RingBuffer(capacity, dtype=np.float32)
        self.accel_x = RingBuffer(capacity, dtype=np.float32)
        self.accel_y = RingBuffer(capacity, dtype=np.float32)
        self.accel_z = RingBuffer(capacity, dtype=np.float32)
        self.gyro_x = RingBuffer(capacity, dtype=np.float32)
        self.gyro_y = RingBuffer(capacity, dtype=np.float32)
        self.gyro_z = RingBuffer(capacity, dtype=np.float32)


# -----------------------------
# Input abstractions
# -----------------------------
class SocketSource:
    """A minimal recv() wrapper to look like a 'source' for ByteSourceThread."""
    def __init__(self, sock: socket.socket):
        self.sock = sock

    def read(self, n: int) -> bytes:
        # recv() returns b'' on close
        return self.sock.recv(n)

    def close(self):
        try:
            self.sock.shutdown(socket.SHUT_RDWR)
        except Exception:
            pass
        try:
            self.sock.close()
        except Exception:
            pass


class TCPListenThread(threading.Thread):
    """
    TCP server that accepts ONE client and then hands the connected socket to main thread.
    """
    def __init__(self, host: str, port: int, stop_event: threading.Event):
        super().__init__(daemon=True)
        self.host = host
        self.port = int(port)
        self.stop_event = stop_event
        self._ready = threading.Event()
        self.client_sock: Optional[socket.socket] = None
        self._err: Optional[str] = None

    def run(self):
        srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            srv.bind((self.host, self.port))
            srv.listen(1)
            srv.settimeout(0.2)
            print(f"[TCP] Listening on {self.host}:{self.port} ...")
            while not self.stop_event.is_set():
                try:
                    c, addr = srv.accept()
                    print(f"[TCP] Client connected from {addr}")
                    c.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
                    self.client_sock = c
                    self._ready.set()
                    return
                except socket.timeout:
                    continue
        except Exception as e:
            self._err = str(e)
            self._ready.set()
        finally:
            try:
                srv.close()
            except Exception:
                pass

    def wait_client(self, timeout: Optional[float] = None) -> Optional[socket.socket]:
        self._ready.wait(timeout=timeout)
        if self._err:
            raise RuntimeError(f"TCP listen failed: {self._err}")
        return self.client_sock


def _list_candidates(watch_dir: str, pattern: str) -> List[str]:
    """Every file in watch_dir matching pattern, oldest to newest (by mtime, filename as
    tiebreaker). The caller wants candidates[-1] - "the current newest file"."""
    entries = []
    for p in glob.glob(os.path.join(watch_dir, pattern)):
        try:
            st = os.stat(p)
        except OSError:
            continue
        entries.append((st.st_mtime, p))
    entries.sort(key=lambda e: (e[0], e[1]))
    return [p for _, p in entries]


class DirWatchSource:
    """A '.read(n) -> bytes' source (same shape as SocketSource) that watches watch_dir
    and always tails whichever matching file is currently newest - never any other file.

    This deliberately does not try to preserve every byte across a rotation: if a newer
    file appears while there's still unread data in the old one, that leftover is simply
    dropped. The point is "always looking at the most recent data" for a live view, not a
    complete byte-exact record - something else (the process actually saving this data,
    outside this repo) already owns that.
    """

    def __init__(
        self,
        watch_dir: str,
        pattern: str,
        poll_interval_s: float,
        stop_event: threading.Event,
        seed_bytes: int,
    ):
        self.watch_dir = watch_dir
        self.pattern = pattern
        self.poll_interval_s = float(poll_interval_s)
        self.stop_event = stop_event
        self.seed_bytes = int(seed_bytes)
        self.current_path: Optional[str] = None
        self.current_fh = None
        self._seeded = False

    def _open(self, path: str, seed: bool):
        fh = open(path, "rb")
        if seed and self.seed_bytes > 0:
            size = os.fstat(fh.fileno()).st_size
            if size > self.seed_bytes:
                fh.seek(size - self.seed_bytes)
        self.current_path = path
        self.current_fh = fh

    def _replaced_on_disk(self) -> bool:
        """True if current_path was truncated+rewritten in place (same name, new inode) or
        removed, rather than rotated to a new filename."""
        if self.current_path is None or self.current_fh is None:
            return False
        try:
            disk_ino = os.stat(self.current_path).st_ino
        except OSError:
            return True
        try:
            open_ino = os.fstat(self.current_fh.fileno()).st_ino
        except OSError:
            return True
        return disk_ino != open_ino

    def read(self, n: int) -> bytes:
        while not self.stop_event.is_set():
            if self.current_fh is None:
                candidates = _list_candidates(self.watch_dir, self.pattern)
                if not candidates:
                    time.sleep(self.poll_interval_s)
                    continue
                self._open(candidates[-1], seed=not self._seeded)
                self._seeded = True
                continue

            data = self.current_fh.read(n)
            if data:
                return data

            # Caught up to this file's current EOF.
            if self._replaced_on_disk():
                self.current_fh.close()
                self.current_fh = None
                continue

            candidates = _list_candidates(self.watch_dir, self.pattern)
            if candidates and candidates[-1] != self.current_path:
                self.current_fh.close()
                self._open(candidates[-1], seed=False)
                continue

            time.sleep(self.poll_interval_s)

        return b""

    def close(self):
        if self.current_fh is not None:
            try:
                self.current_fh.close()
            except Exception:
                pass


# -----------------------------
# Reader thread: file/serial/tcp
# -----------------------------
class ByteSourceThread(threading.Thread):
    def __init__(
        self,
        source,
        is_stream: bool,
        out_queue: queue.Queue,
        stop_event: threading.Event,
        chunk_size: int = 4096,
        start_command: Optional[bytes] = b"som.start\r\n",
        start_delay_s: float = 0.05,
    ):
        """
        source must implement .read(n)->bytes
        If serial-like and you want to send som.start, pass start_command; for TCP/file pass start_command=None.
        """
        super().__init__(daemon=True)
        self.source = source
        self.is_stream = is_stream
        self.out_queue = out_queue
        self.stop_event = stop_event
        self.chunk_size = chunk_size
        self.start_command = start_command
        self.start_delay_s = start_delay_s

    def run(self):
        try:
            # Start command only if source supports .write (serial)
            if self.start_command and hasattr(self.source, "write"):
                try:
                    if hasattr(self.source, "reset_input_buffer"):
                        self.source.reset_input_buffer()
                    self.source.write(self.start_command)
                    if hasattr(self.source, "flush"):
                        self.source.flush()
                    if self.start_delay_s > 0:
                        time.sleep(self.start_delay_s)
                except Exception as e:
                    print(f"[ByteSourceThread] start_command failed: {e}")

            while not self.stop_event.is_set():
                data = self.source.read(self.chunk_size)
                if not data:
                    # file: EOF; stream: remote closed
                    break
                self.out_queue.put(data)

            self.out_queue.put(None)
        except Exception as e:
            print(f"[ByteSourceThread] Error: {e}")
            self.out_queue.put(None)


# -----------------------------
# State-machine parser
# -----------------------------
class EpsiStateMachineParser:
    STATE_SYNC = 0
    STATE_TAG = 1
    STATE_HEADER = 2
    STATE_PAYLOAD = 3
    HEADER_LEN = 1 + 4 + 16 + 8 + 1 + 2

    def __init__(self, record_queue: queue.Queue):
        self.buffer = bytearray()
        self.state = self.STATE_SYNC
        self.current_tag: Optional[bytes] = None
        self.current_timestamp_ms: Optional[int] = None
        self.current_payload_size: Optional[int] = None
        self.record_queue = record_queue

    def feed(self, data: bytes):
        self.buffer.extend(data)
        while True:
            if self.state == self.STATE_SYNC:
                if not self._sync():
                    break
            elif self.state == self.STATE_TAG:
                if not self._tag():
                    break
            elif self.state == self.STATE_HEADER:
                if not self._header():
                    break
            elif self.state == self.STATE_PAYLOAD:
                if not self._payload():
                    break
            else:
                self.state = self.STATE_SYNC

    def _sync(self) -> bool:
        idx = self.buffer.find(b"$")
        if idx == -1:
            self.buffer.clear()
            return False
        if idx > 0:
            del self.buffer[:idx]
        if len(self.buffer) < 1 + 4:
            return False
        self.state = self.STATE_TAG
        return True

    def _tag(self) -> bool:
        if len(self.buffer) < 1 + 4:
            return False
        if self.buffer[0] != ord("$"):
            del self.buffer[0]
            self.state = self.STATE_SYNC
            return True
        tag_bytes = bytes(self.buffer[1:5])
        if tag_bytes not in VALID_TAGS:
            del self.buffer[0]
            self.state = self.STATE_SYNC
            return True
        self.current_tag = tag_bytes
        self.state = self.STATE_HEADER
        return True

    def _header(self) -> bool:
        if len(self.buffer) < self.HEADER_LEN:
            return False
        header = self.buffer[:self.HEADER_LEN]
        if header[0] != ord("$") or header[29] != ord("*"):
            del self.buffer[0]
            self.state = self.STATE_SYNC
            return True

        ts_hex = header[5:21].decode("ascii", errors="ignore")
        size_hex = header[21:29].decode("ascii", errors="ignore")
        cksum_hex = header[30:32].decode("ascii", errors="ignore")
        try:
            timestamp_ms = int(ts_hex, 16)
            payload_size = int(size_hex, 16)
            published = int(cksum_hex, 16)
        except ValueError:
            del self.buffer[0]
            self.state = self.STATE_SYNC
            return True

        computed = 0
        for b in header[0:29]:
            computed ^= b
        computed &= 0xFF

        if computed == published:
            self.current_timestamp_ms = timestamp_ms
            self.current_payload_size = payload_size
            del self.buffer[:self.HEADER_LEN]
            self.state = self.STATE_PAYLOAD
        else:
            del self.buffer[0]
            self.state = self.STATE_SYNC
        return True

    def _payload(self) -> bool:
        if self.current_payload_size is None or self.current_tag is None or self.current_timestamp_ms is None:
            self.state = self.STATE_SYNC
            return True

        needed = self.current_payload_size + 1 + 2
        if len(self.buffer) < needed:
            return False

        payload = self.buffer[: self.current_payload_size]
        star = self.buffer[self.current_payload_size]
        cksum_bytes = self.buffer[self.current_payload_size + 1: self.current_payload_size + 3]
        if star != ord("*"):
            del self.buffer[0]
            self._reset()
            self.state = self.STATE_SYNC
            return True

        try:
            published = int(cksum_bytes.decode("ascii", errors="ignore"), 16)
        except ValueError:
            del self.buffer[0]
            self._reset()
            self.state = self.STATE_SYNC
            return True

        computed = 0
        for b in payload:
            computed ^= b
        computed &= 0xFF

        if computed == published:
            tag_str = self.current_tag.decode("ascii")
            posix = self.current_timestamp_ms / 1000.0
            dnum = posix_to_matlab_dnum(posix)
            self.record_queue.put(
                Record(tag_str, posix, dnum, self.current_payload_size, bytes(payload))
            )

        del self.buffer[:needed]
        self._reset()
        self.state = self.STATE_SYNC
        return True

    def _reset(self):
        self.current_tag = None
        self.current_timestamp_ms = None
        self.current_payload_size = None


# -----------------------------
# Parser thread: bytes -> records
# -----------------------------
class ParserThread(threading.Thread):
    def __init__(self, byte_queue: queue.Queue, parser: EpsiStateMachineParser, record_queue: queue.Queue, stop_event: threading.Event):
        super().__init__(daemon=True)
        self.byte_queue = byte_queue
        self.parser = parser
        self.record_queue = record_queue
        self.stop_event = stop_event

    def run(self):
        while not self.stop_event.is_set():
            chunk = self.byte_queue.get()
            if chunk is None:
                break
            self.parser.feed(chunk)
        self.record_queue.put(None)


# -----------------------------
# Record processing thread: records -> ring buffers
# -----------------------------
class RecordProcessorThread(threading.Thread):
    L = 0.0382
    TTV_PACKET_SIZE = 35

    def __init__(self, record_queue: queue.Queue, stop_event: threading.Event, buffers: Dict[str, object]):
        super().__init__(daemon=True)
        self.record_queue = record_queue
        self.stop_event = stop_event
        self.buffers = buffers

    def run(self):
        while not self.stop_event.is_set():
            rec = self.record_queue.get()
            if rec is None:
                break
            try:
                self._process(rec)
            except Exception as e:
                print(f"[RecordProcessor] Error {rec.inst_tag}: {e}")

    def _process(self, rec: Record):
        tag = rec.inst_tag
        if tag == "EFE4":
            self._parse_efe4(rec, self.buffers["EFE4"])  # type: ignore
        elif tag in ("TTV1", "TTV2", "TTV3"):
            self._parse_ttv(rec, self.buffers[tag], tag)  # type: ignore
        elif tag == "VNAV":
            self._parse_vnav(rec, self.buffers["VNAV"])  # type: ignore

    def _parse_efe4(self, rec: Record, b: EFE4Buffers):
        payload = rec.payload
        sample_bytes = 8 + 3 * 7  # 29
        n = len(payload) // sample_bytes
        if n <= 0:
            return

        ts = np.empty(n, dtype=np.float64)
        c = np.empty((n, 7), dtype=np.uint32)

        off = 0
        for i in range(n):
            ts_ms = int.from_bytes(payload[off:off+8], "little", signed=False)
            ts[i] = ts_ms / 1000.0
            off += 8
            for k in range(7):
                c[i, k] = int.from_bytes(payload[off:off+3], "big", signed=False)
                off += 3

        t1 = counts24_to_volts_unipolar(c[:, 0], ADC_VREF_TEMP)
        t2 = counts24_to_volts_unipolar(c[:, 1], ADC_VREF_TEMP)
        s1 = counts24_to_volts_bipolar(c[:, 2], ADC_VREF_SHEAR)
        s2 = counts24_to_volts_bipolar(c[:, 3], ADC_VREF_SHEAR)
        a1 = volts_to_g(counts24_to_volts_unipolar(c[:, 4], ADC_VREF_ACCEL))
        a2 = volts_to_g(counts24_to_volts_unipolar(c[:, 5], ADC_VREF_ACCEL))
        a3 = volts_to_g(counts24_to_volts_unipolar(c[:, 6], ADC_VREF_ACCEL))

        with b.lock:
            b.t.append_many(ts)
            b.t1.append_many(t1); b.t2.append_many(t2)
            b.s1.append_many(s1); b.s2.append_many(s2)
            b.a1.append_many(a1); b.a2.append_many(a2); b.a3.append_many(a3)

    def _parse_ttv(self, rec: Record, b: TTVBuffers, tag: str):
        payload = rec.payload
        ps = self.TTV_PACKET_SIZE
        n = len(payload) // ps
        if n <= 0:
            return

        ts = np.empty(n, dtype=np.float64)
        tof_up = np.empty(n, dtype=np.float32)
        tof_dn = np.empty(n, dtype=np.float32)
        dtof = np.empty(n, dtype=np.float32)
        err = np.empty(n, dtype=np.float32)
        up_pk = np.empty(n, dtype=np.float32)
        dn_pk = np.empty(n, dtype=np.float32)

        off = 0
        for i in range(n):
            ts_hex = payload[off:off+16]
            try:
                ts_ms = int(ts_hex.decode("ascii", errors="ignore"), 16)
            except ValueError:
                ts[i] = np.nan
                tof_up[i] = np.nan
                tof_dn[i] = np.nan
                dtof[i] = np.nan
                err[i] = np.nan
                up_pk[i] = np.nan
                dn_pk[i] = np.nan
                off += ps
                continue

            ts[i] = ts_ms / 1000.0
            off += 16
            tof_up[i] = struct.unpack(">f", payload[off:off+4])[0]; off += 4
            tof_dn[i] = struct.unpack(">f", payload[off:off+4])[0]; off += 4
            dtof[i] = struct.unpack(">f", payload[off:off+4])[0]; off += 4
            err[i] = float(payload[off]); off += 1
            up_pk[i] = float(struct.unpack(">H", payload[off:off+2])[0]); off += 2
            dn_pk[i] = float(struct.unpack(">H", payload[off:off+2])[0]); off += 2
            off += 2  # CRLF

        good = np.isfinite(ts)
        if not np.any(good):
            return

        ts = ts[good]
        tof_up = tof_up[good]
        tof_dn = tof_dn[good]
        dtof = dtof[good]
        err = err[good]
        up_pk = up_pk[good]
        dn_pk = dn_pk[good]

        denom = (tof_up.astype(np.float64) * tof_dn.astype(np.float64))
        with np.errstate(divide="ignore", invalid="ignore"):
            v = (self.L / 2.0) * (dtof.astype(np.float64) / denom)
        v = v.astype(np.float32)

        b.t.append_many(ts)
        b.tof_up.append_many(tof_up)
        b.tof_down.append_many(tof_dn)
        b.dtof.append_many(dtof)
        b.error.append_many(err)
        b.up_peak.append_many(up_pk)
        b.down_peak.append_many(dn_pk)

        nan_fill = np.full(v.shape, np.nan, dtype=np.float32)
        if tag == "TTV1":
            b.beam1.append_many(v); b.beam2.append_many(nan_fill); b.beam3.append_many(nan_fill)
        elif tag == "TTV2":
            b.beam1.append_many(nan_fill); b.beam2.append_many(v); b.beam3.append_many(nan_fill)
        else:
            b.beam1.append_many(nan_fill); b.beam2.append_many(nan_fill); b.beam3.append_many(v)

    def _parse_vnav(self, rec: Record, b: VNAVBuffers):
        payload = rec.payload
        nbytes = len(payload)
        i = 0
        while i + 16 + 1 < nbytes:
            ts_bytes = payload[i:i+16]
            try:
                ts_ms = int(ts_bytes.decode("ascii", errors="ignore"), 16)
            except ValueError:
                i += 1
                continue
            if payload[i+16] != ord("$"):
                i += 1
                continue
            star = payload.find(b"*", i+16)
            if star == -1 or star + 2 >= nbytes:
                break
            msg = payload[i+16:star].decode("ascii", errors="ignore")
            fields = msg.split(",")
            if len(fields) < 10:
                i = star + 3
                continue
            try:
                magx = float(fields[1]); magy = float(fields[2]); magz = float(fields[3])
                ax = float(fields[4]); ay = float(fields[5]); az = float(fields[6])
                gx = float(fields[7]); gy = float(fields[8]); gz = float(fields[9])
            except ValueError:
                i = star + 3
                continue

            t = ts_ms / 1000.0
            b.t.append_one(t)
            b.mag_x.append_one(magx); b.mag_y.append_one(magy); b.mag_z.append_one(magz)
            b.accel_x.append_one(ax); b.accel_y.append_one(ay); b.accel_z.append_one(az)
            b.gyro_x.append_one(gx); b.gyro_y.append_one(gy); b.gyro_z.append_one(gz)

            i = star + 3
            while i < nbytes and payload[i] in (0x0D, 0x0A):
                i += 1


# -----------------------------
# PSD helper (fast)
# -----------------------------
def psd_fft(t: np.ndarray, y: np.ndarray) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    if t.size < 32:
        return None, None
    dt = np.median(np.diff(t))
    if not np.isfinite(dt) or dt <= 0:
        return None, None
    fs = 1.0 / dt
    n = y.size
    w = np.hanning(n)
    w2 = (w * w).sum()
    yd = y - np.mean(y)
    Y = np.fft.rfft(yd * w)
    psd = (np.abs(Y) ** 2) / (fs * w2)
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    psd[psd <= 0] = 1e-30
    return f, psd


# -----------------------------
# Colors
# -----------------------------
EFE4_COLORS = {
    "t1": (255, 0, 0),
    "t2": (255, 165, 0),
    "s1": (0, 0, 255),
    "s2": (0, 255, 255),
    "a1": (0, 255, 0),
    "a2": (255, 0, 255),
    "a3": (139, 69, 19),
}
EFE4_LEFT_AXIS_WIDTH = 90  # px, EFE4 time-series y-axes
EFE4_PSD_YMIN = 1e-18      # lowest y shown on the EFE4 PSD panel
EPOCH_MIN = 1e9            # timestamps at/after this (2001-09-09) are treated as real UTC
TTV_TAG_COLORS = {"TTV1": (255, 0, 0), "TTV2": (0, 255, 0), "TTV3": (0, 0, 255)}
VNAV_COLORS = {
    "mag_x": (255, 0, 0), "mag_y": (0, 255, 0), "mag_z": (0, 0, 255),
    "accel_x": (255, 0, 255), "accel_y": (0, 255, 255), "accel_z": (255, 255, 0),
    "gyro_x": (255, 165, 0), "gyro_y": (128, 0, 128), "gyro_z": (0, 128, 0),
}

# -----------------------------
# GUI Windows
# -----------------------------
class TimeAxis(pg.AxisItem):
    """Bottom axis for x in seconds after t0 (posix seconds). Labels ticks as UTC
    HH:MM:SS when `absolute` is set, otherwise as plain seconds."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.t0 = 0.0
        self.absolute = False

    def tickStrings(self, values, scale, spacing):
        if not self.absolute:
            return super().tickStrings(values, scale, spacing)
        out = []
        for v in values:
            s = datetime.fromtimestamp(self.t0 + v, tz=timezone.utc).strftime("%H:%M:%S.%f")
            out.append(s[:-7] if spacing >= 1 else s[:-5])  # tenths only if ticks are sub-second
        return out


class EFE4Window(QtWidgets.QMainWindow):
    def __init__(self, buf: EFE4Buffers, epsi_scan_length: int = 1024):
        super().__init__()
        self.buf = buf
        self.epsi_scan_length = int(epsi_scan_length)
        self.setWindowTitle(f"EFE4 (fast) - epsi scan length {self.epsi_scan_length}")

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)
        left = QtWidgets.QVBoxLayout()
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 1)
        layout.addLayout(right, 1)

        # One panel per channel (a2/a3 share one), so each y-axis auto-scales to its own
        # channel instead of being stretched by an offset between channels. Each shows
        # 3 scan lengths; the PSD uses the newest one, to the right of a dashed line.
        # X is linked across panels and set explicitly every update (mouse x zoom off),
        # so all six always have the same limits; only the bottom panel shows time ticks.
        panels = [("t1 [V]", ["t1"]), ("t2 [V]", ["t2"]),
                  ("s1 [V]", ["s1"]), ("s2 [V]", ["s2"]),
                  ("a1 [g]", ["a1"]), (None, ["a2", "a3"])]
        self.curves = {}
        self.scan_lines = []  # one per panel, at the start of the PSD scan
        self.time_axis = TimeAxis(orientation="bottom")
        first = None
        scan_pen = pg.mkPen((200, 200, 200), style=QtCore.Qt.DashLine)
        for i, (ylabel, names) in enumerate(panels):
            last = i == len(panels) - 1
            p = pg.PlotWidget(axisItems={"bottom": self.time_axis} if last else None)
            p.setMouseEnabled(x=False, y=True)
            if len(names) > 1:
                # Color-code the channel names in the axis label instead of a legend,
                # which would cover data in a panel this short.
                ylabel = ", ".join(
                    "<span style='color:#%02x%02x%02x'>%s</span>" % (*EFE4_COLORS[n], n)
                    for n in names
                ) + " [g]"
            p.setLabel("left", ylabel)
            # Fixed axis width keeps panels aligned and leaves room for long tick labels
            # (e.g. "2.500013"); hideOverlappingLabels drops a tick label at the top or
            # bottom edge rather than drawing it half clipped.
            ax = p.getAxis("left")
            ax.setWidth(EFE4_LEFT_AXIS_WIDTH)
            ax.setStyle(hideOverlappingLabels=True)
            p.enableAutoRange(axis="y")
            p.setAutoVisible(y=True)
            if first is None:
                first = p
            else:
                p.setXLink(first)
            if not last:
                p.getAxis("bottom").setStyle(showValues=False)
            for name in names:
                self.curves[name] = p.plot([], [], pen=pg.mkPen(EFE4_COLORS[name]), name=name)
            line = pg.InfiniteLine(angle=90, pen=scan_pen)
            p.addItem(line, ignoreBounds=True)
            self.scan_lines.append(line)
            left.addWidget(p)
        self.p_ts_first = first
        self._time_label = None

        self.scan_info_label = QtWidgets.QLabel("scan center: -- | duration: -- | length: --")
        right.addWidget(self.scan_info_label)

        self.p_psd = pg.PlotWidget()
        self.p_psd.setLogMode(x=True, y=True)
        # Log mode: view coordinates are log10, so this floors the y-axis at 1e-18
        # while leaving the top free to auto-scale.
        self.p_psd.setLimits(yMin=np.log10(EFE4_PSD_YMIN))
        self.p_psd.addLegend()
        right.addWidget(self.p_psd)
        self.psd = {
            "t1": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["t1"]), name="t1"),
            "t2": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["t2"]), name="t2"),
            "s1": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["s1"]), name="s1"),
            "s2": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["s2"]), name="s2"),
            "a1": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["a1"]), name="a1"),
            "a2": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["a2"]), name="a2"),
            "a3": self.p_psd.plot([], [], pen=pg.mkPen(EFE4_COLORS["a3"]), name="a3"),
        }
        pen24 = pg.mkPen((150, 150, 150), style=QtCore.Qt.DashLine)
        pen20 = pg.mkPen((200, 100, 100), style=QtCore.Qt.DashLine)
        pen16 = pg.mkPen((100, 200, 100), style=QtCore.Qt.DashLine)
        self.nf24 = self.p_psd.plot([], [], pen=pen24, name="24-bit NF")
        self.nf20 = self.p_psd.plot([], [], pen=pen20, name="20-bit NF")
        self.nf16 = self.p_psd.plot([], [], pen=pen16, name="16-bit NF")

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plots)
        self.timer.start(100)

    def update_plots(self):
        # One aligned snapshot of 3 scan lengths feeds both the time series and the PSD.
        # The PSD uses a fixed sample count (epsi_scan_length), not a time window, so
        # every update gets the same number of frequency bins.
        n = self.epsi_scan_length
        t, ch = self.buf.snapshot_last_n(3 * n)
        if t.size < 32:
            return
        # The scan is the newest scan length; the 2 before it are context. Until 3 scan
        # lengths have arrived (startup) there's less context to its left.
        i0 = max(0, t.size - n)

        # x is seconds after a whole second, so integer ticks land on whole clock seconds.
        t0 = np.floor(t[0])
        absolute = t0 >= EPOCH_MIN
        self.time_axis.t0 = t0
        self.time_axis.absolute = absolute
        time_label = "time [UTC]" if absolute else "time [s]"
        if time_label != self._time_label:
            self.time_axis.setLabel(time_label)
            self._time_label = time_label

        x = t - t0
        for name, arr in ch.items():
            self.curves[name].setData(x, arr)
        self.p_ts_first.setXRange(x[0], x[-1], padding=0)
        for line in self.scan_lines:
            line.setValue(x[i0])

        pt = t[i0:]
        if pt.size >= 32:
            for name, arr in ch.items():
                f, p = psd_fft(pt, arr[i0:].astype(np.float64))
                if f is not None:
                    self.psd[name].setData(f, p)

            center = (pt[0] + pt[-1]) / 2.0
            if absolute:
                center_str = datetime.fromtimestamp(center, tz=timezone.utc).strftime(
                    "%Y-%m-%d %H:%M:%S.%f")[:-3] + " UTC"
            else:
                center_str = f"{center:.3f} s"
            duration_s = pt[-1] - pt[0]
            self.scan_info_label.setText(
                f"scan center: {center_str} | "
                f"duration: {duration_s:.3f} s | length: {pt.size} samples"
            )

        f = np.logspace(-1, 2, 200)
        self.nf24.setData(f, np.full_like(f, 1e-12))
        self.nf20.setData(f, np.full_like(f, 1e-10))
        self.nf16.setData(f, np.full_like(f, 1e-8))


class TTVWindow(QtWidgets.QMainWindow):
    def __init__(self, ttv: Dict[str, TTVBuffers], window_seconds: float = 5.0):
        super().__init__()
        self.ttv = ttv
        self.window_seconds = float(window_seconds)
        self.setWindowTitle("TTV1/2/3 (fast)")

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)
        left = QtWidgets.QVBoxLayout()
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 1)
        layout.addLayout(right, 1)

        self.ts_plots: Dict[str, pg.PlotWidget] = {}
        self.ts_curves: Dict[Tuple[str, str], pg.PlotDataItem] = {}

        for ch in ["tof_down", "tof_up", "dtof", "error", "up_peak", "down_peak", "beamV"]:
            p = pg.PlotWidget(); p.addLegend()
            left.addWidget(p)
            self.ts_plots[ch] = p
            for tag in ("TTV1", "TTV2", "TTV3"):
                col = TTV_TAG_COLORS[tag]
                self.ts_curves[(tag, ch)] = p.plot([], [], pen=pg.mkPen(col), name=tag)

        self.p_psd = pg.PlotWidget()
        self.p_psd.setLogMode(x=True, y=True)
        self.p_psd.addLegend()
        right.addWidget(self.p_psd)
        self.psd_curves: Dict[Tuple[str, str], pg.PlotDataItem] = {}
        for tag in ("TTV1", "TTV2", "TTV3"):
            col = TTV_TAG_COLORS[tag]
            for ch in ["tof_down", "tof_up", "dtof"]:
                self.psd_curves[(tag, ch)] = self.p_psd.plot([], [], pen=pg.mkPen(col), name=f"{tag}-{ch}")

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plots)
        self.timer.start(150)

    def update_plots(self):
        for tag, b in self.ttv.items():
            t, y_dn = b.tof_down.snapshot_last_seconds(b.t, self.window_seconds)
            if t.size < 32:
                continue
            tr = t - t[0]
            _, y_up = b.tof_up.snapshot_last_seconds(b.t, self.window_seconds)
            _, y_dt = b.dtof.snapshot_last_seconds(b.t, self.window_seconds)
            _, y_er = b.error.snapshot_last_seconds(b.t, self.window_seconds)
            _, y_uk = b.up_peak.snapshot_last_seconds(b.t, self.window_seconds)
            _, y_dk = b.down_peak.snapshot_last_seconds(b.t, self.window_seconds)

            if tag == "TTV1":
                _, vb = b.beam1.snapshot_last_seconds(b.t, self.window_seconds)
            elif tag == "TTV2":
                _, vb = b.beam2.snapshot_last_seconds(b.t, self.window_seconds)
            else:
                _, vb = b.beam3.snapshot_last_seconds(b.t, self.window_seconds)

            self.ts_curves[(tag, "tof_down")].setData(tr, y_dn)
            self.ts_curves[(tag, "tof_up")].setData(tr, y_up)
            self.ts_curves[(tag, "dtof")].setData(tr, y_dt)
            self.ts_curves[(tag, "error")].setData(tr, y_er)
            self.ts_curves[(tag, "up_peak")].setData(tr, y_uk)
            self.ts_curves[(tag, "down_peak")].setData(tr, y_dk)
            self.ts_curves[(tag, "beamV")].setData(tr, vb)

            for chname, arr in [("tof_down", y_dn), ("tof_up", y_up), ("dtof", y_dt)]:
                f, p = psd_fft(t, arr.astype(np.float64))
                if f is not None:
                    self.psd_curves[(tag, chname)].setData(f, p)


class VNAVWindow(QtWidgets.QMainWindow):
    def __init__(self, b: VNAVBuffers, window_seconds: float = 5.0):
        super().__init__()
        self.b = b
        self.window_seconds = float(window_seconds)
        self.setWindowTitle("VNAV (fast)")

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)
        left = QtWidgets.QVBoxLayout()
        right = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 1)
        layout.addLayout(right, 1)

        self.p_mag = pg.PlotWidget(); self.p_mag.addLegend(); left.addWidget(self.p_mag)
        self.p_acc = pg.PlotWidget(); self.p_acc.addLegend(); left.addWidget(self.p_acc)
        self.p_gyr = pg.PlotWidget(); self.p_gyr.addLegend(); left.addWidget(self.p_gyr)

        self.mag = {k: self.p_mag.plot([], [], pen=pg.mkPen(VNAV_COLORS[k]), name=k) for k in ("mag_x","mag_y","mag_z")}
        self.acc = {k: self.p_acc.plot([], [], pen=pg.mkPen(VNAV_COLORS[k]), name=k) for k in ("accel_x","accel_y","accel_z")}
        self.gyr = {k: self.p_gyr.plot([], [], pen=pg.mkPen(VNAV_COLORS[k]), name=k) for k in ("gyro_x","gyro_y","gyro_z")}

        self.p_psd = pg.PlotWidget()
        self.p_psd.setLogMode(x=True, y=True)
        self.p_psd.addLegend()
        right.addWidget(self.p_psd)
        self.psd = {k: self.p_psd.plot([], [], pen=pg.mkPen(VNAV_COLORS[k]), name=k) for k in (
            "mag_x","mag_y","mag_z","accel_x","accel_y","accel_z","gyro_x","gyro_y","gyro_z"
        )}

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plots)
        self.timer.start(150)

    def update_plots(self):
        t, mx = self.b.mag_x.snapshot_last_seconds(self.b.t, self.window_seconds)
        if t.size < 32:
            return
        tr = t - t[0]
        _, my = self.b.mag_y.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, mz = self.b.mag_z.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, ax = self.b.accel_x.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, ay = self.b.accel_y.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, az = self.b.accel_z.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, gx = self.b.gyro_x.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, gy = self.b.gyro_y.snapshot_last_seconds(self.b.t, self.window_seconds)
        _, gz = self.b.gyro_z.snapshot_last_seconds(self.b.t, self.window_seconds)

        self.mag["mag_x"].setData(tr, mx); self.mag["mag_y"].setData(tr, my); self.mag["mag_z"].setData(tr, mz)
        self.acc["accel_x"].setData(tr, ax); self.acc["accel_y"].setData(tr, ay); self.acc["accel_z"].setData(tr, az)
        self.gyr["gyro_x"].setData(tr, gx); self.gyr["gyro_y"].setData(tr, gy); self.gyr["gyro_z"].setData(tr, gz)

        for k, arr in {
            "mag_x": mx, "mag_y": my, "mag_z": mz,
            "accel_x": ax, "accel_y": ay, "accel_z": az,
            "gyro_x": gx, "gyro_y": gy, "gyro_z": gz,
        }.items():
            f, p = psd_fft(t, arr.astype(np.float64))
            if f is not None:
                self.psd[k].setData(f, p)


# -----------------------------
# Window manager
# -----------------------------
class WindowManager(QtCore.QObject):
    def __init__(self, buffers: Dict[str, object], window_seconds: float, epsi_scan_length: int = 1024, parent=None):
        super().__init__(parent)
        self.buffers = buffers
        self.window_seconds = float(window_seconds)
        self.epsi_scan_length = int(epsi_scan_length)

        self.w_efe4: Optional[EFE4Window] = None
        self.w_ttv: Optional[TTVWindow] = None
        self.w_vnav: Optional[VNAVWindow] = None

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.check_open)
        self.timer.start(500)

    @QtCore.pyqtSlot()
    def check_open(self):
        if self.w_efe4 is None:
            b: EFE4Buffers = self.buffers["EFE4"]  # type: ignore
            if b.t.size() > 10:
                self.w_efe4 = EFE4Window(b, self.epsi_scan_length)
                self.w_efe4.show()

        if self.w_ttv is None:
            have = any(self.buffers[t].t.size() > 10 for t in ("TTV1", "TTV2", "TTV3"))  # type: ignore
            if have:
                ttv_dict = {t: self.buffers[t] for t in ("TTV1", "TTV2", "TTV3")}  # type: ignore
                self.w_ttv = TTVWindow(ttv_dict, self.window_seconds)
                self.w_ttv.show()

        if self.w_vnav is None:
            b: VNAVBuffers = self.buffers["VNAV"]  # type: ignore
            if b.t.size() > 10:
                self.w_vnav = VNAVWindow(b, self.window_seconds)
                self.w_vnav.show()


# -----------------------------
# Main
# -----------------------------
def _parse_host_port(two: Tuple[str, str]) -> Tuple[str, int]:
    return two[0], int(two[1])


def main():
    ap = argparse.ArgumentParser(description="Fast MOD-SOM parser + plotting with ring buffers (file/serial/TCP/watch-dir)")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--file", "-f", help="Path to input file")
    group.add_argument("--serial", "-s", help="Serial port (e.g. /dev/tty.usbserial)")
    group.add_argument("--tcp", nargs=2, metavar=("HOST", "PORT"), help="Connect as TCP client to HOST PORT")
    group.add_argument("--tcp-listen", type=int, metavar="PORT", help="Listen as TCP server on PORT (accept 1 client)")
    group.add_argument("--watch-dir", help="Directory another process is writing raw files into; "
                                            "always tails whichever file is currently newest")

    ap.add_argument("--baud", "-b", type=int, default=115200, help="Serial baudrate")
    ap.add_argument("--window", type=float, default=5.0, help="Plot window seconds for TTV/VNAV (default 5); EFE4 uses --epsi-scan-length")
    ap.add_argument("--buffer-seconds", type=float, default=15.0, help="Ring buffer duration (seconds)")
    ap.add_argument("--efe4-fs", type=float, default=320.0, help="EFE4 nominal sample rate for buffer sizing")
    ap.add_argument("--ttv-fs", type=float, default=50.0, help="TTV nominal sample rate for buffer sizing")
    ap.add_argument("--vnav-fs", type=float, default=50.0, help="VNAV nominal sample rate for buffer sizing")
    ap.add_argument("--tcp-bind", default="0.0.0.0", help="Bind address for --tcp-listen (default 0.0.0.0)")
    ap.add_argument("--pattern", default="*.modraw", help="--watch-dir: glob pattern for candidate files (default *.modraw)")
    ap.add_argument("--poll-interval", type=float, default=0.2,
                     help="--watch-dir: seconds between checks for new bytes/new files (default 0.2)")
    ap.add_argument("--seed-bytes", type=int, default=2_000_000,
                     help="--watch-dir: bytes to seed from the end of the newest file on startup, "
                          "so plots aren't blank while waiting for new data (default ~2MB, 0 to disable)")
    ap.add_argument("--epsi-scan-length", type=int, default=1024,
                     help="Sample count used for the EFE4 (epsi) PSD window (default 1024)")
    args = ap.parse_args()

    # EFE4 window shows 3 scan lengths (newest one used for the PSD), so the buffer must hold them.
    cap_efe4 = max(1024, int(np.ceil(args.buffer_seconds * args.efe4_fs)), 3 * args.epsi_scan_length)
    cap_ttv = max(1024, int(np.ceil(args.buffer_seconds * args.ttv_fs)))
    cap_vnav = max(1024, int(np.ceil(args.buffer_seconds * args.vnav_fs)))

    buffers: Dict[str, object] = {
        "EFE4": EFE4Buffers(cap_efe4),
        "TTV1": TTVBuffers(cap_ttv),
        "TTV2": TTVBuffers(cap_ttv),
        "TTV3": TTVBuffers(cap_ttv),
        "VNAV": VNAVBuffers(cap_vnav),
    }

    app = QtWidgets.QApplication([])

    byte_queue: queue.Queue = queue.Queue()
    record_queue: queue.Queue = queue.Queue()
    stop_event = threading.Event()

    parser = EpsiStateMachineParser(record_queue)
    parser_thread = ParserThread(byte_queue, parser, record_queue, stop_event)
    parser_thread.start()

    record_processor = RecordProcessorThread(record_queue, stop_event, buffers)
    record_processor.start()

    # Open chosen source
    source_thread: Optional[ByteSourceThread] = None
    file_obj = None
    ser_obj = None
    sock_src: Optional[SocketSource] = None
    tcp_listen_thread: Optional[TCPListenThread] = None
    dirwatch_src: Optional[DirWatchSource] = None
    using_serial = False

    if args.file:
        file_obj = open(args.file, "rb")
        source_thread = ByteSourceThread(
            file_obj,
            is_stream=False,
            out_queue=byte_queue,
            stop_event=stop_event,
            start_command=None,
        )
        source_thread.start()

    elif args.serial:
        if serial is None:
            raise RuntimeError("pyserial not installed. pip install pyserial")
        using_serial = True
        ser_obj = serial.Serial(port=args.serial, baudrate=args.baud, timeout=0.1)
        source_thread = ByteSourceThread(
            ser_obj,
            is_stream=True,
            out_queue=byte_queue,
            stop_event=stop_event,
            start_command=b"som.start\r\n",
        )
        source_thread.start()

    elif args.tcp:
        host, port = _parse_host_port((args.tcp[0], args.tcp[1]))
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(5.0)
        print(f"[TCP] Connecting to {host}:{port} ...")
        s.connect((host, port))
        s.settimeout(0.5)
        s.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        sock_src = SocketSource(s)
        source_thread = ByteSourceThread(
            sock_src,
            is_stream=True,
            out_queue=byte_queue,
            stop_event=stop_event,
            start_command=None,
        )
        source_thread.start()

    elif args.tcp_listen is not None:
        tcp_listen_thread = TCPListenThread(args.tcp_bind, args.tcp_listen, stop_event)
        tcp_listen_thread.start()
        client_sock = tcp_listen_thread.wait_client(timeout=None)
        if client_sock is None:
            raise RuntimeError("TCP listen ended without a client.")
        client_sock.settimeout(0.5)
        sock_src = SocketSource(client_sock)
        source_thread = ByteSourceThread(
            sock_src,
            is_stream=True,
            out_queue=byte_queue,
            stop_event=stop_event,
            start_command=None,
        )
        source_thread.start()

    elif args.watch_dir:
        dirwatch_src = DirWatchSource(
            args.watch_dir, args.pattern, args.poll_interval, stop_event, args.seed_bytes,
        )
        source_thread = ByteSourceThread(
            dirwatch_src,
            is_stream=True,
            out_queue=byte_queue,
            stop_event=stop_event,
            start_command=None,
        )
        source_thread.start()

    manager = WindowManager(buffers, window_seconds=args.window, epsi_scan_length=args.epsi_scan_length)

    def on_quit():
        stop_event.set()

        # serial stop command only
        if using_serial and ser_obj is not None and getattr(ser_obj, "is_open", False):
            try:
                print("[Main] Sending som.stop to serial...")
                ser_obj.write(b"som.stop\r\n")
                ser_obj.flush()
            except Exception as e:
                print(f"[Main] som.stop failed: {e}")

        # unblock threads
        byte_queue.put(None)

        # join & close
        try:
            if source_thread is not None:
                source_thread.join(timeout=1.0)
        except Exception:
            pass
        try:
            parser_thread.join(timeout=1.0)
        except Exception:
            pass
        try:
            record_queue.put(None)
            record_processor.join(timeout=1.0)
        except Exception:
            pass

        if file_obj is not None:
            try:
                file_obj.close()
            except Exception:
                pass

        if sock_src is not None:
            try:
                sock_src.close()
            except Exception:
                pass

        if ser_obj is not None:
            try:
                ser_obj.close()
            except Exception:
                pass

        if dirwatch_src is not None:
            dirwatch_src.close()

    app.aboutToQuit.connect(on_quit)
    app.exec_()


if __name__ == "__main__":
    main()
