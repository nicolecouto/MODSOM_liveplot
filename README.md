MODSOM Liveplot
MODSOM Liveplot is a Python application designed for real-time data visualization. Built with PyQt and NumPy, this suite of tools allows for the continuous ingestion and live plotting of simulated or real-time data streams (including APEX data).

📂 Core Plotting Modules
The repository contains several iterations of the live plotting script to accommodate different performance needs and architectures:

1. MODSOM_liveplot.py
This is the standard, baseline implementation for real-time plotting. It establishes the primary PyQt interface and handles the basic loop for fetching data and updating the graph. It is ideal for straightforward use cases with moderate data rates.

2. MODSOM_liveplot_v2.py
An updated iteration of the base live plotter. This version typically includes structural improvements, UI refinements, or alternative plotting mechanisms compared to the original, offering a more robust approach to handling the live data stream.

3. MODSOM_liveplot_parallel.py
A high-performance version of the application that utilizes parallel processing (such as threading or multiprocessing). By separating the data ingestion/generation from the UI and plotting updates, this script ensures that the GUI remains responsive even when handling high-frequency data or performing heavier computations.

🛠️ Other Included Tools
modsom_generator.py: A data generation script that can mimic specific profiles (
e.g., P/T/S samples following a 50-150 sinusoidal pattern) to test the plotting capabilities without needing a live sensor connection.

APEX_MODSOM_liveplot.py: A specialized version of the plotter tailored specifically for APEX real-time data integration.

Faster_app.py: An optimized plotting application focused on speed and reduced rendering latency.

🚀 Getting started

These steps assume macOS or Linux with Python 3 installed (tested with Python 3.9 and 3.12). Run them from the cloned `MODSOM_liveplot` folder.

1. Create a virtual environment and install the packages (one time)

```bash
cd MODSOM_liveplot
python3 -m venv .venv
source .venv/bin/activate
pip install numpy PyQt5 pyqtgraph
```

Optional extras, depending on what you run:

```bash
pip install pyserial          # only for --serial (reading the instrument's serial port directly)
pip install netCDF4 pandas    # only for --folder (batch conversion of .modraw files to NetCDF)
```

The serial package is `pyserial`, not `serial` - an unrelated package called `serial` exists on PyPI, and installing it by mistake is a common trap.

Versions known to work: Python 3.9 with numpy 2.0.2, PyQt5 5.15.11, pyqtgraph 0.13.7, netCDF4 1.7.2, pandas 2.3.3; and Python 3.12 with numpy 2.5.3, PyQt5 5.15.11, pyqtgraph 0.14.0 (live plot only).

2. Activate the environment (every new terminal)

```bash
cd MODSOM_liveplot
source .venv/bin/activate
```

💻 Usage

Live epsi spectra (Faster_app.py)

Pick one input source:

```bash
# Watch a directory another logger is writing .modraw files into (the usual ship setup)
python Faster_app.py --watch-dir /path/to/raw --epsi-scan-length 1024

# Read the instrument's serial port directly (needs pyserial)
python Faster_app.py --serial /dev/tty.usbserial-XXXX --baud 115200

# Play back a saved file
python Faster_app.py --file /path/to/file.modraw
```

`python Faster_app.py --help` lists every option.

- `--watch-dir` never opens the serial port itself. It tails whichever file in the directory is currently newest and follows along as new files appear, so you can view live data while another process holds the port and saves the raw `.modraw` files.
- `--epsi-scan-length` (default 1024) is the number of samples in each spectrum. The time series on the left show exactly the samples the spectrum is computed from.
- The EFE4 window has a Light mode / Dark mode button at the top right (light mode is better for screenshots, e.g. on Slack), and a checkbox next to each spectrum in the legend to show or hide it.

Batch NetCDF conversion (needs netCDF4 and pandas)

```bash
python MODSOM_liveplot_parallel.py --folder /path/to/EPSI_PROCESSING/raw
```

Converts every `.modraw` file in the folder to NetCDF, plus a combined `combined_1Hz.nc`, written to a `netcdf` directory next to it (here `/path/to/EPSI_PROCESSING/netcdf`).

Other plotters

```bash
python MODSOM_liveplot.py --file /data/path/filename           # plot the data in that file
python MODSOM_liveplot.py --serial serial_port_name --baud 230400   # real time
```

Testing without an instrument

`modsom_generator.py` writes a simulated `.modraw` file you can play back (about 5 MB per second of data):

```bash
python modsom_generator.py --file test.modraw --duration 10
python Faster_app.py --file test.modraw
```
