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
modsom_generator.py: A data generation script that can mimic specific profiles (e.g., P/T/S samples following a 50-150 sinusoidal pattern) to test the plotting capabilities without needing a live sensor connection.

APEX_MODSOM_liveplot.py: A specialized version of the plotter tailored specifically for APEX real-time data integration.

Faster_app.py: An optimized plotting application focused on speed and reduced rendering latency.

🚀 Requirements
To run this application, you will need Python installed along with the following primary libraries:

numpy

PyQt (or PyQt5/PyQt6 depending on your specific environment)

You can typically install the required dependencies via pip:


pip install numpy PyQt5

💻 Usage
To launch one of the live plotting interfaces, simply run the desired Python script from your terminal:

python MODSOM_liveplot_parallel.py --folder /data/path/.

python MODSOM_liveplot.py --file /data/path/filename (it will plot the data of that file)

For real time 
python MODSOM_liveplot.py --serial serial_port_name -baud 230400 

(If you are testing locally without a live hardware feed, you may need to run modsom_generator.py alongside the plotter to generate the simulated data stream).



