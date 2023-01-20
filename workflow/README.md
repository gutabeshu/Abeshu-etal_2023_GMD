## Xanthos water management
# To reproduce the modeling process
1. Get the required inputs data from https://doi.org/10.5072/zenodo.1075813. Place it in the 'xanthos-wm' folder and unzip.
2. Get xanthos-wm from https://github.com/gutabeshu/xanthos-wm
3. Under xanthos-wm find 'pm_abcd_mrtm_managed.ini' and date the following as needed:
        -'basin_list = 229, 230, 231':provide list of IDs of basins you would like t0 run
        -'set_calibrate=1': this will automatically perform both the first and second stage calibration process
4. The calibration process and model output can be found under 'xanthos-wm/example/output'

# To reproduce the figures
1. Get the model output data can be obtained from https://doi.org/10.5072/zenodo.1075813
2. Open the 'Figure-X.ipyn' you would like to reproduce
3. Update the output data and figure out put directory at the begining
4. Run the script corresponding to the figure of interest.
