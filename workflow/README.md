## Xanthos water management
# To reproduce the modeling process
1. Download the required inputs data from https://doi.org/10.5072/zenodo.1075813
2. Get xanthos-wm from https://github.com/gutabeshu/xanthos-wm
3. Under xanthos-wm find 'pm_abcd_mrtm_managed.ini'. In this file, update the following :
        -provide list of basins you would like t run : 'basin_list'
        -set the 'set_calibrate' to 1, this will automatically perform the first and second stage calibration process

# To reproduce the figures
1. Download model output data can be obtained from https://doi.org/10.5072/zenodo.1075813
2. Open the 'Figure-X.ipyn' you would like to reproduce
3. Update the output data and figure out put directory at the begining
4. Run the script
