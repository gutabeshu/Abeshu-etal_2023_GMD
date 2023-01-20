_your zenodo badge here_?

# Abeshu et al., 2023, GMD
**A new representation of water management for global hydrological models**

Guta Wakbulcho Abeshu<sup>1</sup>, Fuqiang Tian<sup>2</sup>, Hongchang Hu<sup>2</sup>, Yuan Zhuang<sup>2</sup>,
Mohamad Hejazi<sup>3</sup>, Sean Turner<sup>3</sup>, Thomas Wild<sup>3</sup>, Mengqi Zhao<sup>4</sup>,
A F M Kamal Chowdhury<sup>5</sup>, Chris Vernon<sup>3</sup>, and Hong-Yi Li<sup>1\*</sup>

<sup>1 </sup> Department of Civil and Environmental Engineering, University of Houston, Houston, TX, 77204, US
<sup>2 </sup> Department of Hydraulic Engineering, Tsinghua University, Beijing, China
<sup>3 </sup> Joint Global Change Research Institute, Pacific Northwest National Laboratory, College Park, MD 20740, US
<sup>4 </sup> Pacific Northwest National Laboratory, Richland, WA 99354, US
<sup>5 </sup> Earth System Science Interdisciplinary Center (ESSIC), University of Maryland, College Park, MD 20740, US

\* corresponding author:  hongyili.jadison@gmail.com

## Abstract
This study enhances an existing global hydrology model, Xanthos, by adding a new water management module, treating irrigation, hydropower, and flood-control reservoirs differently. We determine a unique operation rule for each of 3790 large reservoirs globally based on their primary purposes, i.e., hydropower, irrigation, flood-control or others. We apply the enhanced Xanthos globally at a 0.5-degree spatial resolution. The model performance is satisfactory in reproducing observed stream flow variability under normal and water scarcity conditions at 91 large river basins with good data availability. Hydropower reservoirs' simulated storage and release patterns are quite different from flood-control reservoirs, suggesting significant implications for freshwater resource management at the regional or larger scales. This new global water management modeling framework will allow for the assessment of future reservoir development and management from a coupled human-natural system perspective.

## Journal reference
Abeshu, G.W., Fuqiang Tian, F., Hongchang Hu, H., Yuan Zhuang, Y., Hejazi, M., Sean Turner, S., Thomas Wild, T., Mengqi Zhao, M., Chowdhury, K., Vernon, C., and Li, H., A new water management module for global hydrologic models. GMD

## Code reference
https://github.com/gutabeshu/xanthos-wm

## Data reference

### Input data
The following are input data required to run the model

- **Forcing Data:**
In the study, we used gridded climatic data from the WATer and global CHange (WATCH; Weedon et al., 2011) dataset from 1971 to 2001. The global monthly datasets include precipitation, maximum temperature, relative humidity and minimum temperature. DataHub: https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/3/

- **Water demand and water use data:**
Monthly water demand and consumptive water use data for various sectors at a 0.5-degree resolution are from Huang et al. (2018), which are available from 1971 to 2010. DataHub: https://doi.org/10.5281/zenodo.1209296

- **Reservoirs Data:**
Global reservoir data are obtained from the GRanD dataset (Lehner et al., 2011). The GRanD database we use here only considers reservoirs with storage capacity values greater than 0.1 km3. Reservoirs with missing storage capacity and those identified with purposes such as tide control are dropped, reducing the total GranD reservoirs from 6862 to 6847. DataHub: https://doi.org/10.5072/zenodo.1075064

- **Streamflow Data:**
Observed stream flow data for model parameter identification and validation are obtained from the Global Runoff Data Center (GRDC)DataHub: https://www.bafg.de/GRDC

### Output data
Abeshu, G.W. (2022). Xanthos-WM output dataset name. DataHub: https://doi.org/10.5072/zenodo.1075879

## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| Xanthos | version |https://github.com/gutabeshu/xanthos_managed steps for installation can be found at  https://github.com/JGCRI/xanthos | https://doi.10.5281/zenodo.5177210  |


## Reproduce my experiment
Fill in detailed info here or link to other documentation that is a thorough walkthrough of how to use what is in this repository to reproduce your experiment.


1. Install the software components required to conduct the experiment from [Contributing modeling software](#contributing-modeling-software)
2. Download and install the supporting input data required to conduct the experiment from [Input data](https://doi.org/10.5072/zenodo.1075813)
3. Run the following scripts in the `workflow` directory to re-create this experiment:

| Experiment Name | Description | How to Run |
| --- | --- | --- |
| `Runoff` | Script to run the first part of my experiment | under shell_script_to_run  **sbatch run_for_91basins.sh** edit directory to **runoff-watch-setup** in 'workflow'|
| `Streamflow` | Script to run the last part of my experiment | under shell_script_to_run  **sbatch run_for_91basins.sh** edit directory to **flow-watch-setup** in 'workflow'|

4. Download and unzip the output data from my experiment [Output data](https://doi.org/10.5072/zenodo.1075879)
5. Run notebook scripts in the `figure` directory to compare my outputs to those from the publication

| Script Name | Description | How to Run |
| ---  | --- | --- |
| `Figure_X.py` | Script to compare my outputs to the original | edit the `dir_in` at the beginning to your output directory|

***File Modification***\
Check each file listed in [Table 5](#table5) and modify every directory within those files to the directory that holds your data.
<a name="table5"></a>
**Table 5:** Files to be modified to run the model.

| Model | Programming Language | Files to be Modified |
|---|:-:|---|
| Xanthos | Python 3.3+ | under workflow/flow-watch-setup/xanthos(basin number):  **pm_abcd_mrtm.ini** under workflow/: **Basin (basin number).py** |


## Reproduce my figures
Use the scripts found in the `figures` directory to reproduce the figures used in this publication.

| Script Name | Description | How to Run |
| --- | --- | --- |
| `Figure_4.ipynb`  | Script to generate Figure 5    | `Figure_4.ipynb change input/output directory path`  & run on notebook  |
| `Figure_5.ipynb`  | Script to generate Figure 5    | `Figure_5.ipynb change input/output directory path`  & run on notebook  |
| `Figure_6.ipynb`  | Script to generate Figure 6    | `Figure_6.ipynb change input/output directory path`  & run on notebook  |
| `Figure_7.ipynb`  | Script to generate Figure 7    | `Figure_7.ipynb change input/output directory path`  & run on notebook  |
| `Figure_8.ipynb`  | Script to generate Figure 8    | `Figure_8.ipynb change input/output directory path`  & run on notebook  |
| `Figure_9.ipynb`  | Script to generate Figure 9    | `Figure_9.ipynb change input/output directory path`  & run on notebook   |
| `Figure_10.ipynb` | Script to generate Figure 10   | `Figure_10.ipynb change input/output directory path` & run on notebook  |
| `Figure_11.ipynb` | Script to generate Figure 11   | `Figure_11.ipynb change input/output directory path` & run on notebook  |
| `Figure_12.ipynb` | Script to generate Figure 12   | `Figure_12.ipynb change input/output directory path` & run on notebook  |
| `Figure_13.ipynb` | Script to generate Figure 13   | `Figure_13.ipynb change input/output directory path` & run on notebook  |
| `Figure_14.ipynb` | Script to generate Figure 14   | `Figure_14.ipynb change input/output directory path` & run on notebook  |
| `Figure_15.ipynb` | Script to generate Figure 15   | `Figure_15.ipynb change input/output directory path` & run on notebook  |
