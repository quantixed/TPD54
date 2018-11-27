# TPD54
Code for TPD54-related projects

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1585668.svg)](https://doi.org/10.5281/zenodo.1585668)

--

### Volcano plot

Proteomic data is available in `data/MSData.tsv`. It can be loaded into Igor and the volcano plot and analysis generated using `VolcanoPlot.ipf`. The output from the paper is available in `output/VolcanoPlot.tsv`, note that the imputation of missing values means that the output of `VolcanoPlot.ipf` will be slightly different each time. 

--

### Vesicle capture

Analysis of mitochondria and vesicle profiles is done using `IMODModelTPD.ipf`. The outputs from IMOD model2point are in `data/EM`. These are IMOD models converted to text files that Igor can read. Load `FileName++.txt` and `TiffName++.txt` and then run IMOD Model Analysis from the macros menu. Point Igor to the dataset you'd like to process. Note that some datasets have thousands of vesicles and the proximity analysis can take a long time to run.

--

### FRAP analysis

Using `FRAP.ijm`, individual movies opened from an mvd2 Ultraview database using BioFormats, can be analysed. The code will specify the FRAP ROI, the user is prompted to add a bg ROI and an ROI encompassing the whole cell. The data from each ROI is saved as CSV and the ROIs are saved for reproducibility.

In Igor, `FRAPKinetics.ipf` and `ParseTimeStampsFromOME.ipf` are used to process this data. Averages, fits and analysis are all generated, including a figure.

--

### Miscellaneous

`bcaFromCSVUneven.r` can be used in R to calculate estimation statistics. The results are output as csv to use in Igor.

`RUSHKinetics.ipf` is a semi-automated procedure to process data from RUSH experiments.
