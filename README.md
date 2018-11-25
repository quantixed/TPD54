# TPD54
Code for TPD54-related projects



--

### FRAP analysis

Using `FRAP.ijm`, individual movies opened from an mvd2 Ultraview database using BioFormats, can be analysed. The code will specify the FRAP ROI, the user is prompted to add a bg ROI and an ROI encompassing the whole cell. The data from each ROI is saved as CSV and the ROIs are saved for reproducibility.

In Igor, `FRAPKinetics.ipf` and `ParseTimeStampsFromOME.ipf` are used to process this data. Averages, fits and analysis are all generated, including a figure.

--

### Miscellaneous

`bcaFromCSVUneven.r` can be used in R to calculate estimation statistics. The results are output as csv to use in Igor.

`RUSHKinetics.ipf` is a semi-automated procedure to process data from RUSH experiments.