**NRVM-DOX-proteomics**
<br>Proteomic response of male and female NRVM exposed to doxorubicin

**Description**
<br>This repository contains code for differential protein abundance analysis using the limma R package. The available dataset was obtained from primary cardiomyocytes isolated from individual hearts of postnatal day 2 rat pups, allowing paired (i.e. within-animal) analysis of doxorubicin (DOX) and vehicle treated cells. The provided script enables quality control assessment and statistical analysis. Log2 transformed, quantile normalised intensities are analysed using a linear model to include treatment, with sex and batch as covariates. Bioreplicate blocking is applied to estimate correlation between paired samples. Moderated t-statistics are computed by linear model fitting, group contrasts and empirical Bayes moderation. 

Code for the manuscript: Sanela Dozic, Michael G. Leeming, Keshava K. Datta, Hitesh Kore, 
Erin J. Howden, Benjamin L. Parker, Lea M.D. Delbridge, Kate L. Weeks, “Sex differences in cardiomyocyte proteomic responses to the cardiotoxic chemotherapy drug doxorubicin.”

**Downloading the raw data**
<br>Mass spectrometry proteomics data can be downloaded from the ProteomeXchange Consortium (http://proteomecentral.proteomexchange.org) via the PRIDE partner repository with the dataset identifier PXD076454. The raw data files carry the project identifier ‘250902_SD-11625’ tagged with the respective sample number (e.g. _01). Sample details are provided in the ‘Experimental_Annotation’ file. Raw data files were processed with FragPipe and the resulting ‘combined_protein’ file used for limma analysis. 

**Usage**
<br>The ‘combined_protein_n22’ and ‘Experimental_Annotation’ data files provided in this repository were used as input for differential protein abundance analysis. The provided script is presented in a form that can be modified to apply different thresholds and imputations for exploratory data analysis. The following libraries are included in the script and installed to run the analysis:
•	  limma
•	  ggplot2
•	  ggrepel
•	  pcaMethods
•	  tidyr
•	  dplyr
•	  tibble
