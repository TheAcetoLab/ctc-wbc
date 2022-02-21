# Neutrophils escort circulating tumour cells to enable cell cycle progression

## Overview
This repository contains the scripts and necessary metadata for the data analysis and figures appearing in the paper **Neutrophils escort circulating tumour cells to enable cell cycle progression** by [Szczerba BM et al. (2019)](https://www.nature.com/articles/s41586-019-0915-y). The analysis scripts are in R Markdown format and they relay in a `SingleCellExperiment` objects that have been deposited on Gene Expression Omnibus under the accession number [GSE109761](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109761).


## Obtaining the data
The code below shows how to download the data and prepare the files in order to run the project scripts.

1.  Clone this repository
```{bash}
git clone https://github.com/TheAcetoLab/ctc-wbc.git
```
2. Download `SingleCellExperiment` deposited on GEO
This can be done by changing your working directory to the one containing the repository and running `./download_data.sh`

3. Compiling Rmd files
Now you can compile the R Markdown using RStudio or rmarkdown::render(). The list of packages required for each analysis is listed at the top of the script, under the section `General Libraries`. Some scripts need the installation of the R-package [RCA (Reference Component Analysis)](GIS-SP-Group/RCA) hosted in github. There is another package named RCA hosted in CRAN that serves a different purpose. In order to install the RCA package from github you can run the following code :
```{R}
library(devtools)
install_github("GIS-SP-Group/RCA")
```



## Analysis scripts

* [rca-patient.Rmd](https://github.com/CMETlab/ctc-wbc/blob/master/code/Rmd/rca-patient.Rmd) : Cell type identification of white blood cells attached to circulating tumor cells in patients 

* [rca-mouse.Rmd](https://github.com/CMETlab/ctc-wbc/blob/master/code/Rmd/rca-mouse.Rmd) : Cell type identification of white blood cells attached to circulating tumor cells in mice

* [cytokines-analysis.Rmd](https://github.com/CMETlab/ctc-wbc/blob/master/code/Rmd/cytokines-analysis.Rmd) : Characterization of cytokine-mediated crosstalk within CTC-neutrophil clusters

* [cam-analysis.Rmd](https://github.com/CMETlab/ctc-wbc/blob/master/code/Rmd/cam-analysis.Rmd) : Characterization of cell-adhesion molecules (CAMs)-receptor pairs on CTC-neutrophil clusters
