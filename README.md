# Tissue-specific transcriptional imprinting and heterogeneity in human innate lymphoid cells revealed by full-length single-cell RNA-sequencing
Repository contacting the code related to the article "Tissue-specific transcriptional imprinting and heterogeneity in human innate lymphoid cells revealed by full-length single-cell RNA-sequencing" [Mazzurana, Czarnewski et al (2020) *Cell Research*](https://doi.org/10.1038/s41422-020-00445-x)


The analysis done herein can be reproduced by installing `conda` and running the `run_workflow.sh` script. It will download the dataset from GEO (GSE150050), the scripts from Sauron.



1. Clone this repository
```
git clone https://github.com/czarnewski/ilc_tissue_imprinting.git
```


2. Create and activate the conda environment
```
cd ilc_tissue_imprinting

conda activate base
conda install -c conda-forge mamba

mamba env create -n ilc_tissue_imprinting -f sauron_environment_20201209.yml
conda activate ilc_tissue_imprinting
```


3. Run the analysis workflow
```
sh run_workflow.sh
```
