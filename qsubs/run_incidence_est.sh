#!/bin/bash
module purge
module load R/4.4.3

Rscript ../code/tcga/incidence_estimation.R -m cnaseg -e BRCA1 -c BRCA