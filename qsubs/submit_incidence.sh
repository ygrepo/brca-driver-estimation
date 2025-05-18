#!/bin/bash
# submit_incidence.sh

# — adjust these if you need to add or remove combos —
# params=(
#   "BRCA BRCA1 snv"
#   "BRCA BRCA2 snv"
#   "OV BRCA1 snv"
#   "OV BRCA2 snv"
#   "BRCA BRCA1 cna"
#   "BRCA BRCA1 cnaseg"
#   "BRCA BRCA2 cna"
#   "BRCA BRCA2 cnaseg"
#   "BRCA BRCA1 amplification"
#   "BRCA BRCA1 amplificationseg"
#   "BRCA BRCA2 amplification"
#   "BRCA BRCA2 amplificationseg"
#   "OV   BRCA1 cna"
#   "OV   BRCA1 cnaseg"
#   "OV   BRCA1 amplification"
#   "OV   BRCA1 amplificationseg"
#   "OV   BRCA2 cna"
#   "OV   BRCA2 cnaseg"
#   "OV   BRCA2 amplification"
#   "OV   BRCA2 amplificationseg"
# )

params=(
  "BRCA BRCA1 deletion"
  "BRCA BRCA2 deletion"
  "OV BRCA1 deletion"
  "OV BRCA2 deletion"
  "BRCA BRCA1 deletionseg"
  "BRCA BRCA2 deletionseg"
  "OV BRCA1 deletionseg"
  "OV BRCA2 deletionseg"
  )
  
# make sure log dir exists
mkdir -p logs

# load modules once (each bsub will re‐load inside its shell)
module purge
module load R/4.4.3

for entry in "${params[@]}"; do
  # split the triple into three variables
  read -r cancer gene mut <<< "$entry"
  jobname="inc_${cancer}_${gene}_${mut}"

  bsub \
    -J "$jobname" \
    -P acc_DiseaseGeneCell \
    -q premium \
    -n 1 \
    -R "rusage[mem=30000]" \
    -W 1:00 \
    -o "logs/${jobname}.%J.out" \
    -eo "logs/${jobname}.%J.err" \
    "module purge; module load R/4.4.3; \
     Rscript /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/brca-driver-estimation/code/tcga/incidence_estimation.R \
       -m '$mut' -e '$gene' -c '$cancer'"
done
