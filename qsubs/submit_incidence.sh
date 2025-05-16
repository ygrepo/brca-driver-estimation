#!/bin/bash
# submit_incidence.sh

# — adjust these if you need to add or remove combos —
params=(
  "BRCA BRCA1 cna"
  "BRCA BRCA1 cnaseg"
  "BRCA BRCA2 cna"
  "BRCA BRCA2 cnaseg"
  "OV   BRCA1 cna"
  "OV   BRCA1 cnaseg"
  "OV   BRCA2 cna"
  "OV   BRCA2 cnaseg"
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
    -R "rusage[mem=10000]" \
    -W 1:00 \
    -o "logs/${jobname}.%J.out" \
    -eo "logs/${jobname}.%J.err" \
    "module purge; module load R/4.4.3; \
     Rscript /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/brca-driver-estimation/code/tcga/incidence_estimation.R \
       -m '$mut' -e '$gene' -c '$cancer'"
done
#!/bin/bash
# submit_incidence.sh

# — adjust these if you need to add or remove combos —
params=(
  "BRCA BRCA1 cna"
  "BRCA BRCA1 cnaseg"
  "BRCA BRCA2 cna"
  "BRCA BRCA2 cnaseg"
  "OV   BRCA1 cna"
  "OV   BRCA1 cnaseg"
  "OV   BRCA2 cna"
  "OV   BRCA2 cnaseg"
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
    -R "rusage[mem=10000]" \
    -W 1:00 \
    -o "logs/${jobname}.%J.out" \
    -eo "logs/${jobname}.%J.err" \
    "module purge; module load R/4.4.3; \
     Rscript /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/brca-driver-estimation/code/tcga/incidence_estimation.R \
       -m '$mut' -e '$gene' -c '$cancer'"
done
