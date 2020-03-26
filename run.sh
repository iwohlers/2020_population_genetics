#!/bin/bash

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="/scratch/programs/miniconda3/bin:$PATH"
if [ ! -d /scratch/programs/miniconda3/envs/workflow_2020_population_genetics ]; then
    conda env create -n workflow_2020_population_genetics --file environment.yaml
fi
source activate workflow_2020_population_genetics

echo "RUNNING SNAKEMAKE WORKFLOW..."

snakemake --unlock --rerun-incomplete -p --use-conda preprocess_egyptgsa preprocess_egyptgsapso preprocess_egyptwgs preprocess_lazaridis preprocess_busby  preprocess_1000g  preprocess_fernandes preprocess_scott preprocess_bergstroem preprocess_rodriguez
#snakemake -p -j 10 --use-conda liftover_egypt_gsa_to_38_all
#--cluster "sbatch -c 8 --mem-per-cpu=30GB --partition=lied"

# If environment.yaml has been changed, the existing environment needs to be removed 
# in order to re-generate the environment using: 
# source ~/.bashrc; conda env remove -n workflow_2020_population_genetics

