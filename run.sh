#!/bin/bash

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="/work/calonga/miniconda3/bin:$PATH"
if [ ! -d /work/calonga/miniconda3/envs/workflow_2020_population_genetics ]; then
    conda env create -n workflow_2020_population_genetics --file environment.yaml
fi
source activate workflow_2020_population_genetics

echo "RUNNING SNAKEMAKE WORKFLOW..."

#snakemake --cluster "sbatch -c 8 --mem-per-cpu=30GB" -j 10 -p --use-conda LAZARIDIS_NEAREAST/HumanOriginsPublic2068.vcf

snakemake --rerun-incomplete -p --use-conda LAZARIDIS_NEAREAST/HumanOriginsPublic2068.vcf


