#!/bin/bash
#SBATCH --job-name=alphafold3
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --partition=gpu
#SBATCH --mem=160G
#SBATCH --output=stdout.%x.%j
#SBATCH --error=stderr.%x.%j

module purge

<<SETUP
    You will need to prepare the working directory using the following steps:

    1. mkdir code_3.0.1
    2. cd code_3.0.1
    3. git clone https://github.com/google-deepmind/alphafold3.git
    4. cd alphafold3
    5. git checkout v3.0.1
    6. cd ../..
    7. mkdir input output models
    8. move your af3.bin file into the models directory
    9. configure your sequence in the file: input/alphafold_input.json
   10. Your starting directory structure will look like the following:
       .
       ├── code_3.0.1
       │   └── alphafold3
       ├── inputs
       │   └── alphafold_input.json
       ├── outputs
       ├── run_alphafold_3.0.1_a100_1job_grace.sh
       └── models
           └── af3.bin
SETUP

################################### VARIABLES ##################################
########## INPUTS ##########
# TODO Edit these variables as needed:
export AF3_INPUT_DIR=/sw/hprc/sw/bio/containers/alphafold/examples
export AF3_INPUT_FILE=alphafold_input.json

######## PARAMETERS ########
num_recycles=3
max_template_date='2025-01-01'
export AF3_CODE_DIR=$PWD/code_3.0.1
export AF3_MODEL_PARAMETERS_DIR=$PWD/models
export AF3_DATABASES_DIR=/scratch/data/bio/alphafold3/2025.03.13/
export AF3_IMAGE=/sw/hprc/sw/bio/containers/alphafold/alphafold_3.0.1_unlhcc.sif

########## OUTPUTS #########
export AF3_OUTPUT_DIR=$PWD/output_$SLURM_JOB_ID

################################### COMMANDS ###################################
mkdir -p $AF3_OUTPUT_DIR
singularity exec \
     --nv \
     --bind $AF3_INPUT_DIR:/root/af_input \
     --bind $AF3_OUTPUT_DIR:/root/af_output \
     --bind $AF3_MODEL_PARAMETERS_DIR:/root/models \
     --bind $AF3_DATABASES_DIR:/root/public_databases \
     $AF3_IMAGE \
     python ${AF3_CODE_DIR}/alphafold3/run_alphafold.py \
     --json_path=/root/af_input/$AF3_INPUT_FILE \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --max_template_date=$max_template_date \
     --num_recycles=$num_recycles \
     --output_dir=/root/af_output

################################################################################
<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - AlphaFold3:
        Abramson, J., Adler, J., Dunger, J. et al. Accurate structure prediction of
        biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024).
        https://doi.org/10.1038/s41586-024-07487-w
CITATIONS

<<EXPECTED_OUTPUT
https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md

    output/2pv7/
    ├── 2pv7_confidences.json
    ├── 2pv7_data.json
    ├── 2pv7_model.cif
    ├── 2pv7_summary_confidences.json
    ├── ranking_scores.csv
    ├── seed-1_sample-0
    │   ├── confidences.json
    │   ├── model.cif
    │   └── summary_confidences.json
    ├── seed-1_sample-1
    │   ├── confidences.json
    │   ├── model.cif
    │   └── summary_confidences.json
    ├── seed-1_sample-2
    │   ├── confidences.json
    │   ├── model.cif
    │   └── summary_confidences.json
    ├── seed-1_sample-3
    │   ├── confidences.json
    │   ├── model.cif
    │   └── summary_confidences.json
    ├── seed-1_sample-4
    │   ├── confidences.json
    │   ├── model.cif
    │   └── summary_confidences.json
    └── TERMS_OF_USE.md
EXPECTED_OUTPUT
