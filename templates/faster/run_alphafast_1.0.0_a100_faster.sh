#!/bin/bash
#SBATCH --job-name=alphafast
#SBATCH --time=02:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=16          # CPUs (threads) per command
#SBATCH --mem=62G                   # total memory per node
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1           # request one GPU
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

module purge
module --ignore-cache load gompi/2025b AlphaFast/1.0.0-CUDA-12.8.0

<<README
    - AlphaFold manual: https://github.com/RomeroLab/alphafast
README

############ SYNOPSIS #############
# this script will run alphafast
# use just one GPU when running one input sequence
# AlphaFast will use multiple GPUs only when there are many sequences

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
input_dir=$PWD/inputs
mkdir -p $input_dir
# example: copy example .json file to the inputs directory with the following commands
cp -n /scratch/data/bio/alphafast/examples/2PV7/2PV7.json $input_dir

# directory containing your weights file af3.bin.zst
weights_dir=$SCRATCH

######## PARAMETERS ########

########## OUTPUTS #########
output_dir=out_${SLURM_JOB_NAME}_${SLURM_JOB_ID}

################################### COMMANDS ###################################

run_alphafast.sh --input_dir $input_dir \
    --output_dir $output_dir \
    --weights_dir $weights_dir \
    --db_dir /scratch/data/bio/alphafast/db \
    --container /sw/hprc/sw/bio/containers/alphafast/alphafast.sif \
    --num_gpus $SLURM_GPUS_ON_NODE

################################################################################
<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - AlphaFast:
        Benjamin C. Perry, Jeonghyeon Kim, Philip A. Romero
        AlphaFast: High-throughput AlphaFold 3 via GPU-accelerated MSA construction.
        bioRxiv 2026.02.17.706409; doi: https://doi.org/10.64898/2026.02.17.706409 

    - AlphaFold3:
        Abramson, J., Adler, J., Dunger, J. et al. Accurate structure prediction of
        biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024).
        https://doi.org/10.1038/s41586-024-07487-w

    - MMseqs2:
        Kallenborn, F., Chacon, A., Hundt, C. et al. GPU-accelerated homology search
        with MMseqs2. Nat Methods 22, 2024–2027 (2025). https://doi.org/10.1038/s41592-025-02819-8

CITATIONS

