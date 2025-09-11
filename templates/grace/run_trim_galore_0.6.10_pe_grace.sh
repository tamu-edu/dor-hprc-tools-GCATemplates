#!/bin/bash
#SBATCH --job-name=trim_galore      # job name
#SBATCH --time=01:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=12G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

module purge
module load GCCcore/12.3.0 Trim_Galore/0.6.10

<<README
    - Trim_Galore manual:
          https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md       
README

######### SYNOPSIS #########
# quality trims paired input files and runs fastqc on trimmed files

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_R1='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R1.fastq.gz'
pe1_R2='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R2.fastq.gz'

######## PARAMETERS ########
cores=$SLURM_CPUS_PER_TASK

########## OUTPUTS #########
out_dir='trimmed_files'

################################### COMMANDS ###################################

trim_galore --cores $cores --fastqc --output_dir $out_dir --paired $pe1_R1 $pe1_R2

################################################################################
<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - Trim_Galore
          Trimgalore (2021), GitHub repository, https://github.com/FelixKrueger/TrimGalore
