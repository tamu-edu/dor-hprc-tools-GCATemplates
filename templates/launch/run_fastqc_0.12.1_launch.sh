#!/bin/bash                                                        
#SBATCH --job-name=fastqc           # job name       
#SBATCH --time=01:00:00             # max job run time dd-hh:mm:ss
#SBATCH --nodes=1                   # number of compute nodes to request
#SBATCH --ntasks-per-node=2         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command  
#SBATCH --mem=6G                    # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to a file named stdout.jobname.jobid
#SBATCH --error=stderr.%x.%j        # save stderr to a file named stderr.jobname.jobid

module load FastQC/0.12.1-Java-11

<<README
    - FASTQC homepage: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
    - FASTQC manual: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_R1='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R1.fastq.gz'
pe1_R2='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R2.fastq.gz'

######## PARAMETERS ########
threads=$SLURM_NTASKS_PER_NODE

########## OUTPUTS #########
output_dir='./'

################################### COMMANDS ###################################
# use -o <directory> to save results to <directory> instead of directory where reads are located
#   <directory> must already exist before using -o <directory> option
# --nogroup will calculate average at each base instead of bins after the first 50 bp
# fastqc runs one thread per file; using 20 threads for 2 files does not speed up the processing

fastqc -t $threads -o $output_dir $pe1_R1 $pe1_R2

################################################################################
<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
CITATION
