#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J masurca_pe             # job name
#BSUB -n 10                     # assigns 10 cores for execution
#BSUB -R "span[ptile=10]"       # assigns 10 cores per node
#BSUB -R "rusage[mem=1000]"     # reserves 1000MB memory per core
#BSUB -M 1000                   # sets to 1000MB (~1GB) the per process enforceable memory limit.
#BSUB -W 2:00                   # sets to 2 hours the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid

module load MaSuRCA/3.1.2

<<README
    - MaSuRCA manual: http://spades.bioinf.spbau.ru/release3.5.0/manual.html
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
config_file='build_masurca_pe.conf'

pe1_1='../../../data/sra/m_tuberculosis/ERR551611_pe_1.fastq.gz'
pe1_2='../../../data/sra/m_tuberculosis/ERR551611_pe_2.fastq.gz'

######## PARAMETERS ########
threads=8                       # make sure this is <= your BSUB -n value

echo "DATA
PE= pe 400 20 $pe1_1 $pe1_2
END

PARAMETERS
GRAPH_KMER_SIZE=auto
USE_LINKING_MATES=1
LIMIT_JUMP_COVERAGE=60
CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.25 ovlMemory=4GB
NUM_THREADS=$threads
# JF_SIZE is 10 * genome_size
JF_SIZE=44000000
DO_HOMOPOLYMER_TRIM=0
END" > $config_file

########## OUTPUTS #########
# use masurca output defaults

################################### COMMANDS ###################################
# command to run with parameters specified in config file
masurca $config_file

./assemble.sh

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - MaSuRCA:
        Zimin, A. et al. The MaSuRCA genome Assembler. Bioinformatics (2013). doi:10.1093/bioinformatics/btt476 
CITATION

<<NOTES
    Estimated run time on E. coli sample data: ~55 minutes; max memory >10 Gb
        genome size 4.4 Mb
        210,924 300bp paired end read pairs 
NOTES
