#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's environment.
#BSUB -J stringtie              # job name
#BSUB -n 20                     # assigns 20 cores for execution
#BSUB -R "span[ptile=20]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2700MB per process enforceable memory limit.
#BSUB -W 24:00                  # sets to 24 hour the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid

module load StringTie/1.3.3-GCCcore-6.3.0

<<README
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
threads=20                       # make sure this is <= your BSUB -n value

# reference annotation_file can be GTF or GFF3 format
annotation_file='/scratch/datasets/GCATemplates/data/rnaseqc/c_parapsilosis/C_parapsilosis_CDC317_version_s01-m03-r11_Mt_features.gtf'

bam_file='/scratch/datasets/GCATemplates/data/rnaseqc/c_parapsilosis/DR18_bowtie.bam'

label='DR18'
out_gtf="${label}_out.gtf"

################################### COMMANDS ###################################
#
stringtie $bam_file -G $annotation_file -l $label -o $out_gtf -p $threads


<<CITATION
    - Acknowledge TAMU HPRC: http://sc.tamu.edu/research/citation.php

    - StringTie:
        Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL Transcript-level expression
        analysis of RNA-seq experiments with HISAT, StringTie and Ballgown,
        Nature Protocols 11, 1650-1667 (2016), doi:10.1038/nprot.2016.095 
CITATION
