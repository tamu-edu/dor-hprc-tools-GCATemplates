#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J repeatmasker           # job name
#BSUB -n 20                     # assigns 20 cores for execution
#BSUB -R "span[ptile=20]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2700MB process enforceable memory limit. (M * n)
#BSUB -W 96:00                  # sets to 96  hour the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid

module load RepeatScout/1.0.5-GCCcore-6.3.0
module load RepeatMasker/4.0.7-GCCcore-6.3.0

<<README
    RepeatMasker homepage: http://www.repeatmasker.org/
    TRF homepage: https://tandem.bu.edu/trf/trf.html
    RepeatScout manual: http://bix.ucsd.edu/repeatscout/readme.1.0.5.txt
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
input_fasta='c_albicans_SC5314.fasta'

# see the wiki page for details on pre-configured RepeatMasker species
# https://sc.tamu.edu/wiki/index.php/Ada:Bioinformatics#RepeatMasker
repeatmasker_species='fungi'

threads=10                      # use #BSUB -n 20 and set threads to 10 or compute node will be overloaded

################################### COMMANDS ###################################
# RepeatMasker

RepeatMasker -s -species $repeatmasker_species -pa $threads $input_fasta

rm_out="${input_fasta}.masked"

################################################################################
# TRF
trf $rm_out 2 7 7 80 10 50 500 -m

rm_trf_masked_fasta="$rm_out.2.7.7.80.10.50.500.mask"

################################################################################
# RepeatScout
build_lmer_table -sequence $rm_trf_masked_fasta -freq output_lmer.frequency

RepeatScout -sequence $rm_trf_masked_fasta -output out_repeatscout_C_tropicalis_MYA-3404_chromosomes.fasta -freq output_lmer.frequency -vv

filter-stage-1.prl out_repeatscout_C_tropicalis_MYA-3404_chromosomes.fasta > output_repeats.fas.filtered_1

RepeatMasker $rm_trf_masked_fasta -lib output_repeats.fas.filtered_1

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - RepeatMasker: Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0.

    - TRF: G. Benson, "Tandem repeats finder: a program to analyze DNA sequences"
         Nucleic Acids Research (1999) Vol. 27, No. 2, pp. 573-580.

    - RepeatScout: Price A.L., Jones N.C. and Pevzner P.A. 2005.  De novo identification of 
        repeat families in large genomes.  To appear in Proceedings of the
        13 Annual International conference on Intelligent Systems for
        Molecular Biology (ISMB-05).  Detroit, Michigan.
CITATION
