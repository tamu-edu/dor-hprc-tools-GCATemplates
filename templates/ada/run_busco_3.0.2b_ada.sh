#BSUB -L /bin/bash              # uses the bash login shell to initialize the job's execution environment.
#BSUB -J busco                  # job name
#BSUB -n 20                     # assigns 20 cores for execution
#BSUB -R "span[ptile=20]"       # assigns 20 cores per node
#BSUB -R "rusage[mem=2700]"     # reserves 2700MB memory per core
#BSUB -M 2700                   # sets to 2700MB per process enforceable memory limit. Total memory = (M * n)
#BSUB -W 24:00                  # sets to 24 hours the job's runtime wall-clock limit.
#BSUB -o stdout.%J              # directs the job's standard output to stdout.jobid
#BSUB -e stderr.%J              # directs the job's standard error to stderr.jobid

module load BUSCO/3.0.2b-intel-2015B-python-2.7.6-generic

<<README
    - BUSCO manual: http://buscos.ezlab.org/files/BUSCO_userguide.pdf
    - Homepage: http://busco.ezlab.org/
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
genome_file='/scratch/datasets/gmod_genomes/C_albicans_SC5314/C_albicans_SC5314_version_A21-s02-m09-r10_chromosomes.fasta'

######## PARAMETERS ########
threads=20                      # make sure this is <= your BSUB -n value

# see available BUSCO lineages in this directory: /scratch/datasets/BUSCO/v3.0.2/
busco_lineage='/scratch/datasets/BUSCO/v3.0.2/fungi_odb9'
busco_mode='genome'             # genome, transcriptome, proteins

# --species:  see available species here: /software/easybuild/software/AUGUSTUS/3.3-intel-2015B/config/species/
augustus_species='candida_albicans'
augustus_model='intronless'     # partial (default), intronless, complete, atleastone, exactlyone

########## OUTPUTS #########
output_prefix="out_busco_${augustus_species}"

################################### COMMANDS ###################################
# copy augustus config to $SCRATCH
if [ ! -d "$SCRATCH/my_augustus_config/config" ]; then
  echo "Copying AUGUSTUS config directories to $SCRATCH/my_augustus_config"
  mkdir $SCRATCH/my_augustus_config
    if [ "-$EBROOTAUGUSTUS" == "-" ]; then
        echo "Augustus module not loaded"; exit 1
    fi
  rsync -r $EBROOTAUGUSTUS/ $SCRATCH/my_augustus_config
  chmod -R 755 $SCRATCH/my_augustus_config
fi
export AUGUSTUS_CONFIG_PATH="$SCRATCH/my_augustus_config/config"

python $EBROOTBUSCO/scripts/run_BUSCO.py -f --in $genome_file --mode $busco_mode --species $augustus_species --cpu $threads --lineage $busco_lineage --out $output_prefix --tmp_path $TMPDIR --augustus_parameters="--genemodel=$augustus_model"

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - BUSCO:
        BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.
        Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
        Bioinformatics, published online June 9, 2015 | doi: 10.1093/bioinformatics/btv351
CITATION
