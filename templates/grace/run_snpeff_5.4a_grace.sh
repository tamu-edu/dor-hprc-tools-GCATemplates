#!/bin/bash
#SBATCH --job-name=snpeff           # job name
#SBATCH --time=01:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=7G                    # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

module purge
module load GCCcore/13.3.0 snpEff/5.4a-Java-21

<<README
    - snpEff manual: http://snpeff.sourceforge.net/SnpEff_manual.html
README


################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
vcf_file='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34.vcf'

# to see available snpEff databases: java -jar $EBROOTSNPEFF/snpEff.jar databases > dbs_snpeff.txt
# to download a snpEff database:     java -jar $EBROOTSNPEFF/snpEff.jar download DATABASE_NAME
snpeff_database='Candida_dubliniensis_cd36_gca_000026945'

######## PARAMETERS ########
sample='DR34'

########## OUTPUTS #########
outfile="${sample}_snpeff_ann.vcf"

################################### COMMANDS ###################################
# NOTE for bacterial genomes you may want to adjust the upstream/downstream option (-ud) default is 5000 bp
#      or use the options -no-upstream -no-downstream -no-intergenic
java -Xmx4g -jar $EBROOTSNPEFF/snpEff.jar ann -v $snpeff_database $vcf_file > $outfile

################################################################################
<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - snpEff:
        Cingolani P1, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM.
        A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff:
        SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.
        Fly (Austin). 2012 Apr-Jun;6(2):80-92. doi: 10.4161/fly.19695.
CITATION


