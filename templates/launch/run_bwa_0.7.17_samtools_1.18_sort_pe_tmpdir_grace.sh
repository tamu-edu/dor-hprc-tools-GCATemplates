#!/bin/bash
#SBATCH --job-name=bwa_pe           # job name
#SBATCH --time=24:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=360G                  # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

module purge
module load GCC/12.3.0 SAMtools/1.18 BWA/0.7.17

<<README
    - BWA manual: http://bio-bwa.sourceforge.net/bwa.shtml
    - SAMtools manual: http://samtools.github.io/hts-specs/SAMv1.pdf
README

######### SYNOPSIS #########
# aligns to a reference; uses $TMPDIR for sorting; output is a sorted bam file

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_R1='/scratch/data/bio/sra_reads/e_coli/SRR10561103_1.fastq.gz'
pe1_R2='/scratch/data/bio/sra_reads/e_coli/SRR10561103_2.fastq.gz'
# look for already indexed genome here /scratch/data/bio/genomes/indexes
ref_genome='/scratch/data/bio/genomes/indexes/e_coli/GCF_000005845.2_ASM584v2_genomic.fna'

######## PARAMETERS ########
readgroup='ecoli_sra'
sample='SRR10561103'
library='pe'
platform='ILLUMINA'         # ILLUMINA, CAPILLARY, LS454, SOLID, HELICOS, IONTORRENT, ONT, PACBIO
cpus=$SLURM_CPUS_PER_TASK

########## OUTPUTS #########
output_bam="${sample}_sorted.bam"

################################### COMMANDS ###################################
# NOTE index genome only if not using already indexed genome
if [ ! -f ${ref_genome}.bwt ]; then
	  bwa index $ref_genome
  fi

  bwa mem -M -t $cpus -R "@RG\tID:$readgroup\tLB:$library\tSM:$sample\tPL:$platform" $ref_genome $pe1_R1 $pe1_R2 -o $TMPDIR/${sample}_aln.sam

  samtools sort $TMPDIR/${sample}_aln.sam -o $output_bam -m 5G -@ $cpus -T $TMPDIR/tmp4sort

  ################################################################################
  <<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/wiki/index.php/HPRC:AckUs

    - BWA:
        Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754-1760.

    - SAMtools:
        Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project
        Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.
CITATION

