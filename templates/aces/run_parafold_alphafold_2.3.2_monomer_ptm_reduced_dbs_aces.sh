#!/bin/bash
#SBATCH --job-name=pf-cpu           # job name
#SBATCH --time=1-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=32G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

<<README
    - AlphaFold manual: https://github.com/deepmind/alphafold
    - ParaFold: https://github.com/Zuricho/ParallelFold
    - AlphaPickle: https://github.com/mattarnoldbio/alphapickle
README

######### SYNOPSIS #########
# this script will run AlphaFold in two jobs one for the CPU step and a second for the GPU steps and graph .pkl files.
# currently AlphaFold supports running on only one GPU

module purge
module load GCC/11.3.0  OpenMPI/4.1.4 AlphaFold/2.3.2-CUDA-11.8.0
module load ParaFold/2.0-CUDA-11.8.0

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
protein_fasta='/scratch/data/bio/alphafold/example_data/1L2Y.fasta'

######## PARAMETERS ########
ALPHAFOLD_DATA_DIR=/scratch/data/bio/alphafold/2.3.2  # 3TB data already downloaded here
max_template_date=2025-1-1
model_preset=monomer_ptm          # monomer, monomer_casp14, monomer_ptm, multimer
db_preset=reduced_dbs             # full_dbs, reduced_dbs

########## OUTPUTS #########
protein_basename=$(basename ${protein_fasta%.*})
pf_output_dir=out_parafold_${protein_basename}_${model_preset}_${db_preset}
pickle_out_dir=${pf_output_dir}/${protein_basename}

################################### COMMANDS ###################################
# First, run CPU-only steps to get multiple sequence alignments:
# Second, run GPU steps as a separate job after the first part completes successfully:

run_alphafold.sh -d $ALPHAFOLD_DATA_DIR -o $pf_output_dir -p $model_preset -i $protein_fasta -t $max_template_date -c $db_preset -f && \
sbatch --job-name=parafold-gpu --time=01:00:00 --ntasks-per-node=1 --cpus-per-task=24 --mem=122G --gres=gpu:h100:1 --partition=gpu_debug --output=stdout.%x.%j --error=stderr.%x.%j --dependency=afterok:$SLURM_JOBID<<EOF
#!/bin/bash
module purge
module load GCC/11.3.0  OpenMPI/4.1.4 AlphaFold/2.3.2-CUDA-11.8.0
module load ParaFold/2.0-CUDA-11.8.0 AlphaPickle/1.4.1
run_alphafold.sh -g -u 0 -d $ALPHAFOLD_DATA_DIR -o $pf_output_dir -p $model_preset -i $protein_fasta -t $max_template_date -c $db_preset
# graph pLDDT and PAE .pkl files
run_AlphaPickle.py -od $pickle_out_dir
EOF

################################################################################
<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - AlphaFold:
        Jumper, John et al. "Highly accurate protein structure prediction with AlphaFold". Nature 596. 7873(2021): 583–589.

        Tunyasuvunakool, Kathryn et al. "Highly accurate protein structure prediction for the human proteome". Nature 596. 7873(2021): 590–596.

    - AlphaFold-multimer:
        Evans, R et al. Protein complex prediction with AlphaFold-Multimer, doi.org/10.1101/2021.10.04.463034

    - ParaFold:
        Bozitao Zhong, Xiaoming Su, Minhua Wen, Sichen Zuo, Liang Hong, James Lin. ParaFold: Paralleling AlphaFold for Large-Scale Predictions. 
        2021. arXiv:2111.06340. doi.org/10.48550/arXiv.2111.06340

    - AlphaPickle
        Arnold, M. J. (2021) AlphaPickle, doi.org/10.5281/zenodo.5708709
CITATIONS
