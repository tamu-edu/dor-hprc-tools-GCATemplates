#!/bin/bash
#SBATCH --job-name=alphafold        # job name
#SBATCH --time=2-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=16          # CPUs (threads) per command
#SBATCH --mem=96G                   # total memory per node
#SBATCH --gres=gpu:a30:1            # request one GPU
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file

<<README
    - AlphaFold manual: https://github.com/deepmind/alphafold

README

######### SYNOPSIS #########
# this script will run the alphafold singularity container
# currently AlphaFold supports running on only one GPU

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
protein_fasta='/scratch/data/bio/alphafold/example_data/IL2Y.fasta'

######## PARAMETERS ########
DOWNLOAD_DIR='/scratch/data/bio/alphafold/2.3.0'  # 3.4TB data already downloaded here
max_template_date='2023-1-1'
model_preset='monomer'              # monomer, monomer_casp14, monomer_ptm, multimer
db_preset='full_dbs'                # full_dbs, reduced_dbs

########## OUTPUTS #########
output_dir='out_alphafold_2.3.1'

################################### COMMANDS ###################################

jobstats &

singularity exec --nv /sw/hprc/sw/bio/containers/alphafold/alphafold_2.3.1.sif python /app/alphafold/run_alphafold.py  \
  --data_dir=$DOWNLOAD_DIR  --use_gpu_relax \
  --uniref90_database_path=$DOWNLOAD_DIR/uniref90/uniref90.fasta  \
  --mgnify_database_path=$DOWNLOAD_DIR/mgnify/mgy_clusters_2022_05.fa  \
  --bfd_database_path=$DOWNLOAD_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  \
  --uniref30_database_path=$DOWNLOAD_DIR/uniclust30/UniRef30_2022_02 \
  --pdb70_database_path=$DOWNLOAD_DIR/pdb70/pdb70  \
  --template_mmcif_dir=$DOWNLOAD_DIR/pdb_mmcif/mmcif_files  \
  --obsolete_pdbs_path=$DOWNLOAD_DIR/pdb_mmcif/obsolete.dat \
  --model_preset=$model_preset \
  --max_template_date=$max_template_date \
  --db_preset=$db_preset \
  --output_dir=$output_dir \
  --fasta_paths=$protein_fasta

jobstats

################################################################################
<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - AlphaFold:
        Jumper, John et al. "Highly accurate protein structure prediction with AlphaFold". Nature 596. 7873(2021): 583–589.
CITATIONS

<<USAGE
Full AlphaFold protein structure prediction script.
flags:

/app/alphafold/run_alphafold.py:
  --[no]benchmark: Run multiple JAX model evaluations to obtain a timing that excludes the compilation time, which should be more indicative of the time required for inferencing many proteins.
    (default: 'false')
  --bfd_database_path: Path to the BFD database for use by HHblits.
  --data_dir: Path to directory of supporting data.
  --db_preset: <full_dbs|reduced_dbs>: Choose preset MSA database configuration - smaller genetic database config (reduced_dbs) or full genetic database config  (full_dbs)
    (default: 'full_dbs')
  --fasta_paths: Paths to FASTA files, each containing a prediction target that will be folded one after another. If a FASTA file contains multiple sequences, then it will be folded as a multimer. Paths should be separated by
    commas. All FASTA paths must have a unique basename as the basename is used to name the output directories for each prediction.
    (a comma separated list)
  --hhblits_binary_path: Path to the HHblits executable.
    (default: '/usr/bin/hhblits')
  --hhsearch_binary_path: Path to the HHsearch executable.
    (default: '/usr/bin/hhsearch')
  --hmmbuild_binary_path: Path to the hmmbuild executable.
    (default: '/usr/bin/hmmbuild')
  --hmmsearch_binary_path: Path to the hmmsearch executable.
    (default: '/usr/bin/hmmsearch')
  --jackhmmer_binary_path: Path to the JackHMMER executable.
    (default: '/usr/bin/jackhmmer')
  --kalign_binary_path: Path to the Kalign executable.
    (default: '/usr/bin/kalign')
  --max_template_date: Maximum template release date to consider. Important if folding historical test sets.
  --mgnify_database_path: Path to the MGnify database for use by JackHMMER.
  --model_preset: <monomer|monomer_casp14|monomer_ptm|multimer>: Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model
    (default: 'monomer')
  --num_multimer_predictions_per_model: How many predictions (each with a different random seed) will be generated per model. E.g. if this is 2 and there are 5 models then there will be 10 predictions per input. Note: this
    FLAG only applies if model_preset=multimer
    (default: '5')
    (an integer)
  --obsolete_pdbs_path: Path to file containing a mapping from obsolete PDB IDs to the PDB IDs of their replacements.
  --output_dir: Path to a directory that will store the results.
  --pdb70_database_path: Path to the PDB70 database for use by HHsearch.
  --pdb_seqres_database_path: Path to the PDB seqres database for use by hmmsearch.
  --random_seed: The random seed for the data pipeline. By default, this is randomly generated. Note that even if this is set, Alphafold may still not be deterministic, because processes like GPU inference are
    nondeterministic.
    (an integer)
  --[no]run_relax: Whether to run the final relaxation step on the predicted models. Turning relax off might result in predictions with distracting stereochemical violations but might help in case you are having issues with
    the relaxation stage.
    (default: 'true')
  --small_bfd_database_path: Path to the small version of BFD used with the "reduced_dbs" preset.
  --template_mmcif_dir: Path to a directory with template mmCIF structures, each named <pdb_id>.cif
  --uniprot_database_path: Path to the Uniprot database for use by JackHMMer.
  --uniref30_database_path: Path to the UniRef30 database for use by HHblits.
  --uniref90_database_path: Path to the Uniref90 database for use by JackHMMER.
  --[no]use_gpu_relax: Whether to relax on GPU. Relax on GPU can be much faster than CPU, so it is recommended to enable if possible. GPUs must be available if this setting is enabled.
  --[no]use_precomputed_msas: Whether to read MSAs that have been written to disk instead of running the MSA tools. The MSA files are looked up in the output directory, so it must stay the same between multiple runs that are
    to reuse the MSAs. WARNING: This will not check if the sequence, database or configuration have changed.
    (default: 'false')

absl.app:
  -?,--[no]help: show this help
    (default: 'false')
  --[no]helpfull: show full help
    (default: 'false')
  --[no]helpshort: show this help
    (default: 'false')
  --[no]helpxml: like --helpfull, but generates XML output
    (default: 'false')
  --[no]only_check_args: Set to true to validate args and exit.
    (default: 'false')
  --[no]pdb: Alias for --pdb_post_mortem.
    (default: 'false')
  --[no]pdb_post_mortem: Set to true to handle uncaught exceptions with PDB post mortem.
    (default: 'false')
  --profile_file: Dump profile information to a file (for python -m pstats). Implies --run_with_profiling.
  --[no]run_with_pdb: Set to true for PDB debug mode
    (default: 'false')
  --[no]run_with_profiling: Set to true for profiling the script. Execution will be slower, and the output format might change over time.
    (default: 'false')
  --[no]use_cprofile_for_profiling: Use cProfile instead of the profile module for profiling. This has no effect unless --run_with_profiling is set.
    (default: 'true')

absl.logging:
  --[no]alsologtostderr: also log to stderr?
    (default: 'false')
  --log_dir: directory to write logfiles into
    (default: '')
  --logger_levels: Specify log level of loggers. The format is a CSV list of 'name:level'. Where 'name' is the logger name used with 'logging.getLogger()', and 'level' is a level name  (INFO, DEBUG, etc). e.g.
    'myapp.foo:INFO,other.logger:DEBUG'
    (default: '')
  --[no]logtostderr: Should only log to stderr?
    (default: 'false')
  --[no]showprefixforinfo: If False, do not prepend prefix to info messages when it's logged to stderr, --verbosity is set to INFO level, and python logging is used.
    (default: 'true')
  --stderrthreshold: log messages at this level, or more severe, to stderr in addition to the logfile.  Possible values are 'debug', 'info', 'warning', 'error', and 'fatal'.  Obsoletes --alsologtostderr. Using
    --alsologtostderr cancels the effect of this flag. Please also note that this flag is subject to --verbosity and requires logfile not be stderr.
    (default: 'fatal')
  -v,--verbosity: Logging verbosity level. Messages logged at this level or lower will be included. Set to 1 for debug logging. If the flag was not set or supplied, the value will be changed from the default of -1 (warning) to
    0 (info) after flags are parsed.
    (default: '-1')
    (an integer)

absl.testing.absltest:
  --test_random_seed: Random seed for testing. Some test frameworks may change the default value of this flag between runs, so it is not appropriate for seeding probabilistic tests.
    (default: '301')
    (an integer)
  --test_randomize_ordering_seed: If positive, use this as a seed to randomize the execution order for test cases. If "random", pick a random seed to use. If 0 or not set, do not randomize test case execution order. This flag
    also overrides the TEST_RANDOMIZE_ORDERING_SEED environment variable.
    (default: '')
  --test_srcdir: Root of directory tree where source files live
    (default: '')
  --test_tmpdir: Directory for temporary testing files
    (default: '/tmp/absl_testing')
  --xml_output_file: File to store XML test results
    (default: '')

tensorflow.python.ops.parallel_for.pfor:
  --[no]op_conversion_fallback_to_while_loop: DEPRECATED: Flag is ignored.
    (default: 'true')

tensorflow.python.tpu.client.client:
  --[no]hbm_oom_exit: Exit the script when the TPU HBM is OOM.
    (default: 'true')
  --[no]runtime_oom_exit: Exit the script when the TPU runtime is OOM.
    (default: 'true')

absl.flags:
  --flagfile: Insert flag definitions from the given file into the command line.
    (default: '')
  --undefok: comma-separated list of flag names that it is okay to specify on the command line even if the program does not define a flag with that name.  IMPORTANT: flags in this list that have arguments MUST use the
    --flag=value format.
    (default: '')
USAGE
