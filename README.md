# Prerequisite 
- Sequence: '{viral/host}.fa' files in 'target_seqs' directory
- ML: 'combined_classifier.pth', 'combined_regressor.pth', and 'standard_scaler.joblib' in 'src' directory
- Parameter file: 'params.txt' in the current directory
- Scripts: 'generate_primers.py', 'prepare_input.py', 'evaluate_primers.py', and 'build_final_output.py' in 'scripts' directory
> Software: snakemake, bowtie2, sam2pairwise, ViennaRNA

# Run 
- snakemake --resources gpu=1 --cores all
