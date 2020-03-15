# My_UMI_tool
Tool for handling Unique Molecular Identifiers in NGS data sets to determine the absolute number of unique molecules by identifying duplicate reads.

Testing data available at: https://drive.google.com/drive/folders/1cqu8up0rY1XYA45hC0bq-P55ssEXXI3B?usp=sharing

To run the pipeline you need to:
1. Download snakemake: https://snakemake.readthedocs.io/en/stable/
2. Use command: snakemake --use-conda  --snakefile .../scripts/snakemake_test.smk -r --nolock
