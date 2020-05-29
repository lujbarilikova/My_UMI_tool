# My_UMI_tool
Tool for handling Unique Molecular Identifiers in NGS data sets to determine the absolute number of unique molecules by identifying duplicate reads.

Testing data available at: https://drive.google.com/drive/folders/1cqu8up0rY1XYA45hC0bq-P55ssEXXI3B?usp=sharing


To run My_UMI_tool from command line:
  1. Install the Miniconda Python3 distribution
  2. Install snakemake:
	    conda install -c conda-forge -c bioconda snakemake
  3. Create the environtment from the yaml file:		
	    conda env create -f env.yaml  (you can find the yaml file in My_UMI_tool/scripts/wraps/fastq2bam_RNA/clust_umi/env.yaml)
  4. Activate the new environment:
	    conda activate myumitool
  5. Run command:
	    Rscript /mnt/ssd/ssd_3/temp/lujza/umi_project/scripts/wraps/fastq2bam_RNA/clust_umi/clust_umi.R {input.fastq} {output.fastq} {output.grp}
    

To run the pipeline:
  1. Install the Miniconda Python3 distribution
  2. Install snakemake:
	    conda install -c conda-forge -c bioconda snakemake
  3. Run command:

