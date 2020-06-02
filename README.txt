My_UMI_tool

Tool for handling Unique Molecular Identifiers in NGS data sets to determine the absolute number of unique molecules by identifying duplicate reads.

Testing data available at: https://drive.google.com/drive/folders/1cqu8up0rY1XYA45hC0bq-P55ssEXXI3B?usp=sharing

To run My_UMI_tool from command line:

1. Install the Miniconda Python3 distribution
2. Install snakemake: conda install -c conda-forge -c bioconda snakemake
3. Create the environtment from the yaml file: conda env create -f env.yaml (you can find the yaml file in /umi_project/scripts/wraps/fastq2bam_RNA/clust_umi/env.yaml)
4. Activate the new environment: conda activate myumitool
5. Run command: Rscript clust_umi.R {input.fastq} {output.fastq} {output.grp} (direct from directory or write down full path of a file clust_umi.R)


To run the pipeline:

1. Install the Miniconda Python3 distribution
2. Install snakemake: conda install -c conda-forge -c bioconda snakemake
3. Download umi_project folder 
4. Download references:
  annot file	to	/umi_project/references/homsap/GRCh38-p10/annot 		
  SAindex file 	to 	/umi_project/references/homsap/GRCh38-p10/index/STAR/SAindex 	
  seq file 	to 	/umi_project/references/homsap/GRCh38-p10/seq 			

  annot file (unrar):			https://drive.google.com/drive/folders/1cHsKabMHTuJvxEPFXM9W4O_9qltEqhke
  SAindex file (unrar):			https://drive.google.com/drive/folders/1BAhhQdTLwJfklCvH5CwZul1moC5BvtC9
  seq file (unzip from fasta.gz to fa): 	https://drive.google.com/drive/folders/1fs2YJ0ecqDWvfNynRtq9kutKG2uWB4Y4

5. Run command: snakemake --use-conda --snakefile /umi_project/scripts/snakemake_test.smk -j (number of available cores) -r
