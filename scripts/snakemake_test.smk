REF_DIR = "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10"
workdir: "/mnt/ssd/ssd_3/temp/lujza/umi_project"
WDIR = "/mnt/ssd/ssd_3/temp/lujza/umi_project"
sample_list = ["L102-Q_SE","L102FFPE-Q_SE"]

rule all_res:
    input: #myumi = expand(WDIR + "/myumi_dedup_fastq/{sample}_grouped.tsv",sample = sample_list)
           bam_umitools = expand(WDIR + "/test/mapped/{sample}_umitools_dedup.bam",sample = sample_list),
           bam_myumi = expand(WDIR + "/test/mapped/{sample}.myumi_dedup.bam",sample = sample_list)

rule mark_duplicates_umi_tool:
    input:  bam = WDIR + "/test/mapped/{sample}.normal.bam"
    output: bam = WDIR + "/test/mapped/{sample}_umitools_dedup.bam",
            bai = WDIR + "/test/mapped/{sample}_umitools_dedup.bam.bai",
            mtx = WDIR + "/test/postQC/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
            grp = WDIR + "/test/mapped/{sample}_umitools_dedup_grouped.bam",
            grptsv = WDIR + "/test/mapped/{sample}_umitools_dedup_grouped.tsv"
    log:    run = WDIR + "/sample_logs/{sample}/mark_duplicates.log"
    threads:  8
    resources:  mem = 15
    params: rmDup = "false"
    conda:  "wraps/fastq2bam_RNA/mark_duplicates/env.yaml"
    shell: """
          umi_tools dedup -I {input.bam} -S {output.bam} --log {output.mtx} \
          --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0 \
          --spliced-is-unique --multimapping-detection-method=NH
          samtools index -@ {threads} {output.bam}  
          umi_tools group -I {input.bam} --group-out={output.grptsv} --output-bam -S {output.grp} 
    """

        # umi_tools dedup -I {input.bam} -S {output.bam} --log {output.mtx} \
        # --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0 \
        # --spliced-is-unique --multimapping-detection-method=NH
        # samtools index -@ {threads} {output.bam}   
   

 
rule alignment_SE:
    input:  #fastq = lambda wildcards: expand("{type}_fastq/{sample}.fastq.gz", type = ["normal", "myumi_dedup"], sample = wildcards.sample),
            fastq = WDIR + "/{type}_fastq/{sample}.fastq.gz",
            genome = REF_DIR + "/seq/GRCh38-p10.fa",
            gtf = REF_DIR + "/annot/GRCh38-p10.gtf",
            index = REF_DIR + "/index/STAR/SAindex",
            gendir = REF_DIR + "/index/STAR",
    output: bam = WDIR + "/test/mapped/{sample}.{type}.bam",
            bai = WDIR + "/test/mapped/{sample}.{type}.bam.bai"
    log:    run = WDIR + "/sample_logs/{sample}/alignment_{type}_SE.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "test/mapped/{sample}/{sample}",
            strandness = "fwd", #"true" or "false" - STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM
            perc_mismatch = 0.1, # percentage of mismatches between reference and MAPPED part of read – 0.05, maybe much strict, should be set to 0.1, default is 0.3 which is too loose.
            num_mismatch = 999, # Maximum number of mismatches; set this to high number (999) to disable and to use only perc_mismatch
            max_intron = 1000000, # Default used by ENCODE is 1000000; to turn this off set 1
            max_mate_dist = 1000000, # Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this
            read_len = 100, # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
    conda:  "wraps/fastq2bam_RNA/alignment_SE_RNA/env.yaml"
    shell: """
        mkdir -p $(dirname {params.prefix}) || echo $? > /dev/null;
        STAR --runMode alignReads --runThreadN {threads} \
               --genomeDir "$(dirname {input.index})" \
               --readFilesIn {input.fastq} \
               --readFilesCommand zcat \
               --sjdbOverhang {params.read_len} \
               --sjdbGTFfile {input.gtf} \
               --outFileNamePrefix  {params.prefix} \
               --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
               --outFilterMismatchNmax {params.num_mismatch} \
               --outFilterMismatchNoverReadLmax 1.0 \
               --outFilterMismatchNoverLmax {params.perc_mismatch} \
               --alignIntronMin 20 --alignIntronMax {params.max_intron} \
               --alignMatesGapMax {params.max_mate_dist} \
               --outFilterMatchNmin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 \
               --outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
               --outSAMattrRGline ID: {wildcards.sample} PL:Illumina PU: {wildcards.sample} SM: {wildcards.sample}\
               --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All \
               extra_flags_star_motif --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic \
               --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate >> {log.run} 2>&1 &&
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} &&
        samtools index -@ {threads} {output.bam} {output.bai}
  """

rule myumi:
    input:  fastq = WDIR + "/normal_fastq/{sample}.fastq.gz",
    output: fastq = WDIR + "/myumi_dedup_fastq/{sample}.fastq.gz",
            grp = WDIR + "/myumi_dedup_fastq/{sample}_grouped.tsv",
    conda:  WDIR + "/scripts/wraps/fastq2bam_RNA/clust_umi/env.yaml"
    shell: """
        Rscript /mnt/ssd/ssd_3/temp/lujza/umi_project/scripts/wraps/fastq2bam_RNA/clust_umi/clust_umi.R {input.fastq} {output.fastq} {output.grp}
  """