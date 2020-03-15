
######################################
# wrapper for rule: alignment_SE_RNA
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

STAR = "STAR"
SAMTOOLS = "samtools"
BEDGRAPH2BIGWIG = "bedGraphToBigWig"

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: alignment_SE_RNA \n##\n")
f.close()


version = str(subprocess.Popen(STAR+" --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

version = str(subprocess.Popen(SAMTOOLS+" --version 2>&1 | grep samtools", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.params.prefix)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# if os.path.isdir(snakemake.params.dir+"_STARtmp") :
#   command = "rm -rf "+snakemake.params.dir+"_STARtmp >> "+snakemake.log.run+" 2>&1"
#   f = open(snakemake.log.run, 'at')
#   f.write("## COMMAND: "+command+"\n")
#   f.close()
#   shell(command)
star_index_dir = snakemake.input.index.replace("/SAindex","")

f = open(snakemake.log.run, 'at')

strandness = False if snakemake.params.strandness == "unstr" else True

if strandness == True:
    extra_flags_star_motif = "" # For STAR outSAMstrandField
    extra_flags_star_wig = " --outWigStrand Stranded" # For START bedGraph
    f.write("Running as Stranded experiment \n")
else:
    extra_flags_star_motif = " --outSAMstrandField intronMotif" # For STAR outSAMstrandField
    extra_flags_star_wig = " --outWigStrand Unstranded" # For START bedGraph
    f.write("Running as Unstranded experiment \n")
f.close()

command = STAR+" --runMode alignReads --runThreadN " + str(snakemake.threads) + \
               " --genomeDir " + star_index_dir + \
               " --readFilesIn " + str(snakemake.input.fastq)  + \
               " --readFilesCommand zcat" + \
               " --sjdbOverhang " + str(snakemake.params.read_len) + \
               " --sjdbGTFfile " + str(snakemake.input.gtf) + \
               " --outFileNamePrefix " + snakemake.params.prefix + \
               " --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1" + \
               " --outFilterMismatchNmax " + str(snakemake.params.num_mismatch) + \
               " --outFilterMismatchNoverReadLmax 1.0" + \
               " --outFilterMismatchNoverLmax "+ str(snakemake.params.perc_mismatch) + \
               " --alignIntronMin 20 --alignIntronMax " + str(snakemake.params.max_intron) + \
               " --alignMatesGapMax " + str(snakemake.params.max_mate_dist) +\
               " --outFilterMatchNmin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 " + \
               " --outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold" + \
               " --outSAMattrRGline ID:"+str(snakemake.wildcards.sample)+" PL:Illumina PU:"+str(snakemake.wildcards.sample)+" SM:"+str(snakemake.wildcards.sample) + \
               " --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All" + \
               extra_flags_star_motif +" --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic " + \
               " --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

shell("mv " + snakemake.params.prefix + "Aligned.sortedByCoord.out.bam " + snakemake.output.bam)
# shell("mv " + snakemake.params.prefix + "Aligned.toTranscriptome.out.bam " + snakemake.output.transcriptom_bam)

# shell("rm -r " + snakemake.params.prefix + "*pass1")

command = SAMTOOLS +" index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam + " " + snakemake.output.bai
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

