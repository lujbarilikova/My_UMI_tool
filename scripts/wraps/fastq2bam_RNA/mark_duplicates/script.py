######################################
# wrapper for rule: mark_duplicates
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: mark_duplicates \n##\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.output.mtx)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


# if snakemake.params.umi[0] == "no":
# TOOL = "picard MarkDuplicates"
# SAMTOOLS = "samtools"

# shell.executable("/bin/bash")

# version = str(subprocess.Popen(TOOL+" --version 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(snakemake.log.run, 'at')
# f.write("## VERSION: Picard "+version+"\n")
# f.close()

# command = "export LD_BIND_NOW=1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# command = TOOL+" INPUT="+snakemake.input.bam+" OUTPUT="+snakemake.output.bam+" METRICS_FILE="+snakemake.output.mtx+" REMOVE_DUPLICATES="+snakemake.params.rmDup+" \
#     ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT -Xmx"+str(snakemake.resources.mem)+"g 2>> "+snakemake.log.run+" "
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# command = SAMTOOLS+" index "+snakemake.output.bam
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# command = "rm "+snakemake.input.bam
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

# command = "rm "+snakemake.input.bam + ".bai"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)



# else:
TOOL = "umi_tools"
SAMTOOLS = "samtools"

shell.executable("/bin/bash")

version = str(subprocess.Popen(SAMTOOLS+" --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = TOOL + " dedup -I " + snakemake.input.bam + " -S " + snakemake.output.bam + " --log " + snakemake.output.mtx \
+ " --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0 \
--spliced-is-unique --multimapping-detection-method=NH"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS + " index -@" + str(snakemake.threads) + " " + snakemake.output.bam
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
