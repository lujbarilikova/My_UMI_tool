# a <- Sys.time()
# fastq <- readLines('/mnt/ssd/ssd_3/temp/lujza/umi_project/normal_fastq/L102-Q_SE.fastq.gz')

# DATA UPLOADING & PREPROCESSING 
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = T)
input_fastqfile <- args[1]      
output_fastqfile <- gsub(".gz", "", args[2])
sample_name <- gsub("_.*","",basename(input_fastqfile))
output_grpfile <- args[3]
fastq <- readLines(input_fastqfile)
output_dir <- dirname(output_fastqfile)  
# preprocessing
head <- fastq[seq(from = 1, to = length(fastq), by = 4)]
head <- head[nchar(head) > 0]
head <- gsub("@", "", head)
seq <- fastq[seq(from = 2, to = length(fastq), by = 4)]
q <- fastq[seq(from = 4, to = length(fastq), by = 4)]
umi <- gsub(".*_([ACGT]+) .*", "\\1", head ,perl=T)
head <- gsub(" .*", "", head ,perl=T)
data <- data.table(head,seq,q,umi)
head <- regmatches(head[seq(from = 1, to = length(head))],regexpr(".*_([ACGT])*\\S",head[seq(from = 1, to = length(head))]))

# CLUSTERING by UMIs
dir.create(paste0(output_dir,paste0("/",sample_name,"/umi_fa")),recursive = TRUE)
data[,N := .N,by = umi] 
data[N > 1,cat(paste0(">",head,"\n",seq,"\n"),file = paste0(output_dir,paste0("/",sample_name,"/umi_fa/"),umi[1],".fa"),sep = ""),by = umi]

# VSEARCH clustering + alignment cutoff
dir.create(paste0(output_dir,paste0("/",sample_name,"/umi_clust")))
dir.create(paste0(output_dir,paste0("/",sample_name,"/consensus")))
dir.create(paste0(output_dir,paste0("/",sample_name,"/alignment")))
dir.create(paste0(output_dir,paste0("/",sample_name,"/alignment/score")))
lapply(list.files(paste0(output_dir,paste0("/",sample_name,"/umi_fa"))),function(i){
system(paste0("vsearch --userout ", paste0(output_dir,paste0("/",sample_name,"/alignment/score/"),gsub(".fa", "",i))  ," --threads 10 --cluster_fast ", paste0(output_dir,paste0("/",sample_name,"/umi_fa/")),i," --id 0.8 --uc ", paste0(output_dir,paste0("/",sample_name,"/umi_clust/"),gsub(".fa", "",i),".txt"), " --consout ", paste0(output_dir,paste0("/",sample_name,"/consensus/"),gsub(".fa", "",i)), " --msaout ", paste0(output_dir,paste0("/",sample_name,"/alignment/"),gsub(".fa", "",i))))
})


# REMOVING seqs according to STARTING_GAP_COUNT
length(grep("^>",seq))
files_df <- rbindlist(lapply(list.files(paste0(output_dir,paste0("/",sample_name,"/umi_clust")), full.names=T), function(x) {
  tab <- fread(x)
  tab <- tab[V1 == "S" | V1 == "H"]
  seq_file <- sub(".txt","",sub("umi_clust","alignment",x))
  seq <- readLines(seq_file)
  seq_tab <- data.table(header = gsub(">\\*?","",seq[grep("^>",seq)][grep(">consensus",seq[grep("^>",seq)],invert = T)]),starting_gap_count = nchar(gsub("(^-*).*","\\1",seq[grep("^>",seq) + 1]))[grep(">consensus",seq[grep("^>",seq)],invert = T)])
  tab <- merge(tab,seq_tab,by.x = "V9",by.y = "header")
  return(tab)
  }))

# files_df <- files_df[V1 == "S" | V1 == "H"]
files_df[,umi := gsub(".*_","",V9)]
clust <- files_df[, list(clust = paste(umi,V2,starting_gap_count,sep = "_"),head = V9)]
clust[,head := as.character(head)]
clust <- unique(clust)
my_umi_grouped <- data.frame(clust)
write.table(my_umi_grouped,file = output_grpfile, row.names = FALSE)


merged <- merge(data,clust,by="head",all.x = T)
save(merged,file = "merged.Rdata")
merged <- unique(merged,by = c("umi","clust"))
# merged[merged$clust == 0 | is.na(merged$clust), cat(paste0(">",head,"\n",seq,"\n"),file = paste0(output_dir,"/final/","final.fa"),sep = "")]
merged[merged$clust == 0 | is.na(merged$clust), cat(paste0("@",head,"\n",seq,"\n+\n",q,"\n"),file = output_fastqfile,sep = "")]
system(paste0("gzip ",output_fastqfile))

