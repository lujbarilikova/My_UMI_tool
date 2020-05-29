
# REQUIRED LIBRARIES
library(data.table)
library(stringdist)

# DATA UPLOADING & PREPROCESSING 
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = T)
input_fastqfile <- args[1]      
output_fastqfile <- gsub(".gz", "", args[2])
sample_name <- gsub("[_/.].*","",basename(input_fastqfile)) 
output_grpfile <- args[3]
fastq <- readLines(input_fastqfile)
output_dir <- dirname(output_fastqfile) 

# preprocessing
head <- fastq[seq(from = 1, to = length(fastq), by = 4)]
head <- head[nchar(head) > 0]
head <- gsub("@", "", head)
seq <- fastq[seq(from = 2, to = length(fastq), by = 4)]
q <- fastq[seq(from = 4, to = length(fastq), by = 4)]
umi <- gsub(".*_([ACGT]+).{0,100}", "\\1", head ,perl=T)
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
lapply(list.files(paste0(output_dir,paste0("/",sample_name,"/umi_fa"))),function(i){
system(paste0("vsearch --threads 10 --cluster_fast ", paste0(output_dir,paste0("/",sample_name,"/umi_fa/")),i," --id 0.8 --iddef 1 --uc ", paste0(output_dir,paste0("/",sample_name,"/umi_clust/"),gsub(".fa", "",i),".txt"), " --consout ", paste0(output_dir,paste0("/",sample_name,"/consensus/"),gsub(".fa", "",i)), " --msaout ", paste0(output_dir,paste0("/",sample_name,"/alignment/"),gsub(".fa", "",i))))
})

# Determination of STARTING_GAP_COUNT
files_df <- rbindlist(lapply(list.files(paste0(output_dir,paste0("/",sample_name,"/umi_clust")), full.names=T), function(x) {
  tab <- fread(x)
  tab <- tab[V1 == "S" | V1 == "H"]  # all H (hit) and S (centroid) sequences
  seq_file <- sub(".txt","",sub("umi_clust","alignment",x)) # loading files from alignment directory
  seq <- readLines(seq_file)
  seq_tab <- data.table(header = gsub(">\\*?","",seq[grep("^>",seq)][grep(">consensus",seq[grep("^>",seq)],invert = T)]),starting_gap_count = nchar(gsub("(^-*).*","\\1",seq[grep("^>",seq) + 1]))[grep(">consensus",seq[grep("^>",seq)],invert = T)])
  tab <- merge(tab,seq_tab,by.x = "V9",by.y = "header") # merging to add starting_gap_count
  return(tab)
}))

files_df[,umi := gsub(".*_","",V9)]   # trimming UMI from V9 and making new row called umi
clust <- files_df[, list(clust = paste(umi,V2,starting_gap_count,sep = "_"),head = V9)]  # head contains umi_clusternumber_startinggapcount
clust[,head := as.character(head)]

merged <- merge(data,clust,by="head",all.x = T)
merged2 <- merged[,Nclust := .N,by = clust]   # couning number of sequences in clusts
merged2 <- merged2[!is.na(merged2$clust),]
hm <- round(mean(merged2$Nclust)*1.7) 
dm <- round(mean(merged2$Nclust)*0.3) 
big_clust <- merged2[merged2$Nclust>=hm,]     # counts of sequences above a particular empirically chosen threshold are selected
small_clust <- merged2[merged2$Nclust <= dm,] # counts of sequences under a particular empirically chosen threshold are selected

# UMI error correction
# all single-nucleotide different UMIs that could be observed for each unique combination of original UMI are generated
ACGT <- c("A","C","G","T")
tab <- CJ(0:3,0:3,0:3,0:3,0:3,0:3)
tab[,UMI := apply(tab,1,function(x) paste(ACGT[x + 1],collapse = ""))]
tab_melt <- melt(tab,id.vars = "UMI")
tab_melt <- tab_melt[,.(change = (value + 1:3) %% 4),by = c("UMI","variable")]
tab <- merge(tab,tab_melt,by = "UMI")

for(col_select in unique(tab$variable)){
  tab[variable == col_select,(col_select) := change]
}
tab[,UMI_near := apply(tab[,2:7,with = F],1,function(x) paste(ACGT[x + 1],collapse = ""))]
umi_near_tab <- tab[,.(UMI,UMI_near)]      #final tab of all single-nucleotide different UMIs 
setnames(umi_near_tab,c("umi","umi_near"))

# to finally cluster sequences with similar barcodes, pairwise string distances between read sequences from big_clust clusters and read sequence from the small_clust are computed
unique_big_clust <- unique(big_clust,by = "clust")  
clust_connect <- merge(umi_near_tab,unique_big_clust,by = "umi",allow.cartesian = T) 
clust_connect <- clust_connect[,.(clust,umi_near,Nclust)]
res <- merge(clust_connect,small_clust,by.x = "umi_near",by.y = "umi",allow.cartesian = T) 
res2 <- res[,.SD[stringdist(unique_big_clust[clust == .BY[[1]]]$seq,seq) < 2],by = clust.x]   # applying stringdist with the number of dissimilarities (2 by default)

# postprocessing
setnames(res2, "clust.x", "clust")
res2 <- res2[,c("clust","head")]
new_clust <- merge(merged2,res2,by="head",all.x = T)
new_clust[!is.na(new_clust$clust.y),]$clust.x = new_clust[!is.na(new_clust$clust.y),]$clust.y
setnames(new_clust, "clust.x", "clust")

# results exported to a TSV file
my_umi_grouped <- new_clust[,c("head","clust")]
my_umi_grouped <- data.frame(new_clust)
write.table(my_umi_grouped,file = output_grpfile, row.names = FALSE)

#generating final deduplicated FASTQ file
final_merged <- unique(new_clust,by = "clust")
cat(paste0("@",final_merged$head,"\n",final_merged$seq,"\n+\n",final_merged$q,"\n"),file = output_fastqfile,sep = "")
system(paste0("gzip ",output_fastqfile))

