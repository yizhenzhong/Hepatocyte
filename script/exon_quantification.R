setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/data/expression/EQP_exon_count/")
library(DEXSeq)
library(edgeR)
library(GenABEL)
library(data.table)
library(preprocessCore)
#There is a bug in the qn function
qn <- function(M){#assume the row is gene and column is sample
  #modified from https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/rnaseqnorm.py
  Q = apply(M, 2, function(x) order(x))
  m = dim(M)[1]
  n = dim(M)[2]
  
  #compute quantile vector
  quantiles = rep(0, n)
  
  for(i in 1:n){
    quantiles = quantiles + M[Q[,i],i]
  }
  quantiles = quantiles/n
  
  
  
  for(i in 1:n){
    dupes = rep(0, m)
    for(j in 1:(m-1)){
      if(M[Q[j,i],i]==M[Q[j+1,i],i]){
        dupes[j+1] = dupes[j]+1
      }
    }
    
    #print(dupes)
    #replace column with quantile ranks
    M[Q[,i],i] = quantiles
    
    #print(M)
    
    j = m
    while(j >= 1){
      if(dupes[j] == 0){
        j = j-1
      }else{
        
        #print(dupes[j])
        idxs = Q[(j-dupes[j]):j,i]
        #print((j-dupes[j]:j)
        M[idxs,i] = mean(M[idxs,i])
        j = j - (dupes[j] + 1)
      }
    }
  } 
  return(M)
}


inverse_normal<-function(M){
  #transform column to a standard normal distribution
  M_norm = apply(M,2,rntransform)
  return(M_norm)
}

M = matrix(c(5,2,3,4,4,1,4,2,3,4,6,8),nrow=4,ncol=3)
normalize.quantiles(M,copy=F)
qn(M)

######################################################
#read the exon count file
#reorder the count files according to the sample order in the covariate file
countFiles = list.files("./", pattern=".out", full.names = F )
basename(countFiles)

temp = sapply(countFiles, function(x) strsplit(x, "-")[[1]][2])
countid = sapply(temp, function(x) substr(x, 1, nchar(x)-26))
countid[grep("MP-1Alig", countFiles)] = "0"  
covariate = read.csv("../../../data/factors_60_GSK.csv", stringsAsFactors = F)
covariate_id = sapply(covariate$ID, function(x) strsplit(x, "-")[[1]][2])
countFiles = countFiles[match(covariate_id, countid)]

###########################
cc.file = NA;
res = NULL
### Loop over files
for(t1 in countFiles ) { 
  
  
  start.time = proc.time()[3];
  tbl = read.table(t1, 
                   header = F, stringsAsFactors=FALSE, colClasses=cc.file);
  end.time = proc.time()[3];
  cat(t1, "loaded in", end.time - start.time, "sec.", nrow(tbl), "genes.", "\n");
  
  ### set colClasses for faster loading of other results
  if(any(is.na(cc.file))) {
    cc.file = sapply(tbl, class);
  }
  if(is.null(res)){res = tbl}else{res = cbind(res, tbl[,2])}

}  


############################
#read the gene length file
gene_pos = read.table("../gencode.v25lift37.exon.annotation.gtf",header = F,stringsAsFactors = F,sep = " ")
gene_pos_sub = gene_pos[match(res$V1,gene_pos$V4),]
all(res$V1 == gene_pos_sub$V4)

rownames(res) = res$V1
res = res[,-1]
colnames(res) = covariate$ID
edge_count <- DGEList(counts=res,genes = gene_pos_sub,samples = covariate)

save(edge_count,file="../../../results/edgeR_exon_count_baseline_all.rda")

##############################
##convert to RPKM, set the normalized lib size to be false because the normalization is then done by
# the quantile normalization
rpkm = rpkm(edge_count,gene.length = (edge_count$genes$V3 - edge_count$genes$V2), 
            normalized.lib.sizes = F)

#use the same filtering criteria as Yilin
index = as.vector(apply(rpkm,1,function(x) sum(x>=0.1) > 10))

rpkm_filter = rpkm[index, ]

######################################
#quantile normalization across samples
rpkm_filter_qn = normalize.quantiles(rpkm_filter)

#inverse quantile normalization for each gene
rpkm_filter_qn_int = t(inverse_normal(t(rpkm_filter_qn)))

#combine with the gene annotation file
colnames(rpkm_filter_qn_int) = colnames(edge_count$counts)

rpkm_filter_qn_int = cbind(edge_count$genes[index,], rpkm_filter_qn_int)
colnames(rpkm_filter_qn_int)[1:6] = c("chr","start","end","pid","gid","strand")

##################
rpkm_filter_qn_int = data.table(rpkm_filter_qn_int)
bed = rpkm_filter_qn_int[chr %in% paste0("chr",1:22)]
bed = bed[order(chr,start)] #order by chr and then by start

colnames(bed)[1] = "#chr"
write.table(bed, gzfile("../baseline_hg19_exon.bed.gz"), sep = "\t", row.names = F, quote = F)
