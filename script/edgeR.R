setwd("/Volumes/fsmresfiles//Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(vsn)
library(edgeR)
library(pheatmap)
library(limma)
library(biomaRt)
library("BiocParallel")
drug = c('Baseline', 'Omeprazole','Phenob`arbital','Dexamethasone','Carbamazepin','Phenytoin','Rifampicin')
condit = paste0("condition",1:7)
#register(SnowParam(4))
#count = read.table("../data/liver_RNAseq_count1.txt",sep = "\t",stringsAsFactors = F,header = T)
###
#mart <- useEnsembl("ensembl")
#ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#gene_pos = getBM(attributes=c('ensembl_gene_id',"chromosome_name",
#                              "start_position","end_position","gene_biotype","hgnc_symbol"), filters =
#                   "ensembl_gene_id", values=gene_pos$Gene.stable.ID, mart=ensembl,
#                 bmHeader=TRUE)
#gene_pos$length = gene_pos$`Gene end (bp)`- gene_pos$`Gene start (bp)`
#gene_pos = gene_pos[!(duplicated(gene_pos$`Gene stable ID`)),]

#write.table(gene_pos,"../data/gene_pos_matrix.txt",col.names = T,row.names = F,quote = F,sep = "\t")
##############
gene_pos = read.table("../data/gene_pos_matrix.txt",header = T,stringsAsFactors = F,sep = "\t")
#################Merge counts
sample = NULL
exp = NULL
for(i in c(1:7)){
  fn = paste0("../data/liver_RNAseq_count",i,".txt")
  f =  read.table(fn,stringsAsFactors = F,sep = "\t",header = T)
  print(paste(i,"condition, expression dimension:", dim(f)[1],  dim(f)[2],sep=" "))
  f$gene = as.vector(sapply(f$gene,function(x) strsplit(x,".",fixed = T)[[1]][1]))
  f = f[match(gene_pos$Gene.stable.ID,f$gene),]
  temp = as.data.frame(sapply(colnames(f)[-1],function(x) substr(x, 1, nchar(x)-6)))
  temp$condition = i
  colnames(temp)[1] = "Hepatocyte"
  id = sapply(colnames(f)[-1],function(x) substr(x, 1, nchar(x)-6))
  sample = rbind(sample,temp)
  if(is.null(exp)){exp = f}else{exp = merge(exp,f,by="gene")}
}


count = data.frame(exp[,-1], row.names=exp[,1])
count = count[-c(1:5),]

#################get gene position info
gene_pos_sub = gene_pos[match(rownames(count),gene_pos$Gene.stable.ID),]
all(rownames(count) == gene_pos_sub$Gene.stable.ID)
edge_count <- DGEList(counts=count,genes = gene_pos_sub,samples = sample)

##############################filtering samples

#fam = read.table("../data/genotype/baseline.genotype.eqtl.excludeoutlier.header.fam",stringsAsFactors = F,header = F)
count_id = sapply(as.character(sample$Hepatocyte),function(x) strsplit(x, ".",fixed = T)[[1]][2])
count_id[sample$Hepatocyte=="Sample_MP."] = "0"
clinic = read.csv("../data/Hepatocytes_clinic_info_all.csv",stringsAsFactors = F)
clinic_id = sapply(as.character(clinic$ID),function(x) strsplit(x, "-",fixed = T)[[1]][2])
edge_count$samples =  cbind(edge_count$samples, clinic[match(count_id, clinic_id),2:4])
edge_count = edge_count[,-which(count_id %in% c("38","56","46","73","48","5")),keep.lib.sizes=F] #filter samples
save(edge_count,file="../results/edgeR_count_filter_429samples.rda")
#edge_count = get(load(file="../results/edgeR_count_filter_429samples.rda"))

###################################
#filter genes

tpm <- do.call(cbind, lapply(1:ncol(edge_count$counts), function(i) {
  rate = log(edge_count$counts[,i]) - log(gene_pos_sub$length)
  denom = log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}))


index1 = apply(tpm,1,function(x) sum(x>0.1) > ncol(tpm)*0.1)
index2 = as.vector(apply(edge_count$counts,1,function(x) sum(x>=6) > ncol(edge_count$counts)*0.1))
index = sapply(1:length(index1), function(x) all(index1[x],index2[x]))

edge_count = edge_count[index,,keep.lib.sizes=F]
edge_count = edge_count[which(edge_count$genes$Gene.type %in% c("protein_coding", "lincRNA")),,keep.lib.sizes=F]
save(edge_count,file="../results/edgeR_count_15575_429.rda")
########################################################normalize library size
edge_count = calcNormFactors(edge_count, method = "TMM")
library_size =  edge_count$samples$lib.size * edge_count$samples$norm.factors
norm_counts= apply(edge_count$counts,1,function(x) x/library_size*1e6)
tmm = list()
tmm[["norm_counts"]] = norm_counts
tmm[["sample"]] = edge_count$samples
tmm[["gene"]] = edge_count$genes
save(tmm,file="../results/edgeR_TMM_15575_429.rda")

marker1 = list(color =brewer.pal(7,"Dark2"))$color
colors = marker1[edge_count$samples$condition]
lcpm = cpm(edge_count,log=T)

pdf("../figure/limma_mds.pdf")
plotMDS(edge_count, labels =as.character(edge_count$samples$condition),
        gene.selection = "pairwise", col = colors, top = 300,cex=0.7)
dev.off()
plotMDS(lcpm, labels =as.character(edge_count$samples$condition))

############################create design and contrast
edge_count$samples$Hepatocyte = factor(edge_count$samples$Hepatocyte)
edge_count$samples$condition = factor(edge_count$samples$condition)

#save(edge_count,file="../results/edgeR_count_15575_429_TMM.rda")
edge_count = get(load("../results/edgeR_count_15575_429_TMM.rda"))
save(edge_count,file="../results/edgeR_count_15575_429_TMM.rda")
design = model.matrix(~0+condition+Hepatocyte,data=edge_count$samples)
contrast = c()
for(i in 2:7){
  for(j in 1:(i-1)){
    contrast = c(contrast,paste0("condition",i,"-condition",j))
  }
}
print(contrast)
contr.matrix = makeContrasts(contrasts=contrast, 
                             levels = colnames(design))

#######################linear modeling
v <- voom(edge_count, design, plot = TRUE)

vfit = lmFit(v,design)
vfit2 = contrasts.fit(vfit,contrasts = contr.matrix)
vfit3 <- eBayes(vfit2)

save(vfit3,file="../results/LIMMA_edgeR_DE_15575_429.rda")

#######################
vdata = voom(edge_count)
save(vdata,file="../results/LIMMA_edgeR_voom_norm_15575_429.rda")




