library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(vsn)
library(edgeR)
library(pheatmap)
library(limma)
library(biomaRt)
setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
drug = c('Baseline', 'Omeprazole','Phenobarbital','Dexamethasone','Carbamazepin','Phenytoin','Rifampicin')

vfit3 = get(load("../results/LIMMA_edgeR_DE.rda"))
dt = decideTests(vfit3)
dt = decideTests(vfit3, p.value = 0.01)
table = summary(dt)

###############pair-wise DE count and plot

contrast = c()
for(i in 2:7){
  for(j in 1:(i-1)){
    contrast = c(contrast,paste0("condition",i,"-condition",j))
  }
}
print(contrast)

out = matrix( rep( 0, len=49), nrow = 7)
all(contrast == colnames(dt))
for(i in 1:ncol(table)){
  down = as.numeric(table[,i][1])
  up = as.numeric(table[,i][3])
  alt = as.numeric(substr(contrast[i],10,10))
  ref = as.numeric(substr(contrast[i],21,21))
  out[alt,ref] = down
  out[ref,alt] = up
}

colnames(out) = drug
rownames(out) = drug
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pdf("../figure/DE_gene_clustering_15575_429_0.01.pdf",onefile=T)
pheatmap(out,display_numbers = T, number_format = "%d",cluster_rows = F,cluster_cols = F)
dev.off()

##########################FDR is directly calculated from p value
index = as.numeric(which(sapply(colnames(dt), function(x) substr(x,21,21)) == "1"))
print(paste0("average down regulated gene " ,mean(table[,index][1,])))
print(paste0("average up regulated gene ", mean(table[,index][3,])))

########################################write table for all up-regulated and all-down regulated genes

upindex = which(apply(dt[,index],1, function(x) all(x==1)))
downindex =  which(apply(dt[,index],1, function(x) all(x==-1)))
allup = vfit3$genes[upindex,]
alldown = vfit3$genes[downindex,]

p_test = vfit3$p.value[,index[1]]
fdr_test = p.adjust(p_test,method = "BH")
length(which(dt[,index[1]] != 0))
all(which(fdr_test<0.05) == which(dt[,index[1]] != 0))
FDR = apply(vfit3$p.value,2,function(x) p.adjust(x,method = "BH"))

table = NULL
for(i in 1:length(index)){
  temp_col = paste0(contrast[index[i]],c("_lfc","_p","_FDR"))
  temp =  cbind(vfit3$coefficients[,index[i]],vfit3$p.value[,index[i]], FDR[,index[i]])
  colnames(temp) = temp_col
  table = cbind(table,temp)
}

table_up = table[upindex,]
table_down = table[downindex,]

table_up = cbind(allup,table_up)
table_down = cbind(alldown,table_down)

col = colnames(table_up)
for(i in 1:7){
  col = gsub(condit[i],drug[i],col)
}

colnames(table_up) = col
colnames(table_down) = col
write.table(table_up,"../results/LIMMA_all_UP_15575_429.txt",sep = "\t",quote = F,col.names = T,row.names = F)
write.table(table_down,"../results/LIMMA_all_DOWN_15575_429.txt",sep = "\t",quote = F,col.names = T,row.names = F)
write.table(vfit3$genes,"../results/LIMMA_background_15575_429.txt", sep = "\t", quote=F, col.names = F,row.names = F)
########################


