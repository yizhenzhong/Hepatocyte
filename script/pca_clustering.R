setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
setwd("/Volumes/fsmresfiles/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
#count = read.table("../data/liver_RNAseq_count1.txt",sep = "\t",stringsAsFactors = F,header = T)

sample = NULL
exp = NULL
for(i in c(1:7)){
  fn = paste0("../data/expression/condition1_tpm.txt")
  f =  read.table(fn,stringsAsFactors = F,sep = "\t",header = T)
  print(paste(i,"condition, expression dimension:", dim(f)[1],  dim(f)[2],sep=" "))
  temp = as.data.frame(colnames(f)[-1])
  temp$condition = i
  colnames(temp)[1] = "sample"
  sample = rbind(sample,temp)
  if(is.null(exp)){exp = f}else{exp = merge(exp,f,by="gene")}
}

sample$sample = as.factor(sample$sample)
sample$condition = as.factor(sample$condition)
###################

plot_pca<-function(input_table){
  voom_counts <- limma::voom(input_table)$E
  pca_voom <- prcomp(t(voom_counts),center = T,scale. = T)
  eigs <- pca_voom$sdev^2
  percentVar = c(round(100*eigs[1] / sum(eigs)),round(100*eigs[2] / sum(eigs)))
  pca_out = pca_voom$x
  data_pca_plot <- data.frame(PC1 = pca_out[,1],
                              PC2 = pca_out[,2],
                              PC3 = pca_out[,3],
                              PC4 = pca_out[,4],
                              PC5 = pca_out[,5],
                              PC6 = pca_out[,6],
                              PC7 = pca_out[,7],
                              PC8 = pca_out[,8],
                              PC9 = pca_out[,9],
                              PC10 = pca_out[,10],
                              Treatment = sample$condition)
  data_pca_plot = list(table = data_pca_plot,pc=percentVar)
  return(data_pca_plot)
}


plot_pca_table = plot_pca(exp[,-1])
#####################################



dds_filter_ro =dds_filter[,-which(data_pca_plot$PC1< (-200))]
plot_pca_table = plot_pca(dds_filter_ro)


tiff("../figure/pca_exp_tmm_voom_condit1_6.tiff",height = 17, width = 17, 
     units = "cm", compression = "lzw", res = 1000)

ggplot(plot_pca_table$table,aes(x=PC1,y=PC2,fill = Treatment))+geom_point(colour = "grey20",pch = 21, alpha = 1,size=2)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())+
  theme(axis.text = element_text(size=20),axis.title = element_text(size = 25),
        legend.key.size = unit(1.5,'lines'),legend.key.width = unit(3,"cm"),legend.text =  element_text(size = 20),
        legend.box="horizontal",legend.title = element_text(colour="black", size=14, face="bold"),
        plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),legend.key = element_rect(fill = "transparent", 
                                                                                      colour = "transparent"))+
  xlab(paste0("PC1: ", plot_pca_table$pc[1], "% variance"))+ylab(paste0("PC2: ", plot_pca_table$pc[2], "% variance"))
#+xlab(paste0("PC1: ", percentVar[1], "% variance"))+ylab(paste0("PC2: ", percentVar[2], "% variance"))

dev.off()
######################################
#hierarchical clustering
#####################################

euclidean_mat <- as.matrix(dist(t(exp[,-1]), method="euclidean"))
hc =hclust(euclidean_mat)
condition = paste0("condition",1:7)
#plot(hc,labels=FALSE, xlab = NULL)
marker1 = list(color =brewer.pal(7,"Dark2"))$color
species_col_adjusted  = marker1[sample$condition][hc$order]

dent = as.dendrogram(hc)
dent <- dent %>% set("leaves_pch", c(15)) %>% set("leaves_cex", c(2)) %>% set("leaves_col", species_col_adjusted)
#dent <- hang.dendrogram(dent,hang_height=0.1)
dent <- set(dent, "labels_cex", 0.5)


tiff("../figure/hierarchical_clustering_17949_416_condition.tiff",height = 17, width = 17, 
     units = "cm", compression = "lzw", res = 1000)

par(mar = c(3,3,3,3))
circlize_dendrogram(dent,labels = F)
#plot(dent,
#     main = "Clustered Hepatocyte",
#     horiz =  F,  nodePar = list(cex = .007),labels=F)


#legend("topleft", legend = condition, fill =marker1)
dev.off()


t(exp[,-1])


###################

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)
species_col_adjusted  = col[sample$sample][hc$order]
dent = as.dendrogram(hc)
dent <- dent %>% set("leaves_pch", c(15)) %>% set("leaves_cex", c(2)) %>% set("leaves_col", species_col_adjusted)
#dent <- hang.dendrogram(dent,hang_height=0.1)
dent <- set(dent, "labels_cex", 0.5)


tiff("../figure/hierarchical_clustering_17949_416_by_sample.tiff",height = 17, width = 17, 
     units = "cm", compression = "lzw", res = 1000)

par(mar = c(3,3,3,3))
#plot(dent,
 #    main = "Clustered Hepatocyte",
#     horiz =  F,  nodePar = list(cex = .007),labels=F)
#legend("topleft", legend = condition, fill =marker1)
circlize_dendrogram(dent,labels = F)
dev.off()

############################################################

library(CountClust)
exp_gom = FitGoM(as.matrix(t(assay(dds_filter))), K=4, tol=0.1)
data("MouseDeng2014.FitGoM")

circlize_dendrogram(dent,labels = F)
