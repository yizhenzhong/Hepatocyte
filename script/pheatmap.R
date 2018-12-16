#the log2 counts per million after voom normalization 
library(sjstats)
library(RColorBrewer)

v=get(load("../results/edgeR_count_15579_429_voom_v.rda"))
cvs = apply(v$E, 1, cv)
gene_cv = v$E[order(cvs, decreasing = T)[1:1000],]



drug = c('Baseline', 'Omeprazole','Phenobarbital','Dexamethasone','Carbamazepin','Phenytoin','Rifampicin')

annotation_col = data.frame(
  drug= factor(as.character(drug[v$targets$condition]))
)

rownames(annotation_col) = colnames(gene_cv)


colors = brewer.pal(7, "Dark2")
ann_colors = list(
  drug = c(Baseline=colors[1],
           Carbamazepin=colors[2],
           Dexamethasone=colors[3],
           Omeprazole=colors[4],
           Phenobarbital=colors[5],
           Phenytoin=colors[6],
           Rifampicin=colors[7])
)
pheatmap(gene_cv, 
         scale = "row", #scale by row
         clustering_distance_rows = "correlation", 
         filename = "../figure/pheatmap2.png",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         show_rownames = F, show_colnames = F,
         height = 14, width = 14,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         fontsize = 15)



#############################
ann_colors = list(
  drug = c(Baseline=colors[1],
           Omeprazole=colors[2]))
annotation_col = data.frame(
  drug= factor(as.character(drug[v$targets$condition][index]))
)
colors = brewer.pal(7, "Dark2")
rownames(annotation_col) = colnames(gene_cv)[index]
index = v$targets$condition %in% c(1,2)
pheatmap(gene_cv[,index], 
         scale = "row", 
         clustering_distance_rows = "correlation", 
         filename = "../figure/pheatmap_1and2.png",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         show_rownames = F, show_colnames = F,
         height = 14, width = 14,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         fontsize = 15)
