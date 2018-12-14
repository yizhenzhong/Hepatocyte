# Download the data from github (click the "raw" button, save as a text file called "results.txt").
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")

drug = c('Baseline', 'Omeprazole','Phenobarbital','Dexamethasone','Carbamazepin','Phenytoin','Rifampicin')

vfit3 = get(load("../results/LIMMA_edgeR_DE_15579_429.rda"))
index = as.numeric(which(sapply(colnames(vfit3$coefficients), function(x) substr(x,21,21)) == "1"))
vfit_b = vfit3[, index]
dt = decideTests(vfit3, p.value = 0.01)[, index] 

#topTable(vfit_b[,1])
#results <- classifyTestsF(vfit3[,index])
#vennCounts(results)

#the coefficient is the log fold change
#confired from vfit_b$coefficients[grep("ENSG00000140465", rownames(vfit_b$coefficients)),]
#topTable(vfit_b[,1])[1:3,]
library(ggrepel)
library(ggplot2)
source("../Yilin/script/Baseline/helper.R")
plotVolcano <- function(i){
  table = data.frame("logFC"=vfit_b$coefficients[,i], 
                     "P.Value"=vfit_b$p.value[,i],
                     "gene" = vfit_b$genes$HGNC.symbol)
  table$label = as.numeric(dt[,i])
  table$label[table$label != 0] <- "DE"
  table$label[table$label == 0] <- "All"
  
  color = RColorBrewer::brewer.pal(9, "Set1")
  table$label = factor(table$label)
  table$logp = -log10(table$P.Value)
  table = table[order(table$P.Value),]
  ggplot(table,aes(x=logFC, y=logp)) +
    geom_point(aes(col = label), size=2.5)+
    theme_hepa()+
    labs(y=expression("-log10(p)"))+
    labs(x=expression("logFC"))+
    scale_color_manual(values=c(color[9], color[1]))+
    ggtitle(drug[i+1])+
    theme(plot.title = element_text(size=20))+
    guides(color=FALSE)+
    geom_text_repel(data=head(table, 20), aes(label=gene))
  
}



p1 = plotVolcano(1)
p2 = plotVolcano(2)
p3 = plotVolcano(3)
p4 = plotVolcano(4)
p5 = plotVolcano(5)
p6 = plotVolcano(6)

library(gridExtra)
combine = grid.arrange(p1, p2, p3, p4, p5, p6,
             nrow = 2)

ggsave("../figure/Volcano_plot.png", height = 10, width = 15, plot = combine)

