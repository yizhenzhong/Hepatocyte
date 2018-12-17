
setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
#save(vdata,file="../results/LIMMA_edgeR_voom_norm_15575_429.rda")
library(variancePartition)
#The variancePartition package implements a
#statistical method to quantify the contribution of multiple sources of variation
#and decouple within/between-individual variation.


vdata = get(load("../results/edgeR_count_15579_429_voom_v.rda"))
#Counts mapping to each gene can be normalized using counts per
#million (CPM), reads per kilobase per million (RPKM) or fragments
#per kilobase per million (FPKM). These count results can be processed
#with limma/voom [6] to model the precision of each observation
#or DESeq2 [7].

#vdata_rm = get(load("../results/edgeR_count_15575_429_voom_v_rmbatch.rda"))
#vdata$targets$cell = as.numeric(vdata$targets$Hepatocyte)
#vdata$targets$drug = as.numeric(vdata$targets$condition)
#write.csv(vdata$targets[,5:ncol(vdata$targets)],"../r  esults/LIMMA_edgeR_voom_norm_15575_429_cov.csv",quote = F)

cols1 <- c(rev(RColorBrewer::brewer.pal(9, "Set1")))

form <- ~ (1|Hepatocyte) + (1|condition) + Age + Sex + Ancestry + Platform + Batch

varPart <- fitExtractVarPartModel(vdata, form, vdata$targets)
save(varPart,file="../results/variancePartition_15575_429samples.rda")
save(varPart,file="../results/variancePartition_15575_429samples_rmbatch.rda")

varPart = get(load("../results/variancePartition_15575_429samples_rmbatch.rda"))

table = cbind(vdata$genes, varPart$condition,varPart$Hepatocyte,varPart$Age, varPart$Sex, varPart$Ancestry)

colnames(table)[(ncol(table)-4): ncol(table)] = c("drug","cell","age","sex","Ancestry")

table = as.data.frame(table)
table$drug =  as.numeric(as.character(table$drug))
table$cell =  as.numeric(as.character(table$cell))

write.table(table, "../results/variancePartition_Summary_table.txt", quote = F, row.names = F, sep = "\t")
pdf("../figure/variancePartition_rmbatch.pdf")


table$Gene.type[grep("^CYP",table$HGNC.symbol)] = "DME"
    
ggplot(table,aes(x=cell,y=drug,color=Gene.type)) + geom_point() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.title.x=element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())+
  theme(axis.text = element_text(size=20),axis.title = element_text(size = 25),legend.title = element_blank(),
        legend.key.size = unit(1.5,'lines'),legend.key.width = unit(3,"cm"),legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text =  element_text(size = 20),legend.position=c(0.75,0.75),legend.box="horizontal",
        plot.margin = unit(c(1,1,1,1), "cm"),axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))+
  scale_color_manual(values=c(cols1[3],cols1[1], cols1[2]))+labs(y="condition variance")+labs(x="cell variance")


dev.off()


############################
library(gProfileR)
var = read.table("../results/variancePartition_Summary_table.txt", header = T, stringsAsFactors = F,
                 sep = "\t")
var = var[order(var$cell, decreasing = T),]
hist(var$cell)
go_ana = gprofiler(var$HGNC.symbol[var$cell>0.9],ordered_query = T, 
          custom_bg = var$HGNC.symbol)

write.csv(go_ana, "../results/variancePartion_cell_biased_gene_goprofier.csv")

####drug
var = var[order(var$drug, decreasing = T),]
hist(var$drug)
go_ana = gprofiler(var$HGNC.symbol[var$drug>0.2],ordered_query = T, 
                   custom_bg = var$HGNC.symbol)

write.csv(go_ana, "../results/variancePartion_drug_biased_gene_goprofier.csv")

####sex

var = var[order(var$sex, decreasing = T),]
hist(var$sex)
go_ana = gprofiler(var$HGNC.symbol[var$sex>0.15],ordered_query = T, 
                   custom_bg = var$HGNC.symbol)

write.csv(go_ana, "../results/variancePartion_sex_biased_gene_goprofier.csv")


####ancestry
var = var[order(var$Ancestry, decreasing = T),]
hist(var$Ancestry)
go_ana = gprofiler(var$HGNC.symbol[var$Ancestry>0.1],ordered_query = T, 
                   custom_bg = var$HGNC.symbol)

write.csv(go_ana, "../results/variancePartion_anceatry_biased_gene_goprofier.csv")



