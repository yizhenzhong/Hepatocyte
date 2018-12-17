
setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")
#save(vdata,file="../results/LIMMA_edgeR_voom_norm_15575_429.rda")
library(variancePartition)
library('doParallel')
library(ggrepel)
source("../Yilin/script/Baseline/helper.R")

#The variancePartition package implements a
#statistical method to quantify the contribution of multiple sources of variation
#and decouple within/between-individual variation.


vdata = get(load("../results/edgeR_count_15579_429_voom_v.rda"))
#Counts mapping to each gene can be normalized using counts per
#million (CPM), reads per kilobase per million (RPKM) or fragments
#per kilobase per million (FPKM). These count results can be processed
#with limma/voom [6] to model the precision of each observation
#or DESeq2 [7].



cl <- makeCluster(4)
registerDoParallel(cl)
vdata$targets$Sex = factor(vdata$targets$Sex)
vdata$targets$Platform = factor(vdata$targets$Platform)
vdata$targets$Batch = factor(vdata$targets$Batch)

form <- ~ (1|Hepatocyte) + (1|condition) + Age + (1|Sex) + Ancestry + (1|Platform) + (1|Batch)
#Categorical variables should (almost) always be modeled as a random effect.
#The difference between modeling a categorical variable as a fixed versus random
#effect is minimal when the sample size is large compared to the number of
#categories (i.e. levels). So variables like disease status, sex or time point will
#not be sensitive to modeling as a fixed versus random effect. However, variables
#with many categories like Individual must be modeled as a random effect in
#order to obtain statistically valid results. So to be on the safe side, categorical
#variable should be modeled as a random effect.

varPart <- fitExtractVarPartModel(vdata$E, form, vdata$targets)
vp <- sortCols( varPart )

plotVarPart( vp)+ 
  theme_hepa()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "None")
ggsave("../figure/variancePartition.png")

save(varPart,file="../results/variancePartition_15579_429samples.rda")


###################
table = cbind(vdata$genes, varPart$condition,varPart$Hepatocyte,varPart$Age, varPart$Sex, varPart$Ancestry)
colnames(table)[(ncol(table)-4): ncol(table)] = c("drug","hepatocyte","age","sex","Ancestry")

table = as.data.frame(table)
table$drug =  as.numeric(as.character(table$drug))
table$hepatocyte =  as.numeric(as.character(table$hepatocyte))
table = table[order(table$drug, decreasing = T),]
#write.table(table, "../results/variancePartition_Summary_table_15579.txt", quote = F, row.names = F, sep = "\t")

cols1 <- c(rev(RColorBrewer::brewer.pal(8, "Dark2")))
ggplot(table,aes(x=hepatocyte,y=drug,color=Gene.type)) + 
  geom_point(size=2) +
  theme_hepa()+  
  scale_color_manual(values=c(cols1[7],cols1[8]))+
  labs(y="variance across drug treatments")+
  labs(x="variance across hepatocytes")+
  theme(legend.position = c(0.75, 0.8),
        legend.text=element_text(size=18))+
  geom_text_repel(data=head(table, 20), aes(label=HGNC.symbol))

ggsave("../figure/variacne_partition.png")

############################
library(gProfileR)
#var = read.table("../results/variancePartition_Summary_table.txt", header = T, stringsAsFactors = F,

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



