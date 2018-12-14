library(data.table)
#install_github('kkdey/maptpx') 
library(maptpx)
library(CountClust)
library(data.table)
library(argparser, quietly=TRUE)

p <- add_argument(p, "--k", help="", default=0.001)
argv <- parse_args(p)

drug = c('Baseline', 'Omeprazole','Phenobarbital','Dexamethasone','Carbamazepin','Phenytoin','Rifampicin')

edge_count = get(load(file="../results/edgeR_count_15575_429.rda"))
edge_count = edge_count[1:10,]

###############topic modeling
Topic_Clus=topics(t(edge_count$counts),argv$k,kill=0,tol=0.1);
save(Topic_Clus,file=paste0("../results/countclust_",argv$k,".rda"))

#################plot heatmap
index = 10
Topic_Clus = get(load(paste0("../results/countclust_", index, ".rda")))
docweights=Topic_Clus$omega;
#drug[edge_count$samples$condition]
samples = drug[edge_count$samples$condition]

docweights_per_tissue_mean <- apply(docweights, 2,
                                    function(x) { tapply(x, samples, mean) })

pdf(paste0("../figure/heatmap_countclust",index,".pdf"))
#pdf(paste0("../figure/heatmap_countclust_6.pdf"))
par(mar=c(5,6,4,6)+.1)
ordering <- heatmap(1-docweights_per_tissue_mean , margins=c(4,10))$rowInd
dev.off()

#############structure plot
# order tissue by hierarhical clustering results
tissue_levels_reordered <- rownames(docweights_per_tissue_mean)[ordering]
annotation <- data.frame(
  sample_id = paste0("X", 1:nrow(docweights)),
  tissue_label = factor(samples,
                        levels = rev(tissue_levels_reordered ) ) )

cols1 <- c(rev(RColorBrewer::brewer.pal(9, "Set1")))

pdf("../figure/structureplot_6cluster.pdf")
StructureGGplot(omega = docweights,
                annotation= annotation,
                palette = cols1,
                yaxis_label = "",
                order_sample = TRUE,
                #yaxis_label = "Drug argument",
                split_line = list(split_lwd = .1,
                                  split_col = "white"),
                legend_text_size = 15,
                legend_title_size = 15, 
                legend_key_size = 0.7,
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 14,
                                 axis_label_face="bold"))
dev.off()

######################extract top expressed genes
topics_theta <- Topic_Clus$theta;

top_features <- ExtractTopFeatures(topics_theta, top_features=100, method="poisson", options="min");
gene_names <- edge_count$genes
gene_names <- substring(gene_names,1,15);
xli  <-  gene_names;
gene_list <- do.call(cbind, lapply(1:dim(top_features$indices)[1], 
                                   function(x) gene_names$HGNC.symbol[top_features$indices[x,]]))

colnames(gene_list) = paste0("cluster",1:6)
write.table(gene_list, "../results/count6clust_top_feature.txt",row.names = F,col.names = T, quote = F,
            sep = "\t")



