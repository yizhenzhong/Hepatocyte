setwd("/Volumes/fsmresfiles/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/script/")


fam = read.csv("../data/factors_60_GSK.csv",stringsAsFactors = F,header = T)
pc = read.table("../data/genotype/Hepatocyte_merge1000genome_PCs.txt",header = F,
                stringsAsFactors = F)
pc$shortid = sapply(pc$V2,function(x) tail(strsplit(x,"-")[[1]],1))
fam$shortid = sapply(fam$ID,function(x) tail(strsplit(x,"-")[[1]],1))
#fam$shortid[grep(pattern = "761",fam$V2)] = "76"
#fam$shortid[grep(pattern = "751",fam$V2)] = "75"
pc$shortid[grep(pattern = "761",pc$V2)] = "76"
pc$shortid[grep(pattern = "751",pc$V2)] = "75"



condit = 1
peer = read.table("../data/expression/PEER/condition1_expression_tmm_normalization.PEER_covariates.txt",
                  header = T,stringsAsFactors = F,check.names = F)
peer_shortid = sapply(colnames(peer),function(x) tail(strsplit(x,".",fixed = T)[[1]],1))


#finalsample = read.table("../data/genotype/test.dosage",header = T,stringsAsFactors = F,check.names = F)
#finalid = sapply(colnames(finalsample),function(x) tail(strsplit(x,"-")[[1]],1))

for(npeer in 3:15){
  covariate = t(peer[1:npeer,c(match(fam$shortid,peer_shortid))])
  colnames(covariate) = paste0("InferredCov",1:npeer)
  fam_sub = as.data.frame(cbind(fam$Sex,fam$Platform, fam$Batch))
  colnames(fam_sub)  = c("Sex", 'Platform', 'Batch') 
  covariate =  cbind(covariate,fam_sub)
  out_nopc = t(covariate)
  
  out_nopc = cbind(rownames(out_nopc), out_nopc)
  colnames(out_nopc) = c('ID', paste0("MP-",fam$shortid))
  write.table(out_nopc,paste0("../data/genotype/eQTL/condition",condit,"_",npeer,"peer_covariate_nopc.txt"),col.names = T,
              row.names = F,quote = F)
  
  pc3 =  as.matrix(pc[match(fam$shortid,pc$shortid),c(3:5)])
  colnames(pc3) = paste0("PC",1:3)
  covariate =  cbind(covariate,pc3)
  
  out = t(covariate)
  out = cbind(rownames(out), out)
  colnames(out) = c('ID', paste0("MP-",fam$shortid))
  write.table(out,paste0("../data/genotype/eQTL/condition",condit,"_",npeer,"peer_covariate.txt"),col.names = T,
              row.names = F,quote = F)
  
}
