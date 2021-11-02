library(GAMBLR)
library(dplyr)
library(ComplexHeatmap)
library(matrixStats)
library(maftools)
library(ggpubr)
library(dendextend)
library(tidyverse)

set.seed(123)
setwd("~/GAMBLR")

sample.names.table <-read.table("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/2021.10.04.all.RNAseqs.txt", stringsAsFactors = FALSE)
sample.names <- strsplit(sample.names.table$V1, split="/")
sample.names <- unlist(lapply(sample.names, tail, 1L))

sample.names.new <- strsplit(sample.names, split="[.]")
sample.names.final <- sapply(sample.names.new, "[", 1 )
num.samples <- length(sample.names.final)

splice.types<-c("SE","RI","A3","A5")

concat.splice<-function(event){
  file.name <- paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,".final.csv")
  splice.final <- read.table(file.name, sep = ",",check.names = FALSE)
  splice.final$avg <- rowMeans(splice.final[,11:1337],na.rm=TRUE)
  splice.final <- splice.final %>% subset(avg > 20)
  
    if ( event == "SE") {
      splice.final$loc<-paste(splice.final$geneSymbol, "SE", splice.final$upstreamEE, splice.final$exonStart_0base, splice.final$exonEnd, splice.final$downstreamES,sep="_")
    } else if ( event == "RI") {
      splice.final$loc<-paste(splice.final$geneSymbol, "RI", splice.final$upstreamEE, splice.final$riExonStart_0base, splice.final$riExonEnd, splice.final$downstreamES,sep="_")
    } else {
      splice.final$loc<-paste(splice.final$geneSymbol, event, splice.final$longExonStart_0base, splice.final$longExonEnd, splice.final$shortES, splice.final$shortEE,sep="_")
    }
  
  splice.final <- splice.final[,2665:ncol(splice.final)]
  splice.final$var <- rowVars(as.matrix(splice.final[,1:num.samples]), na.rm=TRUE)
  
  splice.final <- splice.final %>% subset(var > 0.01)
  splice.loc<-splice.final$loc
  splice.colnames<-colnames(splice.final)
  
  splice.final <- matrix(as.numeric(unlist(splice.final[,1:num.samples])), ncol=(num.samples))
  
  rownames(splice.final)<-splice.loc
  colnames(splice.final)<-splice.colnames[1:1327]
  
  file.name.2<-paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,"final.concat.csv")
  write.table(splice.final, file = file.name.2,sep=",", quote = FALSE)
}

for(splice in splice.types){
  concat.splice(splice)
}



subset.splice<-function(event){
  file.name.2<-paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,"final.concat.csv")
  splice.final <- read.table(file.name.2,sep=",", check.names=FALSE, row.names=NULL)
  col.variance<-colVars(as.matrix(splice.final[,2:(num.samples+1)]), na.rm = TRUE)
  col.subset<-(col.variance > 0)
  col.subset<-replace_na(col.subset,FALSE)
  
  splice.final <- splice.final[,col.subset]
  
  
  rnaseqs.meta<-get_gambl_metadata(seq_type_filter = "mrna")
  samples.splice <- colnames(splice.final)
  samples.splice <- str_remove_all(samples.splice, ("_Inc" ))
  colnames(splice.final) <- samples.splice
  
  splice.meta <- rnaseqs.meta %>% subset(sample_id %in% samples.splice)
  subset.splice <- splice.meta$sample_id
  
  splice.matrix.subset <-  (colnames(splice.final) %in% subset.splice)
  splice.final.subset <-splice.final[,splice.matrix.subset]
  splice.final.subset$row.names<-splice.final$row.names
  file.name.3<-paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,"final.subset.csv")
  write.table(splice.final.subset, file.name.3, quote=FALSE,sep = ",")
}


for (splice in splice.types){
  subset.splice(splice)
}


rnaseqs.meta<-get_gambl_metadata(seq_type_filter = "mrna")
event<-"A3"
cluster.splice.types <- function(event){
  file.name.2 <- paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,"final.subset.csv")
  splice.matrix <- read.table(file.name.2, sep = ",",header=TRUE, check.names = FALSE, row.names=NULL)
  matrix.rows<-splice.matrix[,ncol(splice.matrix)]
  matrix.cols <- colnames(splice.matrix)
  splice.matrix <- splice.matrix %>%
    t()
  colnames(splice.matrix)<-matrix.rows
  rownames(splice.matrix)<-matrix.cols
  splice.matrix<-splice.matrix[-c(1,nrow(splice.matrix)),]
  
  # row.variance<-rowVars(hclust.matrix)
  # row.subset<-(row.variance > 0)
  # row.subset<-replace_na(row.subset,FALSE)
  # 
  # final.matrix <- hclust.matrix[row.subset,]
  splice.matrix <- replace(splice.matrix,is.na(splice.matrix),0)
  
  dend <- splice.matrix %>% dist() %>% hclust(method = "complete") %>% as.dendrogram()
  samples<-rownames(splice.matrix)
  meta.subset<- rnaseq.meta %>% subset(sample_id %in% samples)
  
    colors_to_use <- as.numeric(as.factor(meta.subset[,10]))
    colors_to_use <- colors_to_use[order.dendrogram(dend)]
    
    labels_colors(dend)<-colors_to_use
  # The default `plot()` function can be used to produce a simple dendrogram
  file.name.3 <- paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,"final.subset.cluster.png")
  png(filename=file.name.3, width = 30, height = 12, units = "in", res=300, type="cairo")
  plot(dend, main = paste0("Sample clustering based on ", event))
  
  
  dev.off()
  
  
}

for(splice in splice.types){
  cluster.splice.types(splice)
}

splice.heatmaps<-function(event){
  file.name.3 <- paste0("/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/rmats-4.0.1_grch37/level3/",event,"final.subset.csv")
  splice.final <- read.table(file.name.3, sep = ",", check.names = FALSE)
  rownames(splice.final) <- splice.final$row.names
  sample.names <- colnames(splice.final)[1:num.samples]
  rnaseqs.subset <- rnaseq.meta %>% subset(Tumor_Sample_Barcode %in% sample.names)
  column_ha <- HeatmapAnnotation(pathology = rnaseqs.subset$pathology, dhit_sig = rnaseqs.subset$DHITsig_consensus)
  heatmap<-Heatmap(splice.final, cluster_columns = TRUE, cluster_rows = FALSE, bottom_annotation = column_ha)
  draw(heatmap) 
}
