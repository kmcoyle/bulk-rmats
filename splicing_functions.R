get_splice_data <- function(sample_id, base_dir,  splice_type) {
  countsfile = paste0(base_dir, "/", sample_id, "/", splice_type, ".MATS.JC.txt")
  countdata <- read.table(countsfile, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if(splice_type %in% c("RI", "A5SS","A3SS", "SE")){
    colnames(countdata)[13] = paste0(sample_id, "_IJC_SAMPLE_1")
    colnames(countdata)[14] = paste0(sample_id, "_SJC_SAMPLE1")
    colnames(countdata)[21] = paste0(sample_id, "_IncLevel1")
    # Remove extraneous columns
    counts <- countdata[ ,c(2:11,13:14,21)]
    counts[,13]<-as.numeric(counts[,13])
    counts <- counts[, -c(11:12)]
  } else if(splice_type == "MXE") {
    colnames(countdata)[15] = paste0(sample_id, "_IJC_SAMPLE_1")
    colnames(countdata)[16] = paste0(sample_id, "_SJC_SAMPLE1")
    colnames(countdata)[23] = paste0(sample_id, "_IncLevel1")
    # Remove extraneous columns
    counts <- countdata[ ,c(2:13,15:16,23)]
    counts[,13]<-as.numeric(counts[,15])
    counts <- counts[, -c(15:16)]
  }
  return(counts)
}
