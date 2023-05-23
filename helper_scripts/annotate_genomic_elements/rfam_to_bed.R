args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
species <- args[2]
infile <- read.delim(input, header = F)
ncRNA <- cbind(infile[,c(4,1:3,6)],V7 = species)
colnames(ncRNA) <- c("GeneID","Chr","Start","End","Strand","Species")
ncRNA$Classify <- "other_ncRNA"
ncRNA$Classify[grep("rRNA",ncRNA$GeneID)] <- "rRNA"
ncRNA$Classify[grep("tRNA",ncRNA$GeneID)] <- "tRNA"
ncRNA$Classify[grep("^U\\d|SNOR|ACEA|^DdR\\d",ncRNA$GeneID)] <- "snRNA"
for (classify in unique(ncRNA[ncRNA$Species == species,]$Classify)) {
  write.table(
    x = ncRNA[ncRNA$Species==species&ncRNA$Classify==classify,c(2:4,1,7,5)], 
    file = paste0("../",species,"/annot_", classify, ".bed"), 
    sep = "\t", row.names = F, quote = F, col.names = F)
}