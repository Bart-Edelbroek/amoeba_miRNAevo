args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
annot <- args[2]
tailLength <- as.numeric(args[3])
FTcutoff <- 0.3333 #allow only clusters with reads in same direction
Ccutoff <- 0.5 #max complexity of clusters (precision)
minLength <- 40
raw <- read.delim(file = input, header = T, sep = "\t")
raw$Strand[raw$FracTop>=0.5] <- '+'
raw$Strand[raw$FracTop<0.5] <- '-'
gfile <- read.delim(file = annot, sep = "\t", header = F, skip = 1)
clustSelect <- raw[#(raw$FracTop<FTcutoff|raw$FracTop>(1 - FTcutoff))&
                     raw$MajorRNAReads >= 10&
                     raw$Complexity<Ccutoff&
                     (raw$DicerCall!="N"|(nchar(raw$MajorRNA) >= 20 & nchar(raw$MajorRNA) <= 24))&
                     raw$Length>minLength,]
new_gff <- gfile[grep(paste(clustSelect$Name,collapse=";|"), 
      gfile$V9),]
new_gff$V7 <- clustSelect$Strand
new_gff$V4 <-  new_gff$V4- tailLength
new_gff$V4[new_gff$V4<1] <- 1
new_gff$V5 <- new_gff$V5+ tailLength
write.table(x = new_gff, file = args[4], sep = "\t", quote = F, row.names = F, col.names = F)