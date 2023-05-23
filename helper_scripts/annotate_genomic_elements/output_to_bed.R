args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]
infile <- read.delim(input, header = F)
infile$V8 <- '-'
infile$V8[infile$V2<infile$V3] <- '+'
rev <- infile[infile$V2>infile$V3,c(1,3,2,4,7,8)]
colnames(rev) <- colnames(infile[,c(1:4,7,8)])
outfile <- rbind(infile[infile$V2<infile$V3,c(1:4,7,8)],rev)
outfile$V7 <- outfile$V7/100
outfile <- outfile[order(outfile$V1,outfile$V3),]
rownames(outfile) <- seq(1:nrow(outfile))
newdf <- data.frame()
low <- outfile$V2[1]
high <- outfile$V3[1]
pad <- 300
for (i in seq(2,(nrow(outfile)-1))) {
  if (high+pad > outfile$V2[i] & outfile$V1[i-1] == outfile$V1[i]) {
    low <- min(low, outfile$V2[i])
    high <- max(high, outfile$V3[i])
  } else {
    newdf <- rbind(newdf, c(outfile$V1[i-1],low,high,outfile$V4[i-1],outfile$V7[i-1],outfile$V8[i-1]))
    low <- outfile$V2[i]
    high <- outfile$V3[i]
  }
}

colnames(newdf) <- c("chrom","start","end","ID","score","strand")
#newdf <- newdf[rep(seq_len(nrow(newdf)), each = 2), ]
#newdf$strand <- rep(c("+","-"),nrow(newdf)/2)

write.table(newdf, output, sep = "\t", row.names = F, col.names = F, quote = F)
