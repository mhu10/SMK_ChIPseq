library(dplyr)

out_file1 <- snakemake@output[["bed_ctl"]]
out_file2 <- snakemake@output[["bed_treatment"]]

datainput_ctl <- read.table(snakemake@input[['peaks_ctl']], sep="\t")[,c(2,3,4,1,6,5)]
datainput_treatment <- read.table(snakemake@input[['peaks_treatment']], sep="\t")[,c(2,3,4,1,6,5)]


colnames(datainput_ctl) = c('Chromosome','Start','End','PeakID','Column5','Strand')
colnames(datainput_treatment) = c('Chromosome','Start','End','PeakID','Column5','Strand')


datainput_ctl = arrange(datainput_ctl, Chromosome)
datainput_treatment = arrange(datainput_treatment, Chromosome)


write.table(datainput_ctl, out_file1, quote = FALSE,sep='\t',row.names = FALSE)
write.table(datainput_treatment, out_file2, quote = FALSE,sep='\t',row.names = FALSE)


