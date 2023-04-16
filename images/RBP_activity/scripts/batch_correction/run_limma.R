# load necessary packages
library(limma)

# file path
count.file <- '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.txt'
out.file <- '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.limma.txt'
list.path.batch <- list('batch.shRNA.K562' = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts/batch_correction/batch_shRNA_K562.txt',
                       'batch.shRNA.HepG2' = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts/batch_correction/batch_shRNA_HepG2.txt',
                       'batch.CRISPR.K562' = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts/batch_correction/batch_CRISPR_K562.txt',  
                       'batch.CRISPR.HepG2' = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts/batch_correction/batch_CRISPR_HepG2.txt')

# normalization
data <- read.table(count.file,header=T,sep='\t',quote='',row.names='AltAnalyze_ID')
total <- colSums(data)
sf <- total / 1e6
data.norm <- t(t(data)/sf)

# iterate each category, and combine
list.after.limma = list()
for (i in seq_along(list.path.batch)) {
    name <- names(list.path.batch)[i]
    path <- list.path.batch[[i]]
    batch <- read.table(path,header=T,sep='\t',quote='')
    data.norm.batch <- data.norm[,batch$file]
    batch.batch <- batch$batch
    result <- removeBatchEffect(data.norm.batch,batch.batch)
    list.after.limma[[name]] <- result
}
after.limma <- Reduce(cbind,list.after.limma)
after.limma <- round(after.limma,digits=0)

# clip, upscale, too slow to run, move to python
# clip <- function(x) {
#     pmax(x,0)
# }
# after.limma <- sapply(after.limma,clip)
# after.limma <- after.limma * 25

# write
write.table(after.limma,file=out.file,sep='\t',col.names = NA)

