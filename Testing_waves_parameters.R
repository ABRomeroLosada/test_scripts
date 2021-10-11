library(circacompare)

gene.expression <- read.table(file="../gene_expression.tsv",header = T,sep = "\t")
gene.expression.SD <- gene.expression[,31:48]

gene.expression <- read.table(file="../gene_expression.tsv",header = T,sep = "\t")
gene.expression.LD <- gene.expression[,1:18]

swath.normalized.data.SD <- read.table(file = "../sd_swath_processed_data.tsv",header=T,sep="\t")
rownames(swath.normalized.data.SD)<- swath.normalized.data.SD$X
swath.normalized.data.SD$X <- NULL

swath.normalized.data.LD <- read.table(file = "../swath_processed_data.tsv",header=T,sep="\t")
swath.normalized.data.LD[,"zt20_2"] <- swath.normalized.data.LD[,"zt16_2"]
swath.normalized.data.LD[,"zt16_2"] <- swath.normalized.data.LD[,"zt20_2"]
