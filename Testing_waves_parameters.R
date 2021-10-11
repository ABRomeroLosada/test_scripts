library(circacompare)

gene.expression <- read.table(file="gene_expression.tsv",header = T,sep = "\t")
gene.expression.SD <- gene.expression[,44:85]
gene.expression.SD.DD <- gene.expression[,74:85]
gene.expression.SD.LL <- gene.expression[,62:73]

gene.expression.LD <- gene.expression[,1:43]
gene.expression.LD.DD <- gene.expression[,32:43]
gene.expression.LD.LL <- gene.expression[,20:31]



gene <- "ostta01g00880"
