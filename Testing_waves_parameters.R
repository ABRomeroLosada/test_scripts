library(circacompare)

gene.expression <- read.table(file="gene_expression.tsv",header = T,sep = "\t")
row.names(gene.expression) <- gene.expression$X
gene.expression$X <- NULL
gene.expression.SD <- gene.expression[,43:84]
gene.expression.SD.DD <- gene.expression[,73:84]
gene.expression.SD.LL <- gene.expression[,61:72]

gene.expression.LD <- gene.expression[,1:42]
gene.expression.LD.DD <- gene.expression[,31:42]
gene.expression.LD.LL <- gene.expression[,19:30]

gene <- "ostta01g00880"

time.points <- seq(from=0,by=4,length.out = 12)
circacompare.LL.DD <- matrix(ncol=1,nrow=15)

SD.DD <-gene.expression.SD.DD[gene,]
SD.LL <-gene.expression.SD.LL[gene,]

circacomp.data <- data.frame(time=c(time.points,time.points),
                             measure=c(t(SD.DD/max(SD.DD)),
                                       t(SD.LL/max(SD.LL))),
                             group=c(rep("DD",12),rep("LL",12)))
result.i<- circacompare(x = circacomp.data, 
                        col_time = "time", 
                        col_group = "group", 
                        col_outcome = "measure",
                        alpha_threshold = 1)
png(file="circacompareDDLLostta01g00880.png")
result.i$plot
dev.off()

t <- seq(from=0, to=20, by=0.5)
wave.DD= 5.487292e-01 + 3.977647e-01 * cos( 24 * (t-8.005159))
wave.LL= 7.339324e-01 + 1.429324e-01 * cos( 24 * (t-2.345114e+01))

png(file="DDnew.png")
plot(x= t, 
     y=wave.DD, 
     type="o",col="red",
     lwd=3, ylab="", xlab="", axes=T, lty=3)
# lines(x = time.points, type="o",
#       y=SD.DD/max(SD.DD),
#       lwd=4,col="red",lty=1)
# axis(side=1,lwd=3,
#      at =time.points, labels=colnames(SD.DD),
#      las=2)
# legend("topright", 
#        legend = c("cicacompare DD","real DD"), 
#        lwd=2, col = c("red", "blue", "orange"), 
#        lty=c(3,1,1))
dev.off()

png(file="LLnew.png")
plot(x= t, 
     y=wave.LL, 
     type="o",col="red",
     lwd=3, ylab="", xlab="", axes=T, lty=3)
# lines(x = time.points, type="o",
#       y=SD.DD/max(SD.DD),
#       lwd=4,col="red",lty=1)
# axis(side=1,lwd=3,
#      at =time.points, labels=colnames(SD.DD),
#      las=2)
# legend("topright", 
#        legend = c("cicacompare DD","real DD"), 
#        lwd=2, col = c("red", "blue", "orange"), 
#        lty=c(3,1,1))
dev.off()



