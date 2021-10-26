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
gene <- "ostta03g03470"
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
png(file="circacompareDDLLostta03g03470.png")
result.i$plot
dev.off()

#probando con radianes
t <- seq(from=0, to=10, by=0.5)
wave.DD= 0.4766370219 + 0.4127642249 * cos( 4.190e11 * (t-23.8471404982))
wave.LL= 0.7611205567 + 0.1724573355 * cos( 4.190e11 * (t-23.7558093743))

# ##probando con h xq se supone que cicacompare trata con h
# t<-seq(from=0, to=20, by=0.83)
# wave.DD= 5.487292e-01 + 3.977647e-01 * cos( 24 * (t-8.005159))
# wave.LL= 7.339324e-01 + 1.429324e-01 * cos( 24 * (t-2.345114e+01))

wave.SD <- (wave.DD) + wave.LL
plot(x= t, 
     y=wave.DD/max(wave.DD), 
     type="o",col="red",
     lwd=3, ylab="", xlab="", axes=T, lty=3)
lines(x=t,
      y=wave.LL/max(wave.LL),
      type="o", col="orange",
      lwd=3, lty=3)

lines(x=t,
      y=wave.SD/max(wave.SD),
      type="o", col="red",
      lwd=3, lty=1)
legend("bottomleft", 
       legend = c("DD", "LL", "DD+LL"), 
       lwd=2, col = c("red", "orange", "red"), 
       lty=c(3,3,1))

######Comparing LL_SD vs. LL_LD and DD_SD vs. DD_LD
genes <- read.table(file="genes_two_peaks_sd_one_peak_ll_dd.txt", header = F)
genes <- genes$V1
i<-0
time.points <- seq(from=0,by=4,length.out = 12)
circacompare.LL.LL <- matrix(nrow=length(genes),ncol=15)

for (i in 1:length(genes))
{
gene <- genes[i]
SD.LL <-gene.expression.SD.LL[gene,]
LD.LL <-gene.expression.LD.LL[gene,]

circacomp.data <- data.frame(time=c(time.points,time.points),
                             measure=c(t(SD.LL/max(SD.LL)),
                                       t(LD.LL/max(LD.LL))),
                             group=c(rep("LL_SD",12),rep("LL_LD",12)))
result.i<- circacompare(x = circacomp.data, 
                        col_time = "time", 
                        col_group = "group", 
                        col_outcome = "measure",
                        alpha_threshold = 1)
circacompare.LL.LL[i,] <- result.i$summary[,2]
}
rownames(circacompare.LL.LL) <- as.character(genes)
colnames(circacompare.LL.LL) <- result.i[[2]][,1]
hist(circacompare.LL.LL[,13])
mean(circacompare.LL.LL[,13])


circacompare.DD.DD <- matrix(nrow=length(genes),ncol=15)

for (i in 1:length(genes))
{
   gene <- genes[i]
   SD.DD <-gene.expression.SD.DD[gene,]
   LD.DD <-gene.expression.LD.DD[gene,]
   
   circacomp.data <- data.frame(time=c(time.points,time.points),
                                measure=c(t(SD.DD/max(SD.DD)),
                                          t(LD.DD/max(LD.DD))),
                                group=c(rep("LL_SD",12),rep("LL_LD",12)))
   result.i<- circacompare(x = circacomp.data, 
                           col_time = "time", 
                           col_group = "group", 
                           col_outcome = "measure",
                           alpha_threshold = 1)
   circacompare.DD.DD[i,] <- result.i$summary[,2]
}
rownames(circacompare.DD.DD) <- as.character(genes)
colnames(circacompare.DD.DD) <- result.i[[2]][,1]
phase.diff <- matrix(circacompare.DD.DD[,13], circacompare.LL.LL[,13], ncol=2, nrow=length(genes))
colnames(phase.diff) <- c("DD_SD vs. DD_LD", "LL_SD vs. LL_LD")
boxplot(phase.diff, main="Phase difference between SD and LD")
mean(circacompare.DD.DD[,13])

