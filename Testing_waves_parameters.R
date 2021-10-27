library(circacompare)

gene.expression <- read.table(file="gene_expression.tsv",header = T,sep = "\t")
row.names(gene.expression) <- gene.expression$X
gene.expression$X <- NULL
gene.expression.SD <- gene.expression[,43:84]
gene.expression.SD.DD <- gene.expression[,73:84]
gene.expression.SD.LL <- gene.expression[,61:72]
names(gene.expression.SD.LL)

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


time.points <- seq(from=0,by=4,length.out = 12)
circacompare.LL.LL <- matrix(nrow=length(genes),ncol=15)
for (i in 1:length(genes))
{
gene <- genes[i]
SD.LL <-gene.expression.SD.LL[gene,]
LD.LL <-gene.expression.LD.LL[gene,]

circacomp.data <- data.frame(time=c(time.points,time.points),
                             measure=c(t(LD.LL/max(LD.LL)),
                                       t(SD.LL/max(SD.LL))),
                             group=c(rep("LL_LD",12),rep("LL_SD",12)))
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
boxplot(circacompare.LL.LL[,13])
mean(circacompare.LL.LL[,13])


circacompare.DD.DD <- matrix(nrow=length(genes),ncol=15)

for (i in 1:length(genes))
{
   gene <- genes[i]
   SD.DD <-gene.expression.SD.DD[gene,]
   LD.DD <-gene.expression.LD.DD[gene,]
   
   circacomp.data <- data.frame(time=c(time.points,time.points),
                                measure=c(t(LD.DD/max(LD.DD)),
                                          t(SD.DD/max(SD.DD))),
                                group=c(rep("LL_LD",12),rep("LL_SD",12)))
   result.i<- circacompare(x = circacomp.data, 
                           col_time = "time", 
                           col_group = "group", 
                           col_outcome = "measure",
                           alpha_threshold = 1)
   circacompare.DD.DD[i,] <- result.i$summary[,2]
}
rownames(circacompare.DD.DD) <- as.character(genes)
colnames(circacompare.DD.DD) <- result.i[[2]][,1]

hist(circacompare.DD.DD[,13])
mean(circacompare.DD.DD[,13])

######Datos de marcos
#Function for extracting the strenght and frequency of the waves that compose the resulting one
#as calculated by GeneCycle
#circ.data <- as.numeric(fousdll)
extract.individual.waves <- function(circ.data) #n.samples must be even
{
   perio <- periodogram(circ.data)
   waves <- matrix(nrow = length(circ.data)/2, ncol = 2)
   for (i in 1:length(circ.data)/2)
   {
      waves[i,1] <- perio$freq[i]*length(circ.data)
      waves[i,2] <- perio$spec[i]
   }
   return(waves)
}

# Function for extracting the trajectory from fft output
get.trajectory <- function(X.k,ts,acq.freq) {
   
   N   <- length(ts)
   i   <- complex(real = 0, imaginary = 1)
   x.n <- rep(0,N)           # create vector to keep the trajectory
   ks  <- 0:(length(X.k)-1)
   
   for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
      x.n[n+1] <- sum(X.k * exp(i*2*pi*ks*n/N)) / N
   }
   
   x.n * acq.freq 
}

#Function to draw each of the simple waves that make up the complex waveform
plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
   Xk.h <- rep(0,length(Xk))
   Xk.h[i+1] <- Xk[i+1] # i-th harmonic
   harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
   points(ts, harmonic.trajectory, type="l", col=color)
}

# Function for drawing the periodogram of a complex wave to determine the contribution
#of each wave to the complex pattern
draw.periodogram <- function(circ.data, harmonics=6)
{
   period.data <- periodogram(circ.data)
   plot(period.data$freq[1:harmonics]*length(circ.data), 
        period.data$spec[1:harmonics]/sum(period.data$spec), 
        xlab="Harmonics (Hz)", ylab="Amplitute Density", type="h")
}

###Más funciones
wave.decomposition.scaled <- function(data.fou, n.waves,time=48, acq.freq=0.25, ts=seq(0,44,4))
{
   waves <- extract.individual.waves(data.fou)
   colnames(waves) <- c("freq", "strength")
   data.fft <- fft(data.fou)
   names(data.fft) <- NULL
   wave.order <- order(waves[,2],decreasing = T)
   x.n <- get.trajectory(data.fft,ts,acq.freq)
   x.n <- x.n/max(Mod(x.n))
   amp <- max(Mod(x.n)) - min(Mod(x.n))
   plot(ts,x.n,type="l",lwd=2, ylim=c(0,1))
   wave.colors <- c("red", "blue", "green", "purple", "orange", "gray")
   fft.vector <- rep(0,length(data.fft))
   for (j in n.waves)
   {
      #j<-1
      fft.vector[wave.order[j]+1] <- data.fft[wave.order[j]+1]
      fft.vector[length(ts)+2-(wave.order[j]+1)] <- data.fft[length(ts)+2-(wave.order[j]+1)]
   }
   fft.vector[1] <- data.fft[1]
   harmonic.trajectory <- get.trajectory(fft.vector, ts, acq.freq=acq.freq)
   harmonic.trajectory <- harmonic.trajectory/max(Mod(harmonic.trajectory))
   points(ts, harmonic.trajectory, type="l",lty=2, col="darkred", lwd=3)
   
   return(list(wave.params=waves, wave.points=data.fft,freqs.order=wave.order,comb.wave=harmonic.trajectory))
}
photo.and.skoto <- function(gene.name, genes.expression, sd.ll, sd.cycle)
{
   fou <- genes.expression[gene.name,]
   fousdll =fou[sd.ll]
   wave.order <- wave.decomposition.scaled(as.numeric(fousdll), n.waves = 1:3, time=24, ts=seq(0,20,4))$freqs.order
   order.pos <- which(wave.order == 1) #Extracts the wave with an associated period of 24h
   dec.wave.ll <- wave.decomposition.scaled(as.numeric(fousdll), n.waves = order.pos, time=24, ts=seq(0,20,4))$comb.wave
   
   foutot <- fou[sd.cycle]
   dec.wave.ll.2 <- wave.decomposition.scaled(as.numeric(foutot), n.waves = 1:6, time=48, ts=seq(0,44,4))$comb.wave
   dec.wave.res.cycle <- dec.wave.ll.2 - dec.wave.ll
   ext.dec.wave.cycle.order <- wave.decomposition.scaled(as.numeric(dec.wave.res.cycle), n.waves = 1:6)$freqs.order
   ord.pos <- which(ext.dec.wave.cycle.order == 2) #Extracts the wave with an associated period of 24h
   ext.dec.wave.cycle <- wave.decomposition.scaled(as.numeric(dec.wave.res.cycle), n.waves = ord.pos)$comb.wave
   ext.dec.wave.cycle <- ext.dec.wave.cycle + 1 #Makes all points positive
   ext.dec.wave.cycle <- ext.dec.wave.cycle/max(Mod(ext.dec.wave.cycle))
   dec.wave.ll <- dec.wave.ll - min(Mod(dec.wave.ll))
   dec.wave.ll <- dec.wave.ll/max(Mod(dec.wave.ll))
   dec.wave.ll.2 <- dec.wave.ll.2 - min(Mod(dec.wave.ll.2))
   dec.wave.ll.2 <- dec.wave.ll.2/max(Mod(dec.wave.ll.2))
   plot(seq(0,44,4), ext.dec.wave.cycle, type="l",lty=2, col="darkgrey", lwd=3, ylim=c(-0.2,1),
        main=gene.name, ylab = "Normalized amplitude", xlab="Time")
   points(seq(0,44,4), rep(dec.wave.ll,2), type="l",lty=2, col="darkorange", lwd=3)
   points(seq(0,44,4), dec.wave.ll.2, type="l", col="black", lwd=3) 
   legend("bottomright", legend = c("photoperiod", "skotoperiod"), fill = c("darkorange","darkgrey"))
   return(list(norm.wave = dec.wave.ll.2, norm.photo = dec.wave.ll, skoto.norm = ext.dec.wave.cycle))
}

photo.and.skoto(gene.name=gene, genes.expression = gene.expression,
                sd.ll= names(gene.expression.SD.LL)[1:6],
                sd.cycle = names(gene.expression.SD)[1:12])
