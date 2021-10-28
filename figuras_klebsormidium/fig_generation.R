gene.expression <- read.table(file = "gene_expression_kleb.tsv", header=T, as.is=T)

gene.expression.HL <- (gene.expression[,"HL_1"] + gene.expression[,"HL_2"])/2
names(gene.expression.HL) <- gene.expression$geneID

mean.expression.HL <- (gene.expression[,"HL_1"] + gene.expression[,"HL_2"])/2
names(mean.expression.HL) <- gene.expression$geneID

mean.expression.LL <- (gene.expression[,"LL_1"] + gene.expression[,"LL_2"])/2
names(mean.expression.LL) <- gene.expression$geneID

heatmap.enzyme.id <- function(gene.id,enzyme.name,gene.expression,precision)
{
  pos.colfunc <- colorRampPalette(c("white","gold"))
  pos.colors <- pos.colfunc(precision)
  
  neg.colfunc <- colorRampPalette(c("white","lightslateblue"))
  neg.colors <- neg.colfunc(precision)
  
  fc.max <- 0
  fc.min <- 20
  
  enzyme.mean.HL <- mean.expression.HL[gene.id]
  enzyme.mean.LL <- mean.expression.LL[gene.id]
  fc <- enzyme.mean.HL/enzyme.mean.LL
  
  if(fc >= 1)
  {
    if(fc >= fc.max)
    {
      fc.max <- fc
    }
    
    if(round(fc*10) < precision)
    {
      cols <- pos.colors[round(fc*10)]
    } else
    {
      cols <- pos.colors[precision]
    }
    
  } else
  {
    if(fc <= fc.min)
    {
      fc.min <- fc
    }
    
    if(round(10/fc) < precision)
    {
      cols <- neg.colors[round(10/fc)]
    } else
    {
      cols <- neg.colors[precision]
    }
    
  }
  plot<-plot(x=0,y=0,col="white",
             axes=F,xlab="",ylab="",
             ylim=c(0,10),xlim=c(0,10))
  
  plot<-polygon(x = c(2,8,8,2),y=c(2,2,8,8),lwd=6,col = cols)
  return(plot)
}

gene.id<- "kfl00001_0850"
gene.name<- "kfl00001_0850"
png(filename =paste(c("heatmap_",gene.name,".png"),collapse=""))
heatmap.enzyme.id(gene.id=gene.id, 
                  enzyme.name = "unknow",  
                  gene.expression=gene.expression,
                  precision=50)
dev.off()


## Carotenoids
carotenoids <- read.table(file="carotenoids_klebsormidium.tsv",header=T,sep="\t",as.is=T)
head(carotenoids)
carotenoids.names <- carotenoids$carotenoid

carotenoids <- carotenoids[,2:ncol(carotenoids)]
rownames(carotenoids) <- carotenoids.names
head(carotenoids)

carotenoids.p.value <- vector(mode = "numeric",length=length(carotenoids.names))
carotenoids.fold.change <- vector(mode = "numeric",length=length(carotenoids.names))

names(carotenoids.p.value) <- carotenoids.names
names(carotenoids.fold.change) <- carotenoids.names

for(i in 1:length(carotenoids.names))
{
  current.carotenoid <- carotenoids.names[i]
  
  hl.1 <- unlist(carotenoids[current.carotenoid,1:4])
  ll.1 <- unlist(carotenoids[current.carotenoid,5:7])
  mean.ll.1 <- mean(ll.1)
  norm.hl.1 <- hl.1/mean.ll.1
  norm.ll.1 <- ll.1/mean.ll.1
  
  norm.ll <- c(norm.ll.1)
  norm.hl <- c(norm.hl.1)
  
  means <- c(mean(norm.ll), mean(norm.hl))
  sds <- c(sd(norm.ll), sd(norm.hl))
  
  if(means[1] > means[2])
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "greater")
  } else
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "less")
  }
  
  carotenoids.p.value[i] <- sig$p.value
  carotenoids.fold.change[i] <- means[2]/means[1]
  
  png(filename = paste0(current.carotenoid,".png"),width = 300)
  par(lwd=3)
  xpos <- barplot(means,col=c("lightslateblue","gold"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=current.carotenoid,
                  cex.main=2)
  points(x = rep(xpos[1]+0.2,3),y=norm.ll)
  points(x = rep(xpos[2]+0.2,4),y=norm.hl)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  dev.off()
}


## Gradient
png(file="scale_klebs")
colfunc<-colorRampPalette(c("lightslateblue","white","gold"))
plot(rep(1,1000),col=(colfunc(1000)), pch=15,cex=20,xlim=c(0,730),axes=F,
     ylab="",xlab="")
dev.off()
