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

gene.id<- "kfl00001_0510"
gene.name<- "kfl00001_0510"
png(filename =paste(c("heatmap_",gene.name,".png"),collapse=""))
heatmap.enzyme.id(gene.id=gene.id, 
                  enzyme.name = "unknow",  
                  gene.expression=gene.expression,
                  precision=40)
dev.off()

## Gradient
colfunc<-colorRampPalette(c("lightslateblue","white","gold"))
plot(rep(1,1000),col=(colfunc(1000)), pch=15,cex=20,xlim=c(0,730),axes=F,
     ylab="",xlab="")
