## Load and preprocess metabolomic data from first run
raw.metabolomic.data.1 <- read.table(file="klebso_metabolomic_data_1.tsv",header=T,sep="\t",as.is=T)
head(raw.metabolomic.data.1)
metabolomic.data.1 <- raw.metabolomic.data.1[,3:ncol(raw.metabolomic.data.1)]
weight.1 <- raw.metabolomic.data.1$Peso

for(i in 1:nrow(metabolomic.data.1))
{
  print(i)
  metabolomic.data.1[i,] <- (metabolomic.data.1[i,]/weight.1[i])
}

## Load and preprocess metabolomic data from second run
raw.metabolomic.data.2 <- read.table(file="klebso_metabolomic_data_2.tsv",header=T,sep="\t",as.is=T)
metabolomic.data.2 <- raw.metabolomic.data.2[,2:ncol(raw.metabolomic.data.2)]
paracetamol <- metabolomic.data.2$Paracetamol
weight.2 <- c(20.6, 17.8, 18.2, 19.6, 19.2, 17.5)

for(i in 1:nrow(metabolomic.data.2))
{
  print(i)
  metabolomic.data.2[i,] <- (metabolomic.data.2[i,]/weight.2[i])/paracetamol[i]
}


colnames(metabolomic.data.1)
colnames(metabolomic.data.2)

metabolites <- intersect(colnames(metabolomic.data.1), colnames(metabolomic.data.2))
setdiff(colnames(metabolomic.data.1), colnames(metabolomic.data.2))
setdiff(colnames(metabolomic.data.2), colnames(metabolomic.data.1))

i <- 1

metabolites.p.value <- vector(mode = "numeric",length=length(metabolites))
metabolites.fold.change <- vector(mode = "numeric",length=length(metabolites))

names(metabolites.p.value) <- metabolites
names(metabolites.fold.change) <- metabolites

normalize.metabolites <- matrix(nrow = length(metabolites),ncol = 12)
rownames(normalize.metabolites) <- metabolites
colnames(normalize.metabolites) <- c(paste("HL",1:6,sep="_"),paste("LL",1:6,sep="_"))

for(i in 1:length(metabolites))
{
  current.metabolite <- metabolites[i]
  
  hl.1 <- metabolomic.data.1[1:3,current.metabolite]
  ll.1 <- metabolomic.data.1[4:6,current.metabolite]
  mean.ll.1 <- mean(ll.1)
  norm.hl.1 <- hl.1/mean.ll.1
  norm.ll.1 <- ll.1/mean.ll.1
  
  hl.2 <- metabolomic.data.2[1:3,current.metabolite]
  ll.2 <- metabolomic.data.2[4:6,current.metabolite]
  mean.ll.2 <- mean(ll.2)
  norm.hl.2 <- hl.2/mean.ll.2
  norm.ll.2 <- ll.2/mean.ll.2
  
  norm.ll <- c(norm.ll.1, norm.ll.2)
  norm.hl <- c(norm.hl.1, norm.hl.2)
  
  normalize.metabolites[i,] <- c(norm.hl,norm.ll)
  
  means <- c(mean(norm.ll), mean(norm.hl))
  sds <- c(sd(norm.ll), sd(norm.hl))
  
  if(means[1] > means[2])
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "greater")
  } else
  {
    sig <- t.test(x = norm.ll,norm.hl,alternative = "less")
  }
  metabolites.p.value[i] <- sig$p.value
  
  metabolites.fold.change[i] <- means[2]/means[1]
  
  png(filename = paste0(current.metabolite,".png"),width = 300)
  par(lwd=3)
  xpos <- barplot(means,col=c("blue","firebrick2"),
                  names.arg = c("LL","HL"),las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=current.metabolite,
                  cex.main=2)
  points(x = rep(xpos[1]+0.2,6),y=norm.ll)
  points(x = rep(xpos[2]+0.2,6),y=norm.hl)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  dev.off()
}

head(normalize.metabolites)

activated.metabolites <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change > 1]
repressed.metabolites <- metabolites[metabolites.p.value < 0.05 & metabolites.fold.change < 1]

sort(metabolites.fold.change)

plot(x = log2(metabolites.fold.change),y = -log10(metabolites.p.value), 
     pch=19,xlim = c(-2,2),xlab="log2(Fold-change)",ylab="-log10(p-value)",cex.lab=1.5)

lines(x = c(-4,4),y=c(-log10(0.05),-log10(0.05)),lty=2)
points(x = log2(metabolites.fold.change)[activated.metabolites],
       -log10(metabolites.p.value)[activated.metabolites], col="red",pch=19)
points(x = log2(metabolites.fold.change)[repressed.metabolites],
       -log10(metabolites.p.value)[repressed.metabolites], col="blue",pch=19)

text(x = log2(metabolites.fold.change)["Tryptophan"],
     y = -log10(metabolites.p.value)["Tryptophan"]-0.15,labels = "Tryptophan")

# text(x = log2(metabolites.fold.change)["a.ketoglutaric"],
#      y = -log10(metabolites.p.value)["a.ketoglutaric"]-0.15,
#      labels = expression(paste(alpha,"-ketoglutaric")))

text(x = log2(metabolites.fold.change)["Histidine"],
     y = -log10(metabolites.p.value)["Histidine"]-0.15,labels = "Histidine")

text(x = log2(metabolites.fold.change)["Phenylalanine"],
     y = -log10(metabolites.p.value)["Phenylalanine"]-0.15,labels = "Phenylalanine")

# text(x = log2(metabolites.fold.change)["AMP"],
#      y = -log10(metabolites.p.value)["AMP"]-0.15,labels = "AMP")

# text(x = log2(metabolites.fold.change)["ADPG"],
#      y = -log10(metabolites.p.value)["ADPG"]-0.15,labels = "ADPG")

# text(x = log2(metabolites.fold.change)["Succinic"],
#      y = -log10(metabolites.p.value)["Succinic"]-0.15,labels = "Succinic")


text(x = log2(metabolites.fold.change)["Glutamine"],
     y = -log10(metabolites.p.value)["Glutamine"]-0.15,labels = "Glutamine")

text(x = log2(metabolites.fold.change)["X6phosphogluconic"],
     y = -log10(metabolites.p.value)["X6phosphogluconic"]-0.15,labels = "6-Phosphogluconic")

activated.metabolites
repressed.metabolites
sort(metabolites.fold.change,decreasing = F)

## Carotenoids
carotenoids <- read.table(file="../carotenoids/carotenoids_klebsormidium.tsv",header=T,sep="\t",as.is=T)
head(carotenoids)
carotenoids.names <- carotenoids$carotenoid

carotenoids <- carotenoids[,2:ncol(carotenoids)]
rownames(carotenoids) <- carotenoids.names
head(carotenoids)

carotenoids.p.value <- vector(mode = "numeric",length=length(carotenoids.names))
carotenoids.fold.change <- vector(mode = "numeric",length=length(carotenoids.names))

names(carotenoids.p.value) <- carotenoids.names
names(carotenoids.fold.change) <- carotenoids.names
i <- 1
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
  xpos <- barplot(means,col=c("blue","firebrick2"),
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
i
current.carotenoid

-log10(carotenoids.p.value)
log2(carotenoids.fold.change)


plot(x = c(log2(metabolites.fold.change),log2(carotenoids.fold.change)),
     y = c(-log10(metabolites.p.value),-log10(carotenoids.p.value)), 
     pch=19,xlim = c(-3,3),xlab="log2(Fold-change)",ylab="-log10(p-value)",cex.lab=1.5,ylim=c(0.12,3.5))

lines(x = c(-4,4),y=c(-log10(0.05),-log10(0.05)),lty=2)
points(x = log2(metabolites.fold.change)[activated.metabolites],
       -log10(metabolites.p.value)[activated.metabolites], col="red",pch=19)
points(x = log2(metabolites.fold.change)[repressed.metabolites],
       -log10(metabolites.p.value)[repressed.metabolites], col="blue",pch=19)













library(FactoMineR)
library(factoextra)

pca.metabolites <- data.frame(colnames(normalize.metabolites[,-1]),t(normalize.metabolites[,-1]))
colnames(pca.metabolites)[1] <- "Sample"

res.pca.metabolites <- PCA(pca.metabolites, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.metabolites <- HCPC(res.pca.metabolites, graph=FALSE)   
fviz_dend(res.hcpc.metabolites,k=2,
          cex = 1,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 10      # Augment the room for labels
)

fviz_pca_ind(res.pca.metabolites, col.ind = c(rep("HL",5), rep("LL",6)), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE) 




carotenoids.p.value <- vector(mode = "numeric",length=length(carotenoids.names))
carotenoids.fold.change <- vector(mode = "numeric",length=length(carotenoids.names))

names(carotenoids.p.value) <- carotenoids.names
names(carotenoids.fold.change) <- carotenoids.names
i <- 1
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
  xpos <- barplot(means,col=c("blue","firebrick2"),
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



metabolites.dist <- dist(t(normalize.metabolites),method = "euclidean")
library("cluster")

hierarchical.clustering <- hclust(as.dist(metabolites.dist),method="average")
plot(hierarchical.clustering)
hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)




metabolites.dist <- dist(t(normalize.metabolites[,c("HL_1","HL_4","HL_5","HL_6","LL_2","LL_3","LL_4","LL_6")]),method = "euclidean")
library("cluster")

hierarchical.clustering <- hclust(as.dist(metabolites.dist),method="average")
plot(hierarchical.clustering)

