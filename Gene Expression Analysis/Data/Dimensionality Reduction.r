>FilePath <- "C:/Users/fermi/Documents/Fall 2020/Lab 7/Sotiriou.txt"
>AnnoPath <- "C:/Users/fermi/Documents/Fall 2020/Lab 7/Sotiriou_annotations.txt"

#1)***
dat<- read.table(FilePath, header=T, row.names=1)
anno<-read.table(AnnoPath, header=T, row.names=1)

#2)**
#rename columns simply by site instead of sample?

#ann.dat2 = anno
#newnames = newcols

dat.pca<-prcomp(t(dat),cor=F)
newcols<-anno[,1]
colnames(dat)<- newcols

dat.loadings<-dat.pca$x[,1:3]

>plot(range(dat.loadings[,1]),range(dat.loadings[,2]),xlab='P1', ylab='P2',main='PCA plot of Sotiriou Data')
>points(dat.loadings[,1][anno=="KIU"],dat.loadings[,2][anno=="KIU"],col='blue', pch=16, cex=1.5)
>points(dat.loadings[,1][anno=="OXF"],dat.loadings[,2][anno=="OXF"],col='yellow', pch=16, cex=1.5)
>legend(10,-20,c('KIU', 'OXF'),col=c("blue", "yellow"), pch=15, cex=.7, horiz=F)

#3)***

>dat.pca.var<-round(dat.pca$sdev^2/sum(dat.pca$sdev^2)*100,2)
> plot(c(1:length(dat.pca.var)),dat.pca.var,type="b", xlab="# of components", ylab="% variance", pch=21, col=1, bg=3, cex=1.5)
> title("Scree Plot for % Variability Explained by Eigenvalue")

#4)***

#non-metric MDS
dat.dist<-dist(t(dat))
dat.mds<-isoMDS(dat.dist)
plot(dat.mds$points,type="n")
> points(dat.mds$points[,1][anno=="KIU"], dat.mds$points[,2][anno=="KIU"], col="blue", pch=16, cex=1.5)
> points(dat.mds$points[,1][anno=="OXF"], dat.mds$points[,2][anno=="OXF"], col="yellow", pch=16, cex=1.5)
> title(main="Non-metric MDS Plot of Sotiriou Data of Stress=20%")
> legend(25,-10,c("KIU", "OXF"), col=c('blue', 'yellow'), pch=15, cex=.7,horiz=F)

#clasical MDS
> dat.loc<-cmdscale(dat.dist)
> plot(dat.loc,type="n")
> points(dat.loc[,1][anno=="KIU"],dat.loc[,2][anno=="KIU"],col="blue",pch=16,cex=1.5)
> points(dat.loc[,1][anno=="OXF"],dat.loc[,2][anno=="OXF"],col="yellow",pch=16,cex=1.5)
> title(main="MDS Plot of Sotiriou Data")
> legend(-30,23,c("KIU","OXF"), col=c('blue', 'yellow'),pch=15,cex=.7,horiz=F)
