#1)***
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fibroEset")
	
> library(fibroEset)
data(fibroEset)

> dat<-exprs(fibroEset)
> anno<- fibroEset@phenoData@data[["species"]]
	
	
#2)****
>colnames(dat)<-anno
> set.seed(46)
randat<-dat[sample(nrow(dat), 50),]

#3)***
> randat.dist<-dist(randat, method="manhattan")
> dat.hca<-hclust(randat.dist, method="median")
> plot(dat.hca, main="fibroEset Cluster Dendotram")

> t.randat<-t(randat)
> t.randat.dist<-dist(t.randat, method="manhattan")
> t.dat.hca<-hclust(t.randat.dist, method="median")
> plot(t.dat.hca, main="fibroEset Cluster Dendotram")

#4)***
> heatmap.2(randat, main="HCA on fibroEset - 2D")

#5)***
> pc.randat<-prcomp(randat)
> pc2<-pc.randat[["rotation"]][,1:2]
> cl<-kmeans(pc2, centers=3)

