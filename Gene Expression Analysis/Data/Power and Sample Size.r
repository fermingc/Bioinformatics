>Lab3File <- "C:/Users/fermi/Documents/Fall 2020/Lab 3/eisen.txt"

#2)****
dat<-read.table(Lab3File, header=T, na.strings="NA", blank.lines.skip=F, row.names=1)

#3)****
eClasses<- "C:/Users/fermi/Documents/Fall 2020/Lab 3/eisenClasses.txt"

eCdat<- read.table(eClasses, header=T)

#4)****
> sub.dat<-subset(dat, select = eCdat[,2])
#1-19 classes = 1, 20-39 classes -

#5)****
#picking gene from row 1
> dat.m <- mean(as.numeric(as.matrix(sub.dat[1,])), na.rm=T)
> dat.m

> sub.dat[1,is.na(sub.dat[1,])] <- dat.m   

gc<-sub.dat[,1:19]
act<-sub.dat[,20:39]

> x<-as.numeric(gc[1,])
> y<-as.numeric(act[1,])

#boxplot
> xy.list<-list(x,y)
> boxplot(xy.list,col=c('blue','red'),main='Example Gene from DLBCL cDNA 2-channel dataset',axex=F,ylab="log2(ratio intensity)")
> axis(2)
> axis(1,at=c(1,2),c("GC","ACT"))

#histogram
> par(mfrow=c(2,1))
>hist(x)
>hist(y)

#6)***
> nx<-length(x)
> ny<-length(y)

> pool.var<-(((nx-1)*var(x))+((ny-1)*var(y)))/(nx+ny-2)
>pool.var

> install.packages("pwr")
>library(pwr)

>dif.1.5fold<-log2(1.5)/sqrt(pool.var)
> pl.ss1.5<-pwr.t.test(d=dif.1.5fold,sig.level=.01,power=0.8,type="two.sample")
> pl.ss1.5


#7)***
> dif<-abs(mean(x)-mean(y))/sqrt(pool.var)
> pl.ss<-pwr.t.test(d=dif,sig.level=.01,power=0.8,type="two.sample")
> pl.ss

#8)***
>library(ssize)
>library(gdata)
> install.packages("matrixStats")
> library("matrixStats")

>data(exp.sd)

> dat.matrix<-data.matrix(sub.dat)
> SDs<-rowSds(dat.matrix, na.rm=TRUE)
> hist(SDs,col="cyan",border="blue",main="",xlab="Standard Deviation (for data on the log2 scale)")
> title("Histogram of Standard Deviations for 13,412 Genes")

#9)****
>fold.change=3
>power=.8
>sig.level=.05

> all.size<-ssize(sd=SDs, delta=log2(fold.change),sig.level=sig.level, power=power)
> ssize.plot(all.size,lwd=2,col="magenta",xlim=c(1,20))
> xmax<-par("usr")[2]-1;
> ymin<-par("usr")[3]+.05
> legend(x=xmax,y=ymin,legend=strsplit( paste("fold change=",fold.change,",", "alpha=", sig.level,",", "power=",power,",", "# genes=", length(SDs), sep=''),",")[[1]],xjust=1,yjust=0,cex=1.0)
> title("Sample Size to Detect 2-Fold Change")







