>FilePath <- "C:/Users/fermi/Documents/Fall 2020/Lab 9/lung_cancer.txt"

#1)****
dat<- read.table(FilePath, header=T, row.names=1)

#2)***
> class=c("adeno","adeno","adeno","adeno","adeno","adeno","adeno","adeno","adeno","adeno","SCLC","SCLC","SCLC","SCLC","SCLC","SCLC","SCLC","SCLC","SCLC","Normal","Normal","Normal","Normal","Normal")
> clas<-names(dat)
> datx<-as.data.frame(t(dat))

datx<-data.frame(class,datx)

#3)***
traindat<-datx[1:6,]
traindat<-rbind(traindat, datx[11:16,])
traindat<-rbind(traindat, datx[20:22,])

testdat<-datx[7:10,]
testdat<-rbind(testdat, datx[17:19,])
testdat<-rbind(testdat, datx[23:24,])

> testclass<-testdat[,1]
> testdat<-testdat[,-1]

#4)***
> traindat.lda<-lda(clas~X1007_s_at + X1053_at,traindat)
> testdat.pred<-predict(traindat.lda, testdat[,1:2])

#5)**
#> testdat.pred$x[,1] vs > testdat.pred$x[,2]
#> c(rownames(testdat.pred$x))

> plot(testdat.pred$x[,1], testdat.pred$x[,2], main='Discriminant Functions', xlab='Discriminant Function 1', ylab='Discriminant Function 2', col=1:length(c(rownames(testdat.pred$x))), lwd=3)
> legend("topleft",legend=c(rownames(testdat.pred$x)), fill=1:length(c(rownames(testdat.pred$x))))


#6)***
> traindat.all.lda<-lda(clas~.,traindat)
> testdat.all.pred<-predict(traindat.all.lda, testdat)


#7)***
> plot(testdat.all.pred$x[,1], testdat.all.pred$x[,2], main='Discriminant Functions', xlab='Discriminant Function 1', ylab='Discriminant Function 2', col=1:length(c(rownames(testdat.all.pred$x))), lwd=3)
> legend(0,0,legend=c(rownames(testdat.all.pred$x)), fill=1:length(c(rownames(testdat.all.pred$x))))
