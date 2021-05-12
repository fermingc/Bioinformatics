>HW1 <- "C:/Users/fermi/Documents/Fall 2020/HW1/renal_cell_carcinoma.txt"
>HW1_anno <- "C:/Users/fermi/Documents/Fall 2020/HW1/renal_carcinoma_annotation.txt"

#1)*********
>data<-read.table(HW1, header=T)
>anno<-read.table(HW1_anno)

> dim(data)

#2)**********
>sub.data<-subset(data, select= anno[,1])

#>while(i<22){
names(sub.data)[names(sub.data) == colnames(sub.data[i])] <- colnames(sub.data[i])+"-"+anno[i,9]
i<-i+1
}

#setnames(sub.data, colnames(sub.data[1]) , colnames(sub.data[1])+"_"+anno[1,9])

#names(sub.data)[names(sub.data) == "oldVariableName"] <- "newVariableName"

#> colnames(sub.data[1])

> names(sub.data)[names(sub.data) == "GSM146778"] <- "GSM146778_Normal"
> names(sub.data)[names(sub.data) == "GSM146780"] <- "GSM146780_Normal"
> names(sub.data)[names(sub.data) == "GSM146782"] <- "GSM146782_Normal"
> names(sub.data)[names(sub.data) == "GSM146784"] <- "GSM146784_Normal"
> names(sub.data)[names(sub.data) == "GSM146786"] <- "GSM146786_Normal"
> names(sub.data)[names(sub.data) == "GSM146789"] <- "GSM146789_Normal"
> names(sub.data)[names(sub.data) == "GSM146790"] <- "GSM146790_Normal"
> names(sub.data)[names(sub.data) == "GSM146792"] <- "GSM146792_Normal"
> names(sub.data)[names(sub.data) == "GSM146794"] <- "GSM146794_Normal"
> names(sub.data)[names(sub.data) == "GSM146798"] <- "GSM146798_Normal"
> names(sub.data)[names(sub.data) == "GSM146796"] <- "GSM146796_Normal"
> View(sub.data)
> names(sub.data)[names(sub.data) == "GSM146779"] <- "GSM146779_Tumor"
> names(sub.data)[names(sub.data) == "GSM146781"] <- "GSM146781_Tumor"
> names(sub.data)[names(sub.data) == "GSM146783"] <- "GSM146783_Tumor"
> names(sub.data)[names(sub.data) == "GSM146785"] <- "GSM146785_Tumor"
> names(sub.data)[names(sub.data) == "GSM146787"] <- "GSM146787_Tumor"
> names(sub.data)[names(sub.data) == "GSM146788"] <- "GSM146788_Tumor"
> names(sub.data)[names(sub.data) == "GSM146791"] <- "GSM146791_Tumor"
> names(sub.data)[names(sub.data) == "GSM146799"] <- "GSM146799_Tumor"
> names(sub.data)[names(sub.data) == "GSM146793"] <- "GSM146793_Tumor"
> names(sub.data)[names(sub.data) == "GSM146795"] <- "GSM146795_Tumor"
> names(sub.data)[names(sub.data) == "GSM146797"] <- "GSM146797_Tumor"


#colnames(sub.data) <- paste(colnames(sub.data), "tot_proc", sep = "_")
#3)*********ID outliers using following plots

#4)*********

#Correlation Plot
library(gplots)
dat.cor <- cor(sub.data,use="pairwise.complete.obs")

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"blue","yellow","white"))
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)
image(dat.cor,main="Stage 1/2 cRCC",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)

par(mar=c(1,1,1,1)) 
image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)

#Hierarchical clustering dendrogram
dat<-t(sub.data)
dat.dist<-dist(dat,method="euclidean")
dat.clust<-hclust(dat.dist,method="single")
plot(dat.clust,labels=names(dat),cex=.75)

#CV vs. mean plot
dat.mean <- apply(log2(sub.data),2,mean)		# calculate mean for each sample
dat.sd <- sqrt(apply(log2(sub.data),2,var))		# calculate st.deviation for each sample
dat.cv <- dat.sd/dat.mean			#calculate cv

plot(dat.mean,dat.cv,main="Stage 1/2 cRCC dataset Sample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)

sort(dat.mean) #unable to read crowded labels

# average correlation plot
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))

png(filename="figure.png", width=900, bg="white") #y labels were off page without this
par(mar=c(5,6,4,1)+.1)
barplot(c(1.1, 0.8, 0.7), horiz=TRUE, border="blue", axes=FALSE, col="darkblue")
axis(2, at=1:3, lab=c("elephant", "hippo", "snorkel"), las=1, cex.axis=1.3)
dev.off()

plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of Tumor/Normal samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(sub.data)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

#5)******

#6)*******

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")

#7)**** Remove the outlier samples you identified in the first part of this assignment.

sub.data<-sub.data[,-10] #GSM146798
sub.data<-sub.data[,-18] #GSM146799


#8)******Files Downloaded from: http://biogps.org/#goto=genereport&id=359
KNG1 <- "C:/Users/fermi/Documents/Fall 2020/HW1/KNG1.txt"
AQ2 <- "C:/Users/fermi/Documents/Fall 2020/HW1/AQ2.txt"

KNG1Data<-read.table(KNG1, header=T)
AQ2Data<-read.table(AQ2, header=T)

#First KNG1 Proboset - X206054_at 
#X= KNG1Data[,1]; Y = KNG1Data[,2]
plot(as.numeric(KNG1Data[,2]),type='l',lwd=2,col='blue', main="KNG1 Expression, Proboset X206054_at" ,xlab="Cell Type", ylab="Intensity",axes=F)
axis(1,at=c(1:length(KNG1Data[,1])),labels=as.vector(KNG1Data[,1]),las=2,cex.axis=0.7)
axis(2)

#Second KNG1 Proboset - X217512_at
#X= KNG1Data[,1]; Y = KNG1Data[,3]
plot(as.numeric(KNG1Data[,3]),type='l',lwd=2,col='blue', main="KNG1 Expression, Proboset X217512_at" ,xlab="Cell Type", ylab="Intensity",axes=F)
axis(1,at=c(1:length(KNG1Data[,1])),labels=as.vector(KNG1Data[,1]),las=2,cex.axis=0.7)
axis(2)

#AQ2 Proboset X206672_at
#X= AQ2Data[,1]; Y= AQ2Data[,2]
plot(as.numeric(AQ2Data[,2]),type='l',lwd=2,col='blue', main="AQ2 Expression, Proboset X206672_at" ,xlab="Cell Type", ylab="Intensity",axes=F)
axis(1,at=c(1:length(AQ2Data[,1])),labels=as.vector(AQ2Data[,1]),las=2,cex.axis=0.7)
axis(2)

> AQ2Data[
+     order(AQ2Data[,2] ),
+ ]
 
> KNG1Data[
+     order(KNG1Data[,2] ),
+ ]

> KNG1Data[
+     order(KNG1Data[,3] ),
+ ]
 




#9)******
nine.val<-sub.data['206054_at','GSM146784_Normal']
sub.data['206054_at','GSM146784_Normal']ablinea<-NA

>dataMatrix<-data.matrix(sub.data, rownames.force = NA)

#10)*****
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")

library(impute)

knn<-impute.knn(dataMatrix, k=6, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=362436069)

> knnData<-knn[["data"]]


#11)*****
#relative error = (actual value - predicted)/actual

> knnData['206054_at','GSM146784_Normal']#value = 7559.533


#relative error = (actual value - predicted)/actual
>(nine.val-knnData['206054_at','GSM146784_Normal'])/nine.val #0.09847789


#12)*****
>if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
>BiocManager::install("pcaMethods")
>library(pcaMethods)

> result<-pca(dataMatrix, method="svdImpute", nPcs=9)
> SVDData<-completeObs(result)

> SVDData['206054_at','GSM146784_Normal']


#13)*****
#x = dataMatrix[0,]
#y = dataMatrix['206054_at',] --> Normal data, blue
#  = knnData['206054_at',] --> Knn imputation, green
#  = > SVDData['206054_at',] --> SVD imputation, red

> dataMatrix['206054_at','GSM146784_Normal']<-nine.val

> plot(as.numeric(knnData['206054_at',]),lwd=2,col='green', main="KNG1 Expression, Proboset X206054_at with Imputed Values" ,xlab="Patient Sample", ylab="Intensity",axes=F)
axis(1,at=c(1:length(colnames(dataMatrix))),labels=as.vector(colnames(dataMatrix)),las=2,cex.axis=0.7)
axis(2)

points(SVDData['206054_at',], col= 'red')
points(dataMatrix['206054_at',], col='blue',pch=19)


legend(1, 3000, legend=c("KNN", "SVD", "Actual"),
       col=c("green", "red", "blue"), lty=1:2, cex=0.8)
