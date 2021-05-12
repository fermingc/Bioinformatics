#2)**
>cortexFile<- "C:/Users/fermi/Documents/Fall 2020/Lab 6/agingStudy11FCortexAffy.txt"
>cortexAnn<- "C:/Users/fermi/Documents/Fall 2020/Lab 6/agingStudy1FCortexAffyAnn.txt"

cortexDat<- read.table(cortexFile, header=T, row.names=1)
anno<- read.table(cortexAnn, header=T, row.names=1)

#3)****
#males<-cortexDat[,1:18]
#females<-cortexDat[,19:30]

#> malesM<-unname(as.matrix(males))
#> femalesM<-unname(as.matrix(females))

dat<-unname(as.matrix(cortexDat))

#cortexL50<-cortexDat
#cortexL50<-cortexL50[,-8:-18]
#cortexL50<-cortexL50[,-13:-19]

#> cortexO50<-cortexDat
#> cortexO50<-cortexO50[,-1:-7]
#> cortexO50<-cortexO50[,-12:-16]

#cortexL50<-unname(as.matrix(cortexL50))
#cortexO50<-unname(as.matrix(cortexO50))

#4.) Run the t.test function

> t.test.all.genes<-function(x,s1,s2){
+ x1<-x[s1]
+ x2<-x[s2]
+ x1<-as.numeric(x1)
+ x2<-as.numeric(x2)
+ t.out<-t.test(x1,x2,alternative="two.sided",var.equal=T)
+ out<-as.numeric(t.out$p.value)
+ return(out)
+ }

#Gender Comparison
> ganno<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
rawpGender<-apply(dat,1,t.test.all.genes,s1=ganno==0,s2=ganno==1)

#Age Comparison
> aanno<-c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1)
>rawpAge<-apply(dat,1,t.test.all.genes,s1=aanno==0,s2=aanno==1)

#Holm’s correction
>pAge.cor<-p.adjust(rawpAge,method="holm")
> pGender.cor<-p.adjust(rawpGender,method="holm")


#> t.outt<-t.test(malesM, femalesM, alternative="two.sided", var.equal=T)
#> rawpGender <-as.numeric(t.outt$p.value)
#>rawpGender

#> t.outt<-t.test(cortexL50, cortexO50, alternative="two.sided", var.equal=T)
#> rawpAge<-as.numeric(t.outt$p.value)
#> rawpAge

#5)****
> rawpAge<-sort(rawpAge)
> pAge.cor<-sort(pAge.cor)

> rawpGender<-sort(rawpGender)
> pGender.cor<-sort(pGender.cor)

> plot(rawpGender, type="l", main= "Adjusted p-values",ylab= "Sorted adjusted p-values")
> lines(pGender.cor, col="blue")
> legend("bottomright", legend= c("Raw P", "Adjusted"), fill=c("black", "blue"))

> plot(rawpAge, type="l", main= "Adjusted Age p-values",ylab= "Sorted adjusted p-values")
> lines(pAge.cor, col="blue")
> legend("bottomright", legend= c("Raw P", "Adjusted"), fill=c("black", "blue"))


#6.) Repeat #4 and #5 with the Bonferroni method.

#Bonferroni adjust P
>pAge.bon<-p.adjust(rawpAge,method="bonferroni")
> pGender.bon<-p.adjust(rawpGender,method="bonferroni")

> pAge.bon<-sort(pAge.bon)
> pGender.bon<-sort(pGender.bon)

> plot(rawpAge, type="l", main= "Adjusted Age p-values",ylab= "Sorted adjusted p-values")
> lines(pAge.bon, col="blue")
> legend("bottomright", legend= c("Raw P", "Bonferroni"), fill=c("black", "blue"))

> plot(rawpGender, type="l", main= "Adjusted Gender p-values",ylab= "Sorted adjusted p-values")
> lines(pGender.bon, col="blue")
> legend("bottomright", legend= c("Raw P", "Bonferroni"), fill=c("black", "blue"))

#7)****
>tcgaFile<- "C:/Users/fermi/Documents/Fall 2020/Lab 6/tcga_brca_fpkm.txt"
>tcgaAnFile<- "C:/Users/fermi/Documents/Fall 2020/Lab 6/tcga_brca_fpkm_sam.txt"

tcgaDat<- read.table(tcgaFile, header=T, row.names=1)
tcgaAnno<- read.table(tcgaAnFile, header=T, row.names=1)
