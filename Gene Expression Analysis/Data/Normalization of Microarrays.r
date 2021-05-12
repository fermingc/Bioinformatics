>GSM45<- "C:/Users/fermi/Documents/Fall 2020/HW2/GSM304445.gpr"
>GSM46<- "C:/Users/fermi/Documents/Fall 2020/HW2/GSM304446.gpr"
>GSM47<- "C:/Users/fermi/Documents/Fall 2020/HW2/GSM304447.gpr"
>GSM48<- "C:/Users/fermi/Documents/Fall 2020/HW2/GSM304448.gpr"


#1)****
> library(marray)

sample1<-"GSM304445.gpr"
sample2<-"GSM304446.gpr"
sample3<-"GSM304447.gpr"
sample4<-"GSM304448.gpr"

samples<-c(sample1,sample2,sample3,sample4)

dir.path <- "C:/Users/fermi/Documents/Fall 2020/HW2/"
a.cdna <- read.GenePix(path=dir.path,name.Gf = "F532 Median",name.Gb ="B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median",name.W ="Flags")

#a.cdna@maRf---> calls maRf data
#> str(a.cdna) --> gives the structure of a.cdna

#2.) Normalize each array using median global, loess, and print-tip-group loess methods.  
#Then plot MvA plots of all 4 arrays comparing no normalization to the other 3 normalization approaches. (2 pts)

#Median
> mnorm45<-maNorm(a.cdna[,1], norm="median", echo=TRUE) 
> mnorm46<-maNorm(a.cdna[,2], norm="median", echo=TRUE) 
> mnorm47<-maNorm(a.cdna[,3], norm="median", echo=TRUE) 
> mnorm48<-maNorm(a.cdna[,4], norm="median", echo=TRUE) 

#Loess
> lnorm45<-maNorm(a.cdna[,1], norm="loess", echo=TRUE) 
> lnorm46<-maNorm(a.cdna[,2], norm="loess", echo=TRUE) 
> lnorm47<-maNorm(a.cdna[,3], norm="loess", echo=TRUE) 
> lnorm48<-maNorm(a.cdna[,4], norm="loess", echo=TRUE) 

#PrintTipLoess
> PTnorm45<-maNorm(a.cdna[,1], norm="printTipLoess", echo=TRUE) 
> PTnorm46<-maNorm(a.cdna[,2], norm="printTipLoess", echo=TRUE) 
> PTnorm47<-maNorm(a.cdna[,3], norm="printTipLoess", echo=TRUE) 
> PTnorm48<-maNorm(a.cdna[,4], norm="printTipLoess", echo=TRUE) 

#nonorm, whole set
> a.cdna.no.norm <- maNorm(a.cdna,norm="none", echo=TRUE)

#MvA plots:
> par(mfrow=c(3,1),oma=c(0.1,10,0.1,10))

#Array1 (-45)
> maPlot(a.cdna.no.norm[,1],lines.func=NULL,legend.func=NULL,main='Not Normalized\n Sample #1')
> maPlot(mnorm45,lines.func=NULL,legend.func=NULL,main='Median Global Normalized\n Sample #1')
> maPlot(PTnorm45,lines.func=NULL,legend.func=NULL,main='Location Print-tip Loess Normalized\n Sample #1')
> maPlot(lnorm45,lines.func=NULL,legend.func=NULL,main='Loess Normalized\n Sample #1')

#Array2 (-46)
> maPlot(a.cdna.no.norm[,2],lines.func=NULL,legend.func=NULL,main='Not Normalized\n Sample #2')
> maPlot(mnorm46,lines.func=NULL,legend.func=NULL,main='Median Global Normalized\n Sample #2')
> maPlot(PTnorm46,lines.func=NULL,legend.func=NULL,main='Location Print-tip Loess Normalized\n Sample #2')
> maPlot(lnorm46,lines.func=NULL,legend.func=NULL,main='Loess Normalized\n Sample #2')

#Array3 (-47)
> maPlot(a.cdna.no.norm[,3],lines.func=NULL,legend.func=NULL,main='Not Normalized\n Sample #3')
> maPlot(mnorm47,lines.func=NULL,legend.func=NULL,main='Median Global Normalized\n Sample #3')
> maPlot(PTnorm47,lines.func=NULL,legend.func=NULL,main='Location Print-tip Loess Normalized\n Sample #3')
> maPlot(lnorm47,lines.func=NULL,legend.func=NULL,main='Loess Normalized\n Sample #3')

#Array4 (-48)
> maPlot(a.cdna.no.norm[,4],lines.func=NULL,legend.func=NULL,main='Not Normalized\n Sample #4')
> maPlot(mnorm48,lines.func=NULL,legend.func=NULL,main='Median Global Normalized\n Sample #4')
> maPlot(PTnorm48,lines.func=NULL,legend.func=NULL,main='Location Print-tip Loess Normalized\n Sample #4')
> maPlot(lnorm48,lines.func=NULL,legend.func=NULL,main='Loess Normalized\n Sample #4')

#3.) Plot density plots of the log ratio values for each normalization (and pre normalization) for only array #4.  
#Put them all on the same plot.  Make sure to label the axes and provide a legend. (2 pts)

maM(a.cdna[,4]) 
maM(mnorm48)
maM(PTnorm48)
maM(lnorm48)

>norm48dens<-as.numeric(maM(mnorm48))
>density(norm48dens, na.rm=TRUE)
> plot(density(norm48dens))

nonorm48dens<-na.omit((as.numeric(maM(a.cdna[,4]))))
> lines(density(nonorm48dens))

> PTnorm48dens<-na.omit((as.numeric(maM(PTnorm48))))
> lines(density(PTnorm48dens))

lnorm48dens<-na.omit((as.numeric(maM(lnorm48))))
lines(density(lnorm48dens))

> densnames=c("Median","No norm", "Print-tip","Loess")
> legend("topright",legend=densnames, fill=1:length(densnames))


#4.) Based on the plots generated so far, which normalization do you think is most preferred for this dataset? (2 pts)


#5)**** Cy5 = red
#subtract the background from the foreground values = maRf-maRb
a.cdna[,1]@maRf-(a.cdna[,1]@maRb)
a.cdna[,2]@maRf-(a.cdna[,2]@maRb)
a.cdna[,3]@maRf-(a.cdna[,3]@maRb)
a.cdna[,4]@maRf-(a.cdna[,4]@maRb)

BS45<-a.cdna[,1]@maRf-(a.cdna[,1]@maRb)
BS46<-a.cdna[,2]@maRf-(a.cdna[,2]@maRb)
BS47<-a.cdna[,3]@maRf-(a.cdna[,3]@maRb)
BS48<-a.cdna[,4]@maRf-(a.cdna[,4]@maRb)

#log2
#log2(BSXX)

BS45log2<-log2(BS45)
BS46log2<-log2(BS46)
BS47log2<-log2(BS47)
BS48log2<-log2(BS48)

#calculate global median normalization
> BS45log2M<-BS45log2/median(BS45log2,na.rm=TRUE)
> BS46log2M<-BS46log2/median(BS46log2,na.rm=TRUE)
> BS47log2M<-BS47log2/median(BS47log2,na.rm=TRUE)
> BS48log2M<-BS48log2/median(BS48log2,na.rm=TRUE)

#6)*** Spearman's rank Correlation
#cor(data, method = "spearman")

> AllBS<-matrix(nrow=44290, ncol=4)
> AllBS[,1]<-BS45[,1]
> AllBS[,2]<-BS46[,1]
> AllBS[,3]<-BS47[,1]
> AllBS[,4]<-BS48[,1]

> cor(AllBS, method="spearman")
> BSSpearmans<-cor(AllBS, method="spearman")

> chart.Correlation(AllBS, method="spearman")


#spearman M values from loess (all 4 in 1), scatter plot matrix

> AllLM<-matrix(nrow=44290, ncol=4)

> AllLM[,1]=maM(lnorm45)
> AllLM[,2]=maM(lnorm46)
> AllLM[,3]=maM(lnorm47)
> AllLM[,4]=maM(lnorm48)

> chart.Correlation(AllLM, method="spearman")


#7)*****

BS45<-a.cdna[,1]@maRf-(a.cdna[,1]@maRb)
BS46<-a.cdna[,2]@maRf-(a.cdna[,2]@maRb)
BS47<-a.cdna[,3]@maRf-(a.cdna[,3]@maRb)
BS48<-a.cdna[,4]@maRf-(a.cdna[,4]@maRb)

xanot<-matrix(nrow=44290, ncol=5)
x[,1]=array_names
x[,1]=BS45
x[,2]=BS46
x[,3]=BS47
x[,4]=BS48

x<-matrix(nrow=44290, ncol=4)
x[,1]=BS45
x[,2]=BS46
x[,3]=BS47
x[,4]=BS48

#7.2)
>xsort<-apply(x,2,sort)

#7.3)
>rowMeans(xsort)


#7.4)

xsortRmeans<-matrix(nrow=44290, ncol=4)
xsortRmeans[,1]=rowMeans(xsort)
xsortRmeans[,2]=rowMeans(xsort)
xsortRmeans[,3]=rowMeans(xsort)
xsortRmeans[,4]=rowMeans(xsort)

#7.5)
>x45rank<-rank(x[,1],ties.method="first")
>x46rank<-rank(x[,2],ties.method="first")
>x47rank<-rank(x[,3],ties.method="first")
>x48rank<-rank(x[,4],ties.method="first")

#7.6)

xnorm<-matrix(nrow=44290, ncol=4)

xnorm[,1]<-xsortRmeans[,1][x45rank]
xnorm[,2]<-xsortRmeans[,2][x46rank]
xnorm[,3]<-xsortRmeans[,3][x47rank]
xnorm[,4]<-xsortRmeans[,4][x48rank]

#8)***
#log2(BSXX)

xnormlog2<-log2(xnorm)

> chart.Correlation(xnormlog2, method="spearman")



