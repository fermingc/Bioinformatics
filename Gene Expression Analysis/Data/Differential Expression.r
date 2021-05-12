#1)****
library(Biobase)
library(annotate)
library(golubEsets)
library(multtest)
data(golub)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("golubEsets")

#2)****
#as.data.frame(x, row.names = NULL, optional = FALSE, ...)

seq<- seq(1:3051)
rown<-paste("g", seq, sep="")
dat<-as.data.frame(golub, row.names = rown)

#3.) Get the sample labels (see lecture notes) 
#and set the sample labels to the data frame. (2.5 pts)

colnames(dat)<-golub.cl

#4)****

wilcox.test.all.genes<-function(x,s1,s2){
x1<-x[s1]
x2<-x[s2]
x1<-as.numeric(x1)
x2<-as.numeric(x2)
w.out<-wilcox.test(x1,x2, exact=F, alternative="two.sided", correct=T)
out<-as.numeric(w.out$statistic)
return(out)
}

original.wmw.run<-apply(dat, 1, wilcox.test.all.genes, s1=golub.cl==0, s2=golub.cl==1)

#5)****

#dat_mixed<-dat[ , sample(1:ncol(dat))]

set.seed(123)
mts<-c()
i<-0

for (i in 1:500){
mix_dat<-sample(dat)
shuffled.wmw<-apply(mix_dat, 1, wilcox.test.all.genes, s1=golub.cl==0, s2=golub.cl==1)
mts<-c(mts,(max(shuffled.wmw)))
}

#mts:
# [1] 247 258 246 239 247 256 247 268 242 237 254 237 238 257 252 250 247 244 256 250 231 238 235 240 248 231 246 244 243 245 249 233 238 256 263 261 240 241
# [39] 268 248 243 258 254 235 248 254 243 237 251 241 246 261 253 264 243 253 242 263 266 261 265 246 240 246 234 251 244 234 254 239 240 245 255 272 239 241
# [77] 266 239 244 246 234 248 249 239 260 254 242 242 252 264 236 243 233 244 263 250 256 238 243 266 268 260 250 265 245 241 255 255 253 244 265 248 242 243
#[115] 248 242 248 264 235 245 227 245 262 233 251 245 245 247 238 240 245 251 240 261 244 239 258 240 233 242 243 258 254 274 240 249 252 236 247 248 245 252
#[153] 254 237 255 251 243 248 239 256 251 239 256 245 245 244 243 258 235 225 250 247 244 246 237 244 253 254 258 242 255 247 261 250 249 239 248 246 258 244
#[191] 240 249 234 244 265 237 238 238 250 256 249 248 245 250 264 256 239 242 236 242 257 245 238 254 254 260 251 239 244 265 241 253 253 247 239 252 257 257
#[229] 256 244 262 255 277 243 262 265 249 244 222 250 246 252 238 247 243 238 264 258 245 259 251 242 268 237 249 256 245 251 251 254 253 258 263 245 274 241
#[267] 247 239 251 257 257 251 256 246 249 263 251 238 229 246 249 241 237 254 259 258 246 236 261 243 239 253 252 241 251 240 250 250 237 239 251 243 237 260
#[305] 241 250 251 234 267 233 244 245 247 246 264 259 254 246 240 256 250 231 251 245 265 232 257 244 243 244 254 234 246 258 243 241 246 263 248 251 259 258
#[343] 256 234 246 266 235 252 257 247 244 253 227 249 264 243 254 253 248 270 247 253 242 260 232 243 235 239 233 242 253 246 252 245 253 245 246 237 262 254
#[381] 232 259 235 245 231 246 249 247 252 231 243 254 258 240 255 244 239 259 236 232 263 237 236 250 241 246 263 252 257 244 235 259 254 260 234 258 252 236
#[419] 242 248 236 256 247 255 241 261 252 242 234 261 250 248 251 241 240 237 241 260 243 257 252 232 253 257 248 240 235 244 243 247 237 243 240 237 244 255
#[457] 252 248 245 262 228 249 250 255 256 265 245 264 245 237 252 251 242 254 236 262 252 246 247 256 258 246 234 237 241 252 262 258 257 247 239 235 237 234
#[495] 260 253 232 236 243 253

#6)***
quantile(mts, c(.95))
#264

TSIndex<-which(original.wmv.run>264)
original.wmv.run[TSIndex]


#7)***

> zeros<-golub.cl[1:27]
> ones<-golub.cl[28:38]


design<-cbind(Grp1=1,Grp2vs1=c(rep(1,length(zeros)),rep(0,length(ones))))
fit<- lmFit(dat,design)
fit<-eBayes(fit)

#p-values in: fit$p.value[,2]

#8)***
sort(fit$p.value[,2])

> lowestBayes<-sort(fit$p.value[,2])
> lowestBayes[1:76]

#intersectin of: 
#lowestBayes[1:76] and original.wmv.run[TSIndex]

nBayes<-as.data.frame(lowestBayes[1:76])
> nwmv<- as.data.frame(original.wmv.run[TSIndex])

> intersect(row.names(nBayes),row.names(nwmv))


#9)****
t.test.all.genes<-function(x,s1,s2){
x1<-x[s1]
x2<-x[s2]
x1<-as.numeric(x1)
x2<-as.numeric(x2)
t.out<-t.test(x1,x2, alternative="two.sided", var.equal=TRUE)
out<-as.numeric(t.out$p.value)
return(out)
}

pv<-apply(dat, 1, t.test.all.genes, s1=golub.cl==0, s2=golub.cl==1)

PVIndex<-which(pv<.01)
#pv values < .01 in Student’s t-test = pv[PVIndex]
#same pv values from Bayes = fit$p.value[,2][PVIndex]

> plot(as.numeric(fit$p.value[,2][PVIndex]), as.numeric(pv[PVIndex]), xlab="Empirical Bayes", ylab="Student's T-Test", main="P-value Distribution Comparison Using Golub Data Set", cex=.5, col=3, pch=15)



