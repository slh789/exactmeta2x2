### Directory set-up and packages
rm(list=ls())
### select your directory here (make sure data file is in that directory)
dat <- read.csv("rosiglit.csv")
library(metafor)
library(meta)
library(forestplot)
library(poisbinom)
library("FEMetaBin") # install from github
library("numDeriv")
library(exact2x2)  ### for mid pvalue
library(polynom)
library(exactci)
library(poisbinom)
library(ggplot2)
library(grid)
library(gridExtra)
library(exactmeta)



#########################################
########## Avandia example ##############
#########################################


### Work with non-zero cases only -> select CVD death (.d) or MI (.mi)
# dat_trim <- dat[dat$rosi.d!=0|dat$cont.d!=0,]
# rosi_x <- dat_trim$rosi.d; rosi_n <- dat_trim$rosi.n
# cont_x <- dat_trim$cont.d; cont_n <- dat_trim$cont.n

dat_trim <- dat[dat$rosi.mi+dat$cont.mi!=0,]
rosi_x <- dat_trim$rosi.mi; rosi_n <- dat_trim$rosi.n
cont_x <- dat_trim$cont.mi; cont_n <- dat_trim$cont.n


### Initialized parameters
eta <- seq(0.05,0.95,l=19)
alpha <- 0.025
n <- nrow(dat_trim)
m <- length(eta)


### Find 1-sided confidence intervals for log(OR) using Fisher exact test with midp method
lims <- list()
for(j in 1:m){
  temp_lims <- matrix(rep(NA,2*n),n)
  for(i in 1:n){
    temp_lims[i,1] <- exact2x2(matrix(c(rosi_x[i],cont_x[i],
                               rosi_n[i]-rosi_x[i],cont_n[i]-cont_x[i]),2,2),
                               conf.level=eta[j],alternative="greater",midp=T)$conf.int[1]
    temp_lims[i,2] <- exact2x2(matrix(c(rosi_x[i],cont_x[i],
                               rosi_n[i]-rosi_x[i],cont_n[i]-cont_x[i]),2,2),
                               conf.level=eta[j],alternative="less",midp=T)$conf.int[2]
  }
  lims[[j]] <- log(temp_lims)
}


### Set weights and test set
w_i <- rosi_n+cont_n  ### weight is study size (page 276)
delta <- seq(-3,6,l=10000) ### select a set to test
w_j <- eta*(1-eta) ### page 277


### Write helper function to generate correlated Bernoulli's (page 277) -> if using multiple eta
bern_fun <- function(etas){
  mm <- length(etas)
  probs_eta <- c()
  for(i in 2:mm){probs_eta[i-1] <- etas[i]-etas[i-1]}
  probs_eta <- c(1-etas[mm],etas[1],probs_eta)
  vec_num <- sample(c(0:mm),size=1,replace=F,probs_eta)
  c(rep(0,mm-vec_num),rep(1,vec_num))
}


### Find critical point T(eta)
d_finder <- function(n,m,w_i,w_j,eta){
  mat_b <- t(as.matrix(replicate(n,bern_fun(eta))))
  Tnj <- matrix(rep(NA,m*2),m,2)
  for(j in 1:m){
    Tni <- matrix(rep(NA,n*2),n,2)
    for(i in 1:n){
      Tni[i,1] <- w_i[i]*(mat_b[i,j]-eta[j])*as.numeric(rosi_x[i]!=0)
      Tni[i,2] <- w_i[i]*(mat_b[i,j]-eta[j])*as.numeric(cont_x[i]!=0)
    }
    Tnj[j,1] <- w_j[j]*sum(Tni[,1])
    Tnj[j,2] <- w_j[j]*sum(Tni[,2])
  }
  lTn <- sum(Tnj[,1])
  uTn <- sum(Tnj[,2])
  c(lTn,uTn)
}
set.seed(10) ### choose random seed
d_vals <- t(replicate(10000,d_finder(n,m,w_i,w_j,eta)))


### Use t(eta) (equation 2.3)
tot_sums <- matrix(rep(NA,length(delta)*2),length(delta),2)
for(k in 1:length(delta)){
  tnj <- matrix(rep(NA,m*2),m,2)
  for(j in 1:m){
    tni <- matrix(rep(NA,n*2),n,2)
    for(i in 1:n){
      tni[i,1] <- w_i[i]*(as.numeric(delta[k]>lims[[j]][i,1])-eta[j])*as.numeric(rosi_x[i]!=0)
      tni[i,2] <- w_i[i]*(as.numeric(delta[k]<lims[[j]][i,2])-eta[j])*as.numeric(cont_x[i]!=0)
    }
    tnj[j,1] <- w_j[j]*sum(tni[,1])
    tnj[j,2] <- w_j[j]*sum(tni[,2])
  }
  tot_sums[k,1] <- sum(tnj[,1])
  tot_sums[k,2] <- sum(tnj[,2])
}


### Display results
pt_e <- sum(range(cbind(tot_sums,delta)[abs(tot_sums[,2]-tot_sums[,1])==min(abs(tot_sums[,2]-tot_sums[,1])),3]))/2
lb <- min(delta[tot_sums[,1]>=quantile(d_vals[,1],alpha)]); ub <- max(delta[tot_sums[,2]>=quantile(d_vals[,2],alpha)])
tianetal <- round(rbind(c(pt_e,lb,ub),exp(c(pt_e,lb,ub))),3)
rownames(tianetal) <- c("logOR","OR"); colnames(tianetal) <- c("est","lb","ub")
tianetal


#########################################
############ 2x2 Methods ################
#########################################

### Tian et al
###########################
# CVD D  est     lb    ub
# logOR 0.171 -0.185 0.548
# OR    1.187  0.831 1.730
###########################
# MI      est    lb    ub
# logOR 0.342 0.179 0.559
# OR    1.408 1.196 1.749


rm(list=ls())
### Read in and set-up data (again)
dat <- read.csv("rosiglit.csv")
# xk_mi <- dat$rosi.d
# mk_mi <- dat$rosi.n
# nk_mi <- dat$cont.n
# tk_mi <- dat$rosi.d+dat$cont.d
xk_mi <- dat$rosi.mi
mk_mi <- dat$rosi.n
nk_mi <- dat$cont.n
tk_mi <- dat$rosi.mi+dat$cont.mi


### Functions
multi2x2 <- function(mk, nk, tk){
  K <- length(mk)
  polycoeflist <- sapply(1:K, function(k){
    Re(polyroot( dhyper(0:(tk[k]), mk[k], nk[k], tk[k] ))) }, simplify=FALSE)
  polycoeflist
}

getp <- function(lpsi, lamall, xplus){
  dall <- dbinom( 0:length(lamall), length(lamall), mean(1/(1-lamall/exp(lpsi))) )
  sum(dall[dall<= dall[xplus+1]])	
}

lamstudy <- multi2x2(mk_mi,nk_mi,tk_mi)

### Use mathematica for MI endpoint to get lambda values for two of the larger studies, else skip this part
lamstudy[[42]] <- c(-3.35687, -3.25556, -3.1744, -3.10398, -3.04065, -2.98251, -2.92841, 
                    -2.87759, -2.8295, -2.78374, -2.74002, -2.69807, -2.65772, -2.61878, 
                    -2.58114, -2.54467, -2.50927, -2.47487, -2.44138, -2.40873, -2.37688, 
                    -2.34576, -2.31533, -2.28555, -2.25638, -2.22778, -2.19972, -2.17217, 
                    -2.14511, -2.11849, -2.09231, -2.06654, -2.04115, -2.01614, -1.99146, 
                    -1.96712, -1.94309, -1.91935, -1.89589, -1.87269, -1.84974, -1.82702, 
                    -1.80451, -1.7822, -1.76008, -1.73812, -1.71632, -1.69465, -1.6731, 
                    -1.65165, -1.63029, -1.60898, -1.58772, -1.56646, -1.5452, -1.52389, 
                    -1.5025, -1.48098, -1.4593, -1.43737, -1.41513, -1.39248, -1.36926, 
                    -1.3453, -1.32027, -1.29369, -1.26456, -1.23027)
lamstudy[[41]] <- c(-1.26469, -1.22776, -1.19768, -1.17118, -1.14699, -1.12446, 
                    -1.10319, -1.08291, -1.06342, -1.04459, -1.02629, -1.00842, 
                    -0.990897, -0.973644, -0.956585, -0.939643, -0.922737, -0.905774, 
                    -0.88864, -0.871185, -0.853192, -0.83431, -0.813869, -0.790103)



### Solving for the cMLE log OR
lamindiv <- lamstudy
lamstudy <- unlist(lamstudy)
psi0=1
negclik <- function(lpsi){
  if(sum(lamstudy==0)!=0){
    -dpoisbinom(sum(xk_mi)-1, 1/(1-lamstudy/exp(lpsi))[-which(lamstudy==0)], log_d=TRUE) } 
  else{
    -dpoisbinom(sum(xk_mi), 1/(1-lamstudy/exp(lpsi)), log_d=TRUE)  	  }
}
negclik <- Vectorize(negclik)
nlm1 <- nlm(negclik, p=0)
hess1 <- hessian( negclik, nlm1$est)

datd <- dat
# datd <- as.data.frame(cbind(datd$rosi.d,datd$rosi.n-datd$rosi.d,datd$cont.d,datd$cont.n-datd$cont.d))
datd <- as.data.frame(cbind(datd$rosi.mi,datd$rosi.n-datd$rosi.mi,datd$cont.mi,datd$cont.n-datd$cont.mi))
colnames(datd) <- c("ai","bi","ci","di")

for(i in 1:nrow(datd)){
  if(datd$ai[i]+datd$ci[i]==0){
    datd[i,c("ai","bi","ci","di")] <- datd[i,c("ai","bi","ci","di")]
  }
  else if(datd$ai[i]*datd$ci[i]==0){
    datd[i,c("ai","bi","ci","di")] <- datd[i,c("ai","bi","ci","di")] + 0.5
  }
  else {
    datd[i,c("ai","bi","ci","di")] <- datd[i,c("ai","bi","ci","di")]
  }
}

naive_ind <- rep(NA,nrow(datd))
for(i in 1:nrow(datd)){
  naive_ind[i] <- (datd$ai[i]*datd$di[i])/(datd$bi[i]*datd$ci[i])
}

phet <- c()
for(i in 1:length(naive_ind[!is.na(naive_ind)])){
  ptemp <- 1/(1-lamindiv[!is.na(naive_ind)][[i]]/naive_ind[!is.na(naive_ind)][i])
  phet <- c(phet,ptemp)
}

phom <- c()
for(i in 1:length(naive_ind[!is.na(naive_ind)])){
  ptemp <- 1/(1-lamindiv[!is.na(naive_ind)][[i]]/exp(nlm1$est))
  phom <- c(phom,ptemp)
}

deltaest <- sqrt(mean((phet-mean(phom))^2))  ### delta hat under heterogeneity
xiest <- sqrt(1/(1-deltaest^2/(mean(phom)*(1-mean(phom)))))
toolarge <- 1-2*pnorm(qnorm(.025)*xiest)-0.95

deltaesthom <- sqrt(mean((phom-mean(phom))^2))  ### delta hat under homogeneity
xiesthom <- sqrt(1/(1-deltaesthom^2/(mean(phom)*(1-mean(phom)))))
toolargehom <- 1-2*pnorm(qnorm(.025)*xiesthom)-0.95

rbind(c(deltaest,xiest,toolarge),c(deltaesthom,xiesthom,toolargehom))


### Forest plot of individual studies

#### For CVD mortality
# yyi <- log(naive_ind)
# vvi <- escalc(ai=xk_mi,bi=mk_mi-xk_mi,ci=tk_mi-xk_mi,di=nk_mi-tk_mi+xk_mi,measure="OR",drop00=T)$vi
# vcil <- yyi+qnorm(.025)/sqrt(vvi)
# vciu <- yyi+qnorm(.975)/sqrt(vvi)
# dtt <- as.data.frame(cbind(yyi,vcil,vciu))
# 
# imain.names <- paste0("Study ",c(1,5,8,11,13,15:16,18:20,25:28,30,32:35,39:42))
# iatvals <- c(23:1)
# iat.main <- c(23:1)
# 
# pdf("indiv_d.pdf",width=6,height=3.75,onefile=FALSE)
# par(mar=c(4,2,2,2)+0.1)
# with(dtt[!is.na(yyi),], plot(x=yyi[!is.na(yyi)], y=iatvals, pch=19, xlim=range(dtt[,2:3], na.rm=TRUE), 
#                              ylim=range(iatvals, na.rm=TRUE), axes=FALSE,
#                              xlab="log OR estimate", ylab=""))
# abline(v=0, lty=2, col="grey50")
# with(dtt[!is.na(yyi),], segments(x0=vcil, x1=vciu, y0=iatvals, y1=iatvals))
# axis(side=1, cex.axis=0.85)
# mtext(side=2, at = iat.main, imain.names, las=1, line=2, adj=0, cex=.5)
# dev.off()


#### For MI
yyi <- log(naive_ind)
vvi <- escalc(ai=xk_mi,bi=mk_mi-xk_mi,ci=tk_mi-xk_mi,di=nk_mi-tk_mi+xk_mi,measure="OR",drop00=T)$vi
vcil <- yyi+qnorm(.025)/sqrt(vvi)
vciu <- yyi+qnorm(.975)/sqrt(vvi)
dtt <- as.data.frame(cbind(yyi,vcil,vciu))

imain.names <- paste0("Study ",c(1:19,21:30,32,34:37,39:42))
iatvals <- c(38:1)
iat.main <- c(38:1)

pdf("indiv_mi.pdf",width=6,height=4.25,onefile=FALSE)
par(mar=c(4,2,2,2)+0.1)
with(dtt[!is.na(yyi),], plot(x=yyi[!is.na(yyi)], y=iatvals, pch=19, xlim=range(dtt[,2:3], na.rm=TRUE), 
              ylim=range(iatvals, na.rm=TRUE), axes=FALSE,
              xlab="log OR estimate", ylab=""))
abline(v=0, lty=2, col="grey50")
with(dtt[!is.na(yyi),], segments(x0=vcil, x1=vciu, y0=iatvals, y1=iatvals))
axis(side=1, cex.axis=0.85)
mtext(side=2, at = iat.main, imain.names, las=1, line=2, adj=0, cex=.5)
dev.off()


# MI   [,1]     [,2]        [,3]
# Het 0.1774393 1.070229 0.014060463
# Hom 0.1381158 1.040824 0.008647509

# D    [,1]     [,2]        [,3]
# Het 0.1746012 1.073477 0.014619531
# Hom 0.1065543 1.025568 0.005576888


### Getting the Blaker confidence intervals
### For CVD mortality
# bin.lb <- uniroot( function(lpsi){ getp( lpsi, lamstudy, sum(xk_mi))- 0.05}, c(-1,0))
# bin.ub <- uniroot( function(lpsi){ getp( lpsi, lamstudy, sum(xk_mi))- 0.05}, c(0,3))
### For MI
bin.lb <- uniroot( function(lpsi){ getp( lpsi, lamstudy, sum(xk_mi))- 0.05}, c(-1,.5))
bin.ub <- uniroot( function(lpsi){ getp( lpsi, lamstudy, sum(xk_mi))- 0.05}, c(.5,1))


### cMLE Blaker
cMLE_Blaker <- c( nlm1$est, bin.lb$root, bin.ub$root, getp(lpsi=0, lamall=lamstudy, xplus=sum(xk_mi)))

### cMLE Fisher Information
cMLE_Fish <- c(est=nlm1$est, ci.lb=nlm1$est-qnorm(0.975)/sqrt(hess1), 
              ci.ub=nlm1$est+qnorm(0.975)/sqrt(hess1), p=pchisq((nlm1$est-log(psi0))^2*hess1, df=1, lower=FALSE))

### Common-effect MH
MH=unlist( rma.mh(ai=xk_mi,bi=mk_mi-xk_mi,ci=tk_mi-xk_mi,di=nk_mi-tk_mi+xk_mi,measure="OR")[c("beta","ci.lb","ci.ub","pval")])

### Peto
Peto=unlist( rma.peto(ai=xk_mi,bi=mk_mi-xk_mi,ci=tk_mi-xk_mi,di=nk_mi-tk_mi+xk_mi)[c("beta","ci.lb","ci.ub","pval")])

### Common-effect Woolf
Woolf=as.numeric(unlist( rma.uni(ai=xk_mi,bi=mk_mi-xk_mi,ci=tk_mi-xk_mi,di=nk_mi-tk_mi+xk_mi,measure="OR",method="FE",drop00=TRUE) )[c("beta","ci.lb","ci.ub","pval")])

### Note that Li and Rice's method labels the cells differently with bi being the number of cases in the control arm
# datt <- dat[dat$rosi.d+dat$cont.d>0,] ### drop 00's
# datt <- as.data.frame(cbind(datt$rosi.d,datt$rosi.n-datt$rosi.d,datt$cont.d,datt$cont.n-datt$cont.d))

datt <- dat[dat$rosi.mi+dat$cont.mi>0,] ### drop 00's
datt <- as.data.frame(cbind(datt$rosi.mi,datt$rosi.n,datt$cont.mi,datt$cont.n))

colnames(datt) <- c("ai","n1i","bi","n2i")

### Li Hom MH
KendrickHomMH=unlist(femeta(ai=datt$ai,n1i=datt$n1i,bi=datt$bi,n2i=datt$n2i, null.value=log(psi0), vtype="Hom", drop00=TRUE)[c("estimate","conf.int","p.value")])

### Li Het MH
KendrickHetMH=unlist(femeta(ai=datt$ai,n1i=datt$n1i,bi=datt$bi,n2i=datt$n2i, null.value=log(psi0), drop00=TRUE)[c("estimate","conf.int","p.value")])

### Li Hom Woolf
KendrickHomWoolf=unlist(femeta(ai=datt$ai,n1i=datt$n1i,bi=datt$bi,n2i=datt$n2i, estimator="Woolf",null.value=log(psi0), vtype="Hom", drop00=TRUE)[c("estimate","conf.int","p.value")])

### Li Het Woolf
KendrickHetWoolf=unlist(femeta(ai=datt$ai,n1i=datt$n1i,bi=datt$bi,n2i=datt$n2i, estimator="Woolf",null.value=log(psi0), drop00=TRUE)[c("estimate","conf.int","p.value")])

### Li Hom MLE
KendrickHomMLE=unlist(femeta(ai=datt$ai,n1i=datt$n1i,bi=datt$bi,n2i=datt$n2i, estimator="CEMLE",null.value=log(psi0), vtype="Hom", drop00=TRUE)[c("estimate","conf.int","p.value")])

### Li Het MLE
KendrickHetMLE=unlist(femeta(ai=datt$ai,n1i=datt$n1i,bi=datt$bi,n2i=datt$n2i, estimator="CEMLE",null.value=log(psi0), drop00=TRUE)[c("estimate","conf.int","p.value")])

### Plain MLE
datr <- rbind(data.frame(aa = datt$ai, bb = datt$n1i-datt$ai, xval = 1, lab = 1:nrow(datt)), 
              data.frame(aa = datt$bi, bb = datt$n2i-datt$bi, xval = 0, lab = 1:nrow(datt)))
mod <- glm(cbind(aa, bb) ~ 0 + as.factor(lab) + xval, data = datr, family = "binomial")
PlainMLE <- c(summary(mod)$coef[nrow(datt)+1,][1],confint.default(mod)[nrow(datt)+1,],summary(mod)$coef[nrow(datt)+1,][4])


### For CVD death (input Tian et al numbers manually)
# tianetal <- c(0.171,-0.185,0.548,1) # p-value just a place holder
# dd <- as.data.frame(
#   round(rbind(cMLE_Fish,cMLE_Blaker,
#               tianetal,Peto,
#               MH,KendrickHomMH,KendrickHetMH,
#               Woolf,KendrickHomWoolf,KendrickHetWoolf,
#               PlainMLE,KendrickHomMLE,KendrickHetMLE)[,-4],3))
# tabletext_x <- cbind(c("0.498","0.498","0.171","0.495","0.529","0.183","0.183","0.269","0.206","0.206","0.512","0.511","0.511"),
#                      c("[-0.038,1.034]","[-0.072,1.082]","[-0.185,0.548]","[-0.020,1.010]","[-0.016,1.075]","[-0.589,0.506]","[-0.506,0.872]",
#                        "[-0.217,0.754]","[-0.417,0.828]","[-0.426,0.837]","[-0.023,1.047]","[-0.023,1.044]","[-0.004,1.025]"))
# main.names <- c("cMLE","Tian et al","Peto+0.5","Mantel-Haenszel","Woolf/IVW","Full MLE")
# sub.names <- c("plain","fixed effects w/Blaker",
#                "plain+0.5","common-effect L&R","fixed-effects L&R",
#                "plain+0.5","common-effect L&R","fixed-effects L&R",
#                "plain","common-effect L&R","fixed-effects L&R") 
# atvals <- c(18:17, 15, 13, 11:9, 7:5, 3:1)
# at.main <- c(18,15,13,11,7,3)
# at.sub  <- c(18:17, 11:9, 7:5, 3:1)
# 
# pdf("rosiglit_d.pdf",width=6,height=4,onefile=FALSE)
# par(mar=c(4,12,2,8.5)+0.1)
# with(dd, plot(x=est, y=atvals, pch=19, xlim=range(dd[,2:3], na.rm=TRUE), 
#               ylim=range(atvals, na.rm=TRUE), axes=FALSE,
#               xlab="log OR estimate", ylab=""))
# abline(v=0, lty=2, col="grey50")
# with(dd, segments(x0=ci.lb, x1=ci.ub, y0=atvals, y1=atvals))
# axis(side=1, cex.axis=0.85)
# mtext(side=2, at = at.main, main.names, las=1, line=12, adj=0)
# mtext(side=2, at = at.sub ,  sub.names, las=1, line=0.5, adj=1)
# mtext(side=4, at = atvals, tabletext_x[,1], las=1, line=1, adj=0.5)
# mtext(side=4, at = atvals, tabletext_x[,2], las=1, line=5, adj=0.5)
# mtext(side=2, at=19.5, "CVD mortality data: Method", las=1, adj=0, line=12, font=2)
# mtext(side=4, at=19.5, "estimate", line=1, las=1, font=2, adj=0.5)
# mtext(side=4, at=19.5, "95% CI", line=5, las=1, font=2, adj=0.5)
# dev.off()

### For MI (input Tian et al numbers manually)
tianetal <- c(0.342,0.179,0.559,0) # p-value just a place holder
dd <- as.data.frame(
        round(rbind(cMLE_Fish,cMLE_Blaker,
        tianetal,Peto,
        MH,KendrickHomMH,KendrickHetMH,
        Woolf,KendrickHomWoolf,KendrickHetWoolf,
        PlainMLE,KendrickHomMLE,KendrickHetMLE)[,-4],3))
tabletext_x <- cbind(c("0.355","0.355","0.342","0.356","0.356","0.346","0.346","0.251","0.259","0.259","0.355","0.355","0.355"),
                      c("[0.010,0.694]","[0.029,0.681]","[0.179,0.559]","[0.030,0.683]","[0.029,0.682]","[-0.174,0.866]","[-0.051,0.743]",
                        "[-0.062,0.565]","[-0.103,0.621]","[-0.106,0.623]","[0.029,0.681]","[0.029,0.681]","[0.030,0.680]"))
main.names <- c("cMLE","Tian et al","Peto+0.5","Mantel-Haenszel","Woolf/IVW","Full MLE")
sub.names <- c("plain","fixed effects w/Blaker",
               "plain+0.5","common-effect L&R","fixed-effects L&R",
               "plain+0.5","common-effect L&R","fixed-effects L&R",
               "plain","common-effect L&R","fixed-effects L&R") 
atvals <- c(18:17, 15, 13, 11:9, 7:5, 3:1)
at.main <- c(18,15,13,11,7,3)
at.sub  <- c(18:17, 11:9, 7:5, 3:1)

pdf("rosiglit_mi.pdf",width=6,height=4,onefile=FALSE)
par(mar=c(4,12,2,8.5)+0.1)
with(dd, plot(x=est, y=atvals, pch=19, xlim=range(dd[,2:3], na.rm=TRUE), 
              ylim=range(atvals, na.rm=TRUE), axes=FALSE,
              xlab="log OR estimate", ylab=""))
abline(v=0, lty=2, col="grey50")
with(dd, segments(x0=ci.lb, x1=ci.ub, y0=atvals, y1=atvals))
axis(side=1, cex.axis=0.85)
mtext(side=2, at = at.main, main.names, las=1, line=12, adj=0)
mtext(side=2, at = at.sub ,  sub.names, las=1, line=0.5, adj=1)
mtext(side=4, at = atvals, tabletext_x[,1], las=1, line=1, adj=0.5)
mtext(side=4, at = atvals, tabletext_x[,2], las=1, line=5, adj=0.5)
mtext(side=2, at=19.5, "MI data: Method", las=1, adj=0, line=12, font=2)
mtext(side=4, at=19.5, "estimate", line=1, las=1, font=2, adj=0.5)
mtext(side=4, at=19.5, "95% CI", line=5, las=1, font=2, adj=0.5)
dev.off()





#########################################
############ Enumeration ################
#########################################

rm(list=ls())



### Write functions
multi2x2 <- function(mk, nk, tk){
  K <- length(mk)
  polycoeflist <- sapply(1:K, function(k){
    Re(polyroot( dhyper(max(0,tk[k]-nk[k]):min(tk[k],mk[k]), mk[k], nk[k], tk[k] ))) }, simplify=FALSE)
  polycoeflist
}

crity <- function(p, n, lims){
  c( max( which(lims[1,]<p)) - 1, min(which(lims[2,]>p)) - 1 )
}

pjk_fun <- function(lams,lpsif,lpsi1,lpsi2){
  sum(1/(1-lams[[1]]/exp(lpsi1)))+sum(1/(1-lams[[2]]/exp(lpsi2)))-sum(1/(1-unlist(lams)/exp(lpsif)))
}

get_pvecs <- function(lams,lpsi1,lpsi2){
  p1 <- 1/(1-lams[[1]]/exp(lpsi1))
  p2 <- 1/(1-lams[[2]]/exp(lpsi2))              
  pf <- mean(c(p1,p2))
  delta2 <- mean( (c(p1,p2)-pf)^2 )
  list(p1=p1, p2=p2, pf=pf, delta2=delta2)
}

get_psif_delta <- function(lams,lpsi1,lpsi2){
  lpsif <- uniroot( function(lpsif){ pjk_fun(lams, lpsif, lpsi1, lpsi2)}, c(min(lpsi1,lpsi2),max(lpsi1,lpsi2)))$root
  u1 <- length(lams[[1]])
  u2 <- length(lams[[2]])
  pf <- mean(1/(1-unlist(lams)/exp(lpsif)))
  delta2 <- sum( 1/(u1+u2) * ( 1/(1-lams[[1]]/exp(lpsi1)) - pf )^2 ) +
    sum( 1/(u1+u2) * ( 1/(1-lams[[2]]/exp(lpsi2)) - pf )^2 )
  list(lpsif=lpsif, delta2=delta2, pf=pf)
}

getcov.exc <- function(lams,lims.bl,lpsi1,lpsi2){
  lpsif <- uniroot( function(lpsif){ pjk_fun(lams, lpsif, lpsi1, lpsi2)}, c(min(lpsi1,lpsi2),max(lpsi1,lpsi2)))$root
  pp <- get_pvecs(lams,lpsi1,lpsi2)
  delta2 <- get_psif_delta(lams,lpsi1,lpsi2)$delta2
  uk <- length(unlist(lams))
  cc <- crity(pp$pf, uk, lims.bl)
  covhom <- if(pp$pf==0 | pp$pf==1) return(1) else if( cc[2] == 0 ) pbinom(cc[1], uk, pp$pf) else pbinom(cc[1], uk, pp$pf) - pbinom(cc[2]-1, uk, pp$pf)
  covhet <- if( cc[2] == 0 ) ppoisbinom(cc[1], c(pp$p1,pp$p2)) else if( cc[2] > 1) ppoisbinom(cc[1], c(pp$p1,pp$p2)) - ppoisbinom(cc[2]-1, c(pp$p1,pp$p2)) else ppoisbinom(cc[1], c(pp$p1,pp$p2)) - dpoisbinom(0, c(pp$p1,pp$p2))
  c(lpsif=lpsif,homcov=covhom,hetcov=covhet,exccov=covhet-covhom,delta2=delta2)
}



### Run functions
mk_fix <- c(500,500); nk_fix <- c(500,500); tk_fix <- c(25,25)
lams.ex <- multi2x2(mk=mk_fix,nk=nk_fix,tk=tk_fix)
lims.bl <- sapply(0:length(unlist(lams.ex)),function(x){binom.exact(x,length(unlist(lams.ex)),tsmethod="blaker",control=exactci:::binomControl(tol=1E-7))$conf.int})
size.g <- 501
lpsi.g <- expand.grid(lpsi1=seq(-1,1,l=size.g), lpsi2=seq(-1,1,l=size.g)+0.001)
lpsi.g <- subset(lpsi.g, lpsi1 != lpsi2 )
lpsi.cov <- as.data.frame(t(apply( lpsi.g, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])})))
lpsi.cov <- as.data.frame(cbind(lpsi.g,lpsi.cov))
lpsi.p <- lpsi.cov[lpsi.cov$lpsi1<lpsi.cov$lpsi2,]
lpsi.p <- lpsi.p[order(lpsi.p$lpsif),]
line15 <- contourLines(x=seq(-1,1,l=size.g),
                       y=seq(-1,1,l=size.g)+0.001,
                       z=matrix(sqrt(lpsi.cov$delta2),size.g,size.g), levels=c(0.15))
line20 <- contourLines(x=seq(-1,1,l=size.g),
                       y=seq(-1,1,l=size.g)+0.001,
                       z=matrix(sqrt(lpsi.cov$delta2),size.g,size.g), levels=c(0.20))
lpsi15 <- data.frame(lpsi1=line15[[1]]$x, lpsi2=line15[[1]]$y)
lpsi15 <- cbind(lpsi15, t( apply(lpsi15, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])}) ))
lpsi15 <- lpsi15[lpsi15$lpsi1<lpsi15$lpsi2,]
lpsi20 <- data.frame(lpsi1=line20[[1]]$x, lpsi2=line20[[1]]$y)
lpsi20 <- cbind(lpsi20, t( apply(lpsi20, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])}) ))
lpsi20 <- lpsi20[lpsi20$lpsi1<lpsi20$lpsi2,]
p15 <- as.vector( apply( lpsi15, 1, function(x){get_psif_delta(lams.ex, x[1], x[2])$pf} ) )
p15 <- sapply(p15,function(x){1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-0.15^2)))-0.95})
p20 <- as.vector( apply( lpsi20, 1, function(x){get_psif_delta(lams.ex, x[1], x[2])$pf} ) )
p20 <- sapply(p20,function(x){1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-0.20^2)))-0.95})

pdf("abs_cov_bal.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = hetcov)) + geom_contour_filled() + theme_bw() + ggtitle("Absolute coverage under heterogeneity") + 
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("exc_cov_bal.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = exccov)) + geom_contour_filled() + theme_bw() + ggtitle("Excess coverage under heterogeneity") +
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("delta_bal.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = sqrt(delta2))) + geom_contour_filled() + theme_bw() + ggtitle(expression(delta)) + 
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("homog_cov_bal.pdf",width=6,height=4)
par(mar=c(4,4,1,0.1))
plot(lpsi.p$lpsif,lpsi.p$homcov,type="l",ylim=c(0.95,0.975),
     xlab=expression(paste("log(", psi[F],")")),
     ylab="Coverage under homogeneity")
legend("topleft",bty='n',legend=expression(paste("M"[k],"=N"[k],"=(500,500), T"[k],"=(25,25)")))
dev.off()

pdf("exc_lines_bal.pdf",width=6,height=3)
par(mar=c(4,4,1,0.1))
plot(lpsi15$lpsif,lpsi15$exccov,type="l",lty=1,xlim=c(-1,1),ylim=c(0,.02),
     xlab=expression(paste("log(", psi[F],")")),
     ylab="Excess coverage")
lines(lpsi15$lpsif,p15,lty=3,col="blue")
lines(lpsi20$lpsif,lpsi20$exccov,lty=2)
lines(lpsi20$lpsif,p20,lty=3,col="blue")
text(.57,.0095,expression(paste(delta,"=.15")))
text(.29,.017,expression(paste(delta,"=.20")))
dev.off()


### Run functions
mk_fix <- c(500,500); nk_fix <- c(500,500); tk_fix <- c(15,35)
lams.ex <- multi2x2(mk=mk_fix,nk=nk_fix,tk=tk_fix)
lims.bl <- sapply(0:length(unlist(lams.ex)),function(x){binom.exact(x,length(unlist(lams.ex)),tsmethod="blaker",control=exactci:::binomControl(tol=1E-7))$conf.int})
size.g <- 501
lpsi.g <- expand.grid(lpsi1=seq(-1,1,l=size.g), lpsi2=seq(-1,1,l=size.g)+0.001)
lpsi.g <- subset(lpsi.g, lpsi1 != lpsi2 )
lpsi.cov <- as.data.frame(t(apply( lpsi.g, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])})))
lpsi.cov <- as.data.frame(cbind(lpsi.g,lpsi.cov))
lpsi.p <- lpsi.cov[lpsi.cov$lpsi1<lpsi.cov$lpsi2,]
lpsi.p <- lpsi.p[order(lpsi.p$lpsif),]
line15 <- contourLines(x=seq(-1,1,l=size.g),
                       y=seq(-1,1,l=size.g)+0.001,
                       z=matrix(sqrt(lpsi.cov$delta2),size.g,size.g), levels=c(0.15))
line20 <- contourLines(x=seq(-1,1,l=size.g),
                       y=seq(-1,1,l=size.g)+0.001,
                       z=matrix(sqrt(lpsi.cov$delta2),size.g,size.g), levels=c(0.20))
lpsi15 <- data.frame(lpsi1=line15[[1]]$x, lpsi2=line15[[1]]$y)
lpsi15 <- cbind(lpsi15, t( apply(lpsi15, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])}) ))
lpsi15 <- lpsi15[lpsi15$lpsi1<lpsi15$lpsi2,]
lpsi20 <- data.frame(lpsi1=line20[[1]]$x, lpsi2=line20[[1]]$y)
lpsi20 <- cbind(lpsi20, t( apply(lpsi20, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])}) ))
lpsi20 <- lpsi20[lpsi20$lpsi1<lpsi20$lpsi2,]
p15 <- as.vector( apply( lpsi15, 1, function(x){get_psif_delta(lams.ex, x[1], x[2])$pf} ) )
p15 <- sapply(p15,function(x){1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-0.15^2)))-0.95})
p20 <- as.vector( apply( lpsi20, 1, function(x){get_psif_delta(lams.ex, x[1], x[2])$pf} ) )
p20 <- sapply(p20,function(x){1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-0.20^2)))-0.95})

pdf("abs_cov_unbal.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = hetcov)) + geom_contour_filled() + theme_bw() + ggtitle("Absolute coverage under heterogeneity") + 
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("exc_cov_unbal.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = exccov)) + geom_contour_filled() + theme_bw() + ggtitle("Excess coverage under heterogeneity") +
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("delta_unbal.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = sqrt(delta2))) + geom_contour_filled() + theme_bw() + ggtitle(expression(delta)) + 
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("homog_cov_unbal.pdf",width=6,height=4)
par(mar=c(4,4,1,0.1))
plot(lpsi.p$lpsif,lpsi.p$homcov,type="l",ylim=c(0.95,0.975),
     xlab=expression(paste("log(", psi[F],")")),
     ylab="Coverage under homogeneity")
legend("topleft",bty='n',legend=expression(paste("M"[k],"=N"[k],"=(500,500), T"[k],"=(15,35)")))
dev.off()

pdf("exc_lines_unbal.pdf",width=6,height=3)
par(mar=c(4,4,1,0.1))
plot(lpsi15$lpsif,lpsi15$exccov,type="l",lty=1,xlim=c(-1,1),ylim=c(0,.02),
     xlab=expression(paste("log(", psi[F],")")),
     ylab="Excess coverage")
lines(lpsi15$lpsif,p15,lty=3,col="blue")
lines(lpsi20$lpsif,lpsi20$exccov,lty=2)
lines(lpsi20$lpsif,p20,lty=3,col="blue")
text(.75,.01,expression(paste(delta,"=.15")))
text(.65,.0175,expression(paste(delta,"=.20")))
dev.off()


### Run functions
mk_fix <- c(400,600); nk_fix <- c(350,650); tk_fix <- c(25,25)
lams.ex <- multi2x2(mk=mk_fix,nk=nk_fix,tk=tk_fix)
lims.bl <- sapply(0:length(unlist(lams.ex)),function(x){binom.exact(x,length(unlist(lams.ex)),tsmethod="blaker",control=exactci:::binomControl(tol=1E-7))$conf.int})
size.g <- 501
lpsi.g <- expand.grid(lpsi1=seq(-1,1,l=size.g), lpsi2=seq(-1,1,l=size.g)+0.001)
lpsi.g <- subset(lpsi.g, lpsi1 != lpsi2 )
lpsi.cov <- as.data.frame(t(apply( lpsi.g, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])})))
lpsi.cov <- as.data.frame(cbind(lpsi.g,lpsi.cov))
lpsi.p <- lpsi.cov[lpsi.cov$lpsi1<lpsi.cov$lpsi2,]
lpsi.p <- lpsi.p[order(lpsi.p$lpsif),]
line15 <- contourLines(x=seq(-1,1,l=size.g),
                       y=seq(-1,1,l=size.g)+0.001,
                       z=matrix(sqrt(lpsi.cov$delta2),size.g,size.g), levels=c(0.15))
line20 <- contourLines(x=seq(-1,1,l=size.g),
                       y=seq(-1,1,l=size.g)+0.001,
                       z=matrix(sqrt(lpsi.cov$delta2),size.g,size.g), levels=c(0.20))
lpsi15 <- data.frame(lpsi1=line15[[1]]$x, lpsi2=line15[[1]]$y)
lpsi15 <- cbind(lpsi15, t( apply(lpsi15, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])}) ))
lpsi15 <- lpsi15[lpsi15$lpsi1<lpsi15$lpsi2,]
lpsi20 <- data.frame(lpsi1=line20[[1]]$x, lpsi2=line20[[1]]$y)
lpsi20 <- cbind(lpsi20, t( apply(lpsi20, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2])}) ))
lpsi20 <- lpsi20[lpsi20$lpsi1<lpsi20$lpsi2,]
p15 <- as.vector( apply( lpsi15, 1, function(x){get_psif_delta(lams.ex, x[1], x[2])$pf} ) )
p15 <- sapply(p15,function(x){1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-0.15^2)))-0.95})
p20 <- as.vector( apply( lpsi20, 1, function(x){get_psif_delta(lams.ex, x[1], x[2])$pf} ) )
p20 <- sapply(p20,function(x){1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-0.20^2)))-0.95})

pdf("abs_cov_mn.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = hetcov)) + geom_contour_filled() + theme_bw() + ggtitle("Absolute coverage under heterogeneity") + 
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("exc_cov_mn.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = exccov)) + geom_contour_filled() + theme_bw() + ggtitle("Excess coverage under heterogeneity") +
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("delta_mn.pdf",width=6,height=4)
ggplot(lpsi.cov, aes(lpsi1, lpsi2, z = sqrt(delta2))) + geom_contour_filled() + theme_bw() + ggtitle(expression(delta)) + 
  xlab(expression(paste("log(", psi[1],")"))) + ylab(expression(paste("log(", psi[2],")")))
dev.off()

pdf("homog_cov_mn.pdf",width=6,height=4)
par(mar=c(4,4,1,0.1))
plot(lpsi.p$lpsif,lpsi.p$homcov,type="l",ylim=c(0.95,0.975),
     xlab=expression(paste("log(", psi[F],")")),
     ylab="Coverage under homogeneity")
legend("topleft",bty='n',legend=expression(paste("M"[k],"=(400,600), N"[k],"=(350,650), T"[k],"=(25,25)")))
dev.off()

pdf("exc_lines_mn.pdf",width=6,height=3)
par(mar=c(4,4,1,0.1))
plot(lpsi15$lpsif,lpsi15$exccov,type="l",lty=1,xlim=c(-1,1),ylim=c(0,.02),
     xlab=expression(paste("log(", psi[F],")")),
     ylab="Excess coverage")
lines(lpsi15$lpsif,p15,lty=3,col="blue")
lines(lpsi20$lpsif,lpsi20$exccov,lty=2)
lines(lpsi20$lpsif,p20,lty=3,col="blue")
text(.57,.01,expression(paste(delta,"=.15")))
text(.29,.0175,expression(paste(delta,"=.20")))
dev.off()




### Write functions for five strata
multi2x2 <- function(mk, nk, tk){
  K <- length(mk)
  polycoeflist <- sapply(1:K, function(k){
    Re(polyroot( dhyper(max(0,tk[k]-nk[k]):min(tk[k],mk[k]), mk[k], nk[k], tk[k] ))) }, simplify=FALSE)
  polycoeflist
}

crity <- function(p, n, lims){
  c( max( which(lims[1,]<p)) - 1, min(which(lims[2,]>p)) - 1 )
}

pjk_fun <- function(lams,lpsif,lpsi1,lpsi2,lpsi3,lpsi4,lpsi5){
  sum(1/(1-lams[[1]]/exp(lpsi1)))+sum(1/(1-lams[[2]]/exp(lpsi2)))+
    sum(1/(1-lams[[3]]/exp(lpsi3)))+sum(1/(1-lams[[4]]/exp(lpsi4)))+
    sum(1/(1-lams[[5]]/exp(lpsi5)))-sum(1/(1-unlist(lams)/exp(lpsif)))
}

get_pvecs <- function(lams,lpsi1,lpsi2,lpsi3,lpsi4,lpsi5){
  p1 <- 1/(1-lams[[1]]/exp(lpsi1))
  p2 <- 1/(1-lams[[2]]/exp(lpsi2))
  p3 <- 1/(1-lams[[3]]/exp(lpsi3))
  p4 <- 1/(1-lams[[4]]/exp(lpsi4))
  p5 <- 1/(1-lams[[5]]/exp(lpsi5))
  pf <- mean(c(p1,p2,p3,p4,p5))
  delta2 <- mean( (c(p1,p2,p3,p4,p5)-pf)^2 )
  list(p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, pf=pf, delta2=delta2)
}

getcov.exc <- function(lams,lims.bl,lpsi1,lpsi2,lpsi3,lpsi4,lpsi5){
  lpsif <- uniroot( function(lpsif){ pjk_fun(lams, lpsif, lpsi1, lpsi2, lpsi3, lpsi4, lpsi5)}, c(min(lpsi1, lpsi2, lpsi3, lpsi4, lpsi5),max(lpsi1, lpsi2, lpsi3, lpsi4, lpsi5)))$root
  pp <- get_pvecs(lams, lpsi1, lpsi2, lpsi3, lpsi4, lpsi5)
  uk <- length(unlist(lams))
  cc <- crity(pp$pf, uk, lims.bl)
  covhom <- if(pp$pf==0 | pp$pf==1) return(1) else if( cc[2] == 0 ) pbinom(cc[1], uk, pp$pf) else pbinom(cc[1], uk, pp$pf) - pbinom(cc[2]-1, uk, pp$pf)
  covhet <- if( cc[2] == 0 ) ppoisbinom(cc[1], c(pp$p1,pp$p2,pp$p3,pp$p4,pp$p5)) else if( cc[2] > 1) ppoisbinom(cc[1], c(pp$p1,pp$p2,pp$p3,pp$p4,pp$p5)) - ppoisbinom(cc[2]-1, c(pp$p1,pp$p2,pp$p3,pp$p4,pp$p5)) else ppoisbinom(cc[1], c(pp$p1,pp$p2,pp$p3,pp$p4,pp$p5)) - dpoisbinom(0, c(pp$p1,pp$p2,pp$p3,pp$p4,pp$p5))
  c(lpsif=lpsif,homcov=covhom,hetcov=covhet,exccov=covhet-covhom,delta2=pp$delta2)
}

### Run functions for five strata
mk_fix <- c(200,200,200,200,200); nk_fix <- c(200,200,200,200,200); tk_fix <- c(5,10,20,10,5)
lams.ex <- multi2x2(mk=mk_fix,nk=nk_fix,tk=tk_fix)
lims.bl <- sapply(0:length(unlist(lams.ex)),function(x){binom.exact(x,length(unlist(lams.ex)),tsmethod="blaker",control=exactci:::binomControl(tol=1E-7))$conf.int})
expander_funct <- function(center,spread){
  c(center-spread,center-spread/2,center,center+spread/2,center+spread)
}
lpsi.0 <- as.data.frame(t(sapply(seq(.01,1,l=501),expander_funct,center=0)))
lpsi.cov0 <- as.data.frame(t(apply( lpsi.0, 1, function(x){getcov.exc(lams.ex, lims.bl, x[1], x[2], x[3], x[4], x[5])})))

p0 <- as.vector( apply( lpsi.0, 1, function(x){get_pvecs(lams.ex, x[1], x[2], x[3], x[4], x[5])$pf} ) )
p0del <- sqrt( as.vector( apply( lpsi.0, 1, function(x){get_pvecs(lams.ex, x[1], x[2], x[3], x[4], x[5])$delta2} ) ) )
p0.0 <- c()
for(i in 1:length(p0)){
  x <- p0[i]
  p0.0[i] <- 1-2*pnorm(qnorm(0.025) * sqrt(x*(1-x))/sqrt((x*(1-x)-lpsi.cov0$delta2[i])))-0.95
}

pdf("5_delta_excess.pdf",width=6,height=3)
par(mar=c(4,4,1,0.1))
plot(sqrt(lpsi.cov0$delta2),lpsi.cov0$exccov,type="l",ylim=c(0,.015),
     xlab=expression(delta),
     ylab="Excess coverage")
lines(sqrt(lpsi.cov0$delta2),p0.0,lty=3,col="blue")
dev.off()




