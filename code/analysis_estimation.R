###########################################################################
##                                                                       ##
## New estimators for conditional false discovery rates improve power    ##
##  and smoothness                                                       ##
##                                                                       ##
## Analyse results of simulations, draw plots, and analyse TWAS data     ##
##                                                                       ##
## James Liley and Chris Wallace                                         ##
##                                                                       ##
###########################################################################
#
# This R script will deterministically generate all plots and outputs used
#  in the paper above. Outputs are generated in the same order that they 
#  appear in the paper, with supplementary material generated according to
#  the first referenced in the main  text. 
#
# This script requires all requisite packages to be installed and the 
#  working directory to be set to the file one level up from where this 
#  script is saved (it should contain subdirectories ./code, ./data, 
#  ./simulations and ./plots). All outputs are written to ./outputs.
#
# Simulations are generated using the R script ./code/run_simulation.R. 
#  Each run of the simulation is dependent on a random seed. A matrix of 
#  completed simulation outputs with associated random seeds is stored in 
#  ./data/cfdrsimmatrix.RData.
#
# Input data from the TWAS analysis is included in the R object 
#  ./data/twas_data.RData. This file contains XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#
# The analysis of TWAS data runs quite slowly. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX It will only run if the 
#  file ./data/gwas_cfdr_data.RData does not already exist. A copy of
#  this file is already included in the ./data directory. To re-run the 
#  GWAS analysis, delete this file and run this script. The file 
#  generated will be identical to that included in the directory.
#  

###########################################################################
## Packages, directories and switches #####################################
###########################################################################

# Packages

library(cfdr)
library(mnormt)
library(mgcv)
library(pbivnorm)
library(MASS)
library(fields)
library(matrixStats)
library(plotrix)


# Directories

# This variable gives the location of the folder cfdr_estimation. It should 
#  contain a document README, and four folders  named ./code, ./data, 
#  ./simulations, and ./outputs. Here it is assumed to be the current working 
#  directory
cfdr_dir="./"

# This variable gives the location to save outputs (plots, tables) to.
output_dir=paste0(cfdr_dir,"outputs/")
data_dir=paste0(cfdr_dir,"data/")


# Switches

# set to F to draw plots rather than saving them as PDFs
save_pdf=T 




###########################################################################
## Read data ##############################################################
###########################################################################

# Load file from directory ./data/. See file ./data/cfdrsimmatrix_README 
#  for details.
load(paste0(data_dir,"cfdrsimmatrix.RData"))



###########################################################################
## Figure showing chaoticity of cFDR1 #####################################
###########################################################################

# This figure does not require any data

if (save_pdf) pdf(paste0(output_dir,"figure_1.pdf"),width=5,height=5)

x1=0.001; y1=0.0025;
x2=0.0025; y2=0.001;
m=max(c(x1,y1,x2,y2))*1.5
m0=min(c(x1,y1,x2,y2))*0.75
m2=m0/2

plot(0,type="n",xlim=c(0,m),ylim=c(0,m),bty="n",xaxs="i",yaxs="i",xlab="P",ylab="Q",xaxt="n",yaxt="n")
axis(1,at=c(0,x1,x2,1),labels=c("0",expression(paste(p[1], " (min P)")),expression(p[2]),""))
axis(2,at=c(0,x1,x2,1),labels=c("0",expression(paste(q[2]," (min Q)")),expression(q[1]),""))


s=seq(0,m0,length.out=50); rs=50:1
z=outer(s,s,function(x,y) sqrt(x^2 + y^2))
cx=function(gl) c(gray(1 - gl*(100:1)/101),rep(gray(1),100))
image(x1+s,y1+s,z,col=cx(0.4),add=T)
image(rev(x1-s),rev(y1-s),z[rs,rs],col=cx(0.9),add=T)
image(rev(x1-s),y1+s,z[rs,],col=cx(1),add=T)
image(x1+s,rev(y1-s),z[,rs],col=cx(0.9),add=T)

image(x2+s,y2+s,z,col=cx(0.2),add=T)
image(rev(x2-s),rev(y2-s),z[rs,rs],col=cx(0.2),add=T)
image(rev(x2-s),y2+s,z[rs,],col=cx(0.6),add=T)
image(x2+s,rev(y2-s),z[,rs],col=cx(0.2),add=T)



points(c(x1,x2),c(y1,y2),pch=16)

np=200; Z=cbind(rnorm(np,mean=max(x1,x2),sd=m0*0.4),rnorm(np,mean=max(y1,y2),sd=m0*0.4)); 
#Z=Z[which(Z[,1]> min(x1,x2)+m2 & Z[,1]< max(x1,x2)+m2 & Z[,2]>min(y1,y2)+ m2 & Z[,2]<max(y1,y2)- m2),]
ccol=c(gray(0.8),gray(0.2))[1+ (Z[,2] < max(y1,y2))]
points(Z,pch=16,col=ccol,cex=0.5)

arrows(m*0.9,y1,m*0.9,y2,code=3); text(m*0.95,(y1+y2)/2,expression(N[Q]))
#arrows(x1,m*0.9,x2,m*0.9,code=3); text((x1+x2)/2,m*0.95,expression(N[P]))

segments(c(x1,x2,x1,x2),c(y1,y2,y1,y2),c(x1,x2,0,0),c(0,0,y1,y2),lty=3)

points(x1+ c(-1,-1,1,1)*m2/2,y1+c(-1,1,-1,1)*m2/2,col="red",pch=16)
points(x2+ c(-1,-1,1,1)*m2/2,y2+c(-1,1,-1,1)*m2/2,col="red",pch=16)

sc=1.2
text(x1+sc*m2/2,y1+sc*m2/2,expression(paste("p'(N"[Q],"+1)/2")),col="red",adj=c(0,0))
text(x1-sc*m2/2,y1+sc*m2/2,expression(paste("p'(N"[Q],"+1)")),col="red",adj=c(1,0))
text(x1+sc*m2/2,y1-sc*m2/2,expression(paste("p'N"[Q])),col="red",adj=c(0,1))
text(x1-sc*m2/2,y1-sc*m2/2,expression(paste("p'N"[Q])),col="red",adj=c(1,1))


text(x2+sc*m2/2,y2+sc*m2/2,"p'",col="red",adj=c(0,0))
text(x2-sc*m2/2,y2+sc*m2/2,"2p'",col="red",adj=c(1,0))
text(x2+sc*m2/2,y2-sc*m2/2,"p'",col="red",adj=c(0,1))
text(x2-sc*m2/2,y2-sc*m2/2,"p'",col="red",adj=c(1,1))


if (save_pdf) dev.off()



###########################################################################
## Different types of L-curves with different cFDR types ##################
###########################################################################

# Simulation
set.seed(1) # reproducibility

N=10000; n1p=100; n1q=100; n1pq=100; sp=2; sq=2

zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(N-n1p-n1pq,sd=1))
zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(N-n1p-n1pq-n1q,sd=1))

# P-values
p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

# True parameters for mixture Gaussian model
pars=c((N-n1p-n1q-n1pq)/N,n1p/N,n1q/N,sp,sp,sq,sq)

# Subset of curves at which to plot curves/points
xmin=5e-4
sub=which(p< xmin)


# Draw L curves
if (save_pdf) pdf(paste0(output_dir,"Lplot_setcut_F.pdf"),width=9,height=3)

cx=seq(0.1,0.8,length.out=8)
v1=vl(p,q,at=cx,adj=F);
v2=vlx(p,q,pars=pars,at=cx,adj=F);
v3=vly(p,q,at=cx,adj=F); 

ex1=expression(atop({L^1}[S](alpha),alpha %in% ("0.1,0.2,...,0.8")))
ex2=expression(atop({L^2}[S](alpha),alpha %in% ("0.1,0.2,...,0.8")))
ex3=expression(atop({L^3}[S](alpha),alpha %in% ("0.1,0.2,...,0.8")))

par(mfrow=c(1,3))
plot(0,0,type="n",xlab="P",ylab="Q",xlim=c(0,0.01),ylim=c(0,1),main=ex1)
for (i in 1:length(cx)) lines(v1$x[i,],v1$y,col="gray")
points(p,q,col="red",cex=0.5)

plot(0,0,type="n",xlab="P",ylab="Q",xlim=c(0,0.01),ylim=c(0,1),main=ex2)
for (i in 1:length(cx)) lines(v2$x[i,],v2$y,col="gray")
points(p,q,col="red",cex=0.5)

plot(0,0,type="n",xlab="P",ylab="Q",xlim=c(0,0.01),ylim=c(0,1),main=ex3)
for (i in 1:length(cx)) lines(v3$x[i,],v3$y,col="gray")
points(p,q,col="red",cex=0.5)

if (save_pdf) dev.off()








###########################################################################
## Comparison of power by cFDR types on method 3b #########################
###########################################################################

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
if ("rx" %in% search()) detach(rx)
attach(rx)

for (plot_type in 0:3) {

xdist=plot_type
if (plot_type==0) xdist=1:3

normalise=T # set to F to draw absolute power rather than relative to t2r_p
ltx=1

if (save_pdf) pdf(paste0(output_dir,"power_fdp3_",paste(xdist,collapse=""),".pdf"),width=4,height=4)

if (xdist[1]==1) ptitle="Normal"
if (xdist[1]==2) ptitle="t (3df)"
if (xdist[1]==3) ptitle="Cauchy"
if (length(xdist)==3) ptitle="All"


# shorthands
fp=t2r_p
f1=t2r_cf1_fdr3b_adj0_dist1 
f2=t2r_cf2_fdr3b_adj0_dist1; 
f3=t2r_cf3_fdr3b_adj0_dist1; 
xx=(n1p+n1pq)

# subset of points to plot at. Use simulations from continuous distributions to better show trends with n1p+n1pq
w=which(pmin(fp,f1,f2,f3)> -0.5 & # remove simulations which failed to finish in time
          xx>0 & # t2r is indeterminate if xx==0
          dist1 %in% xdist &
          !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions


fp=fp[w]; f1=f1[w]; f2=f2[w]; f3=f3[w]; xx=xx[w]

sdx=0.15 # kernel for gaussian smoothing
sdx=sdx*(max(xx)-min(xx)) # scale kernel

ux=sort(unique(xx))
y1=rep(0,length(ux)); y2=rep(0,length(ux)); y3=rep(0,length(ux)); yp=rep(0,length(ux))

for (i in 1:length(ux)) {
  wts=dnorm(xx-ux[i],sd=sdx); wts=wts/sum(wts)
  y1[i]=sum(f1*wts);     y2[i]=sum(f2*wts);     y3[i]=sum(f3*wts); yp[i]=sum(fp*wts)
}

if (normalise) {
  y1=y1-yp; y2=y2-yp; y3=y3-yp; yp=yp-yp
}

if (normalise) ylab=expression(paste(Delta,"T2R")) else ylab="T2R"

plot(0,type="n",cex=0.5,col="gray",main=ptitle,
     ylim=range(c(y1,y2,y3,yp)), xlim=range(xx), 
     xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=ylab)
lines(sort(ux),y1[order(ux)],col="black",lwd=2,lty=ltx);
lines(sort(ux),y2[order(ux)],col="darkgray",lwd=2,lty=ltx);
lines(sort(ux),y3[order(ux)],col="lightgray",lwd=2,lty=ltx);
lines(sort(ux),yp[order(ux)],col="black",lwd=2,lty=2);

abline(h=0.1,col="black",lty=3)

lp1=range(c(y1,y2,y3))
legend( "bottomright", #0.45*max(xx), lp1[1]+ 0.4*(lp1[2]-lp1[1]),
        c(expression(widehat(cFDR)^1),expression(widehat(cFDR)^2),expression(widehat(cFDR)^3)),
        lty=c(1,1,1),col=c("black","darkgray","lightgray"),bg="white")


if (save_pdf) dev.off()

}

detach(rx)



###########################################################################
## Comparison of power between adjusted and unadjusted cFDR ###############
###########################################################################


suppressWarnings(rm(list=intersect(colnames(rx),ls())))
if ("rx" %in% search()) detach(rx)
attach(rx)

normalise=T
ltx=1 # type of lines for unadjusted
ltx2=2 # type of lines for adjusted

if (save_pdf) pdf(paste0(output_dir,"power_fdp3_adj_123.pdf"),width=4,height=4)

# shorthands
fp=t2r_p
f1=t2r_cf1_fdr3b_adj0_dist1 
f2=t2r_cf2_fdr3b_adj0_dist1; 
f3=t2r_cf3_fdr3b_adj0_dist1; 
f1s=t2r_cf1_fdr3b_adj1_dist1 
f2s=t2r_cf2_fdr3b_adj1_dist1; 
f3s=t2r_cf3_fdr3b_adj1_dist1; 
xx=(n1p+n1pq)

# subset of points to plot at. Use simulations from continuous distributions to better show trends with n1p+n1pq
w=which(pmin(fp,f1,f2,f3,f1s,f2s,f3s)> -0.5 & # remove simulations which failed to finish in time
          xx>0 & # t2r is indeterminate if xx==0
          !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions


fp=fp[w]; f1=f1[w]; f2=f2[w]; f3=f3[w]; 
f1s=f1s[w]; f2s=f2s[w]; f3s=f3s[w];
xx=xx[w]; 

sdx=0.15 # kernel for gaussian smoothing
sdx=sdx*(max(xx)-min(xx)) # scale kernel

ux=sort(unique(xx))
y1=rep(0,length(ux)); y2=rep(0,length(ux)); y3=rep(0,length(ux)); yp=rep(0,length(ux))
y1s=y1; y2s=y2; y3s=y3

for (i in 1:length(ux)) {
  wts=dnorm(xx-ux[i],sd=sdx); wts=wts/sum(wts)
  y1[i]=sum(f1*wts);     y2[i]=sum(f2*wts);     y3[i]=sum(f3*wts); 
  y1s[i]=sum(f1s*wts);  y2s[i]=sum(f2s*wts); y3s[i]=sum(f3s*wts);
  yp[i]=sum(fp*wts)
}

if (normalise) {
  y1=y1-yp; y2=y2-yp; y3=y3-yp; 
  y1s=y1s-yp; y2s=y2s-yp; y3s=y3s-yp
  yp=yp-yp
}

if (normalise) ylab=expression(paste(Delta,"T2R")) else ylab="T2R"

plot(0,type="n",cex=0.5,col="gray",
     ylim=range(c(y1,y2,y3,y1s,y2s,y3s,yp)), xlim=range(xx), 
     xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=ylab)
lines(sort(ux),y1[order(ux)],col="black",lwd=2,lty=ltx);
lines(sort(ux),y2[order(ux)],col="darkgray",lwd=2,lty=ltx);
lines(sort(ux),y3[order(ux)],col="lightgray",lwd=2,lty=ltx);
lines(sort(ux),y1s[order(ux)],col="black",lwd=2,lty=ltx2);
lines(sort(ux),y2s[order(ux)],col="darkgray",lwd=2,lty=ltx2);
lines(sort(ux),y3s[order(ux)],col="lightgray",lwd=2,lty=ltx2);
lines(sort(ux),yp[order(ux)],col="black",lwd=2,lty=2);

abline(h=0.1,col="black",lty=3)

lp1=range(c(y1,y2,y3))
legend( "bottomright", #0.45*max(xx), lp1[1]+ 0.4*(lp1[2]-lp1[1]),
        c(expression(widehat(cFDR)^1),expression(widehat(cFDR)^2),expression(widehat(cFDR)^3)),
        lty=c(1,1,1),col=c("black","darkgray","lightgray"),bg="white")


if (save_pdf) dev.off()

detach(rx)




###########################################################################
## Local vs standard cFDR #################################################
###########################################################################


suppressWarnings(rm(list=intersect(colnames(rx),ls())))
if ("rx" %in% search()) detach(rx)
attach(rx)


normalise=T
ltx=1

# Shorthands
fp=t2r_p
f1=t2r_cf2_fdr3b_adj0_dist1; 
f2=t2r_cf2_fdr4_adj0_dist1; 
xx=(n1p+n1pq)


# subset of points to plot at. Use simulations from continuous distributions to better show trends with n1p+n1pq
w=which(pmin(f1,f2,fp) > -0.5 & # remove simulations which failed to finish in time
          xx>0 & # t2r is indeterminate if xx==0
          !(N %in% c(1000,10000))) # use simulations where variables are sampled from 'continuous' distributions

fp=fp[w]; f1=f1[w]; f2=f2[w]; xx=xx[w]

sdx=0.15 # gaussian smoothing
ux=sort(unique(xx))
y1=rep(0,length(ux));  y2=rep(0,length(ux)); yp=rep(0,length(ux))

sdx=sdx*(max(xx)-min(xx))
for (i in 1:length(ux)) {
  wts=dnorm(xx-ux[i],sd=sdx); wts=wts/sum(wts)
  y1[i]=sum(f1*wts);     y2[i]=sum(f2*wts);    yp[i]=sum(fp*wts)
}

if (normalise) {
    y1=y1-yp; y2=y2-yp; yp=yp-yp
}


if (save_pdf) pdf(paste0(output_dir,"power_local.pdf"),width=4,height=4)

if (normalise) ylab=expression(paste(Delta,"T2R")) else ylab="T2R"
plot(0,type="n",cex=0.5,col="gray",ylim=range(c(y1,y2,yp)), xlim=range(xx), xlab=expression(paste(n[1]^p, "+ n"[1]^{pq})),ylab=ylab)
lines(sort(ux),y1[order(ux)],col="black",lwd=2,lty=ltx);
lines(sort(ux),y2[order(ux)],col="darkgray",lwd=2,lty=ltx);
lines(sort(ux),yp[order(ux)],col="black",lwd=2,lty=2);

abline(h=0.1,col="black",lty=3)

lp1=range(c(y1,y2))
legend( "bottomright", #0.45*max(xx), lp1[1]+ 0.4*(lp1[2]-lp1[1]),
        c(expression(widehat(cFDR)^2),expression(widehat(cfdr))),
        lty=c(1,1,1),col=c("black","darkgray"),bg="white")


if (save_pdf) dev.off()

detach(rx)





###########################################################################
## Numerical T2R comparison between methods ###############################
###########################################################################

# The output of this section of code is saved in ./outputs as a text file
#  ./outputs/T2R_comparison.txt.

suppressWarnings(rm(list=intersect(colnames(rx),ls())))
if ("rx" %in% search()) detach(rx)
attach(rx)



# Shorthands. 
tp=t2r_p; # T2R for p-values (BH procedure, for reference)
t1=t2r_cf1_fdr3b_adj0_dist1; # T2R using cFDR1, no adjustment
t1s=t2r_cf1_fdr3b_adj1_dist1; # T2R using cFDR1, adjusted
t2=t2r_cf2_fdr3b_adj0_dist1; # T2R using cFDR2, no adjustment
t2s=t2r_cf2_fdr3b_adj1_dist1; # T2R using cFDR2, adjusted
t3=t2r_cf3_fdr3b_adj0_dist1; # T2R using cFDR3, no adjustment
t3s=t2r_cf3_fdr3b_adj1_dist1; # T2R using cFDR3, adjusted
tx=t2r_cf2_fdr4_adj0_dist1; # T2R using local cfdr
xx=n1p+n1pq

# We only consider values with n1p+n1pq>0, since T2R is indeterminate if 
#  n1p+n1pq==0
x1=which(xx>0)


# T2R comparison between methods, using paired Wilcoxon rank sum tests.
#  The sign of the difference is given by the sign of the pseudomedian;
#  that is, if the pseudomedian of wilcox.test(x,y,...) is >0, then 
#  in general x>y.

# Comparisons between cFDR type and p-value
wilcox.test(t1[x1],tp[x1],paired=T,conf.int=T) # cFDR1 is more powerful than p-value
wilcox.test(t2[x1],tp[x1],paired=T,conf.int=T) # cFDR2 is more powerful than p-value
wilcox.test(t3[x1],tp[x1],paired=T,conf.int=T) # cFDR3 is more powerful than p-value
wilcox.test(t1s[x1],tp[x1],paired=T,conf.int=T) # cFDR1s is more powerful than p-value
wilcox.test(t2s[x1],tp[x1],paired=T,conf.int=T) # cFDR2s is more powerful than p-value
wilcox.test(t3s[x1],tp[x1],paired=T,conf.int=T) # cFDR3s is more powerful than p-value
wilcox.test(tx[x1],tp[x1],paired=T,conf.int=T) # cfdr is NOT clearly more powerful than p-value

# Comparisons between cFDR1, cFDR2, cFDR3
wilcox.test(t1[x1],t2[x1],paired=T,conf.int=T) # cFDR1 is more powerful than cFDR2
wilcox.test(t1[x1],t3[x1],paired=T,conf.int=T) # cFDR1 is more powerful than cFDR3
wilcox.test(t2[x1],t3[x1],paired=T,conf.int=T) # cFDR2 is more powerful than cFDR3



# Comparisons between cFDR1s, cFDR2s, cFDR3s
wilcox.test(t1s[x1],t2s[x1],paired=T,conf.int=T) # cFDR1s is more powerful than cFDR2s
wilcox.test(t1s[x1],t3s[x1],paired=T,conf.int=T) # cFDR1s is more powerful than cFDR3s
wilcox.test(t2s[x1],t3s[x1],paired=T,conf.int=T) # cFDR2s is more powerful than cFDR3s


# Comparisons between cFDR1s, cFDR1, cFDR2s, cFDR2, cFDR3s, cFDR3.
wilcox.test(t1s[x1],t1[x1],paired=T,conf.int=T) # cFDR1s is more powerful than cFDR1
wilcox.test(t2s[x1],t2[x1],paired=T,conf.int=T) # cFDR2s is more powerful than cFDR2
wilcox.test(t3s[x1],t3[x1],paired=T,conf.int=T) # cFDR3s is more powerful than cFDR3


# It is difficult to sensibly calculate the power of a Wilcoxon test; however,
#  the test is powerful. Approximately 25% of T2Rs are equal in each comparison. We show below
#  that we have approximately 90% power to detect a 2% difference in 
#  P(FDP(method A) > FDP(method B)) and P(FDP(method A) < FDP(method B))
delta=0.04 # P(FDP(method A) > FDP(method B)) - P(FDP(method A) < FDP(method B))
equal=0.25 # P(FDP(method A)=FDP(method B))
n=length(x1); np=0; ntrial=1000; P=0.05; set.seed(1)
for (i in 1:ntrial) {
 s=sample(c(-1,0,1),n,prob=c((1-equal-delta)/2,equal,(1-equal+delta)/2),rep=T)
 pw=wilcox.test(s,rep(0,n),paired=T)$p.value; if (pw< P) np=np+1
}
np/ntrial # power to detect difference

detach(rx)




###########################################################################
## Demonstrate chaoticity of cfdr #########################################
###########################################################################


#### PLOT FOR LOCAL cfdr ####

if (save_pdf) pdf(paste0(output_dir,"local_vs_cdf1.pdf"),width=4,height=4)

N=10000; n1p=200; n1q=200; n1pq=200; sp=3; sq=3 # parameters for simulation
pars=c((N-n1p-n1q-n1pq)/N,n1p/N,n1q/N,sp,sp,sq,sq) # true parameters of mixture-Gaussian

p0=2*pnorm(-4); q0=1 # force curve through point (p0,q0) but don't include it in points used to fit parameters

v1=vlxl(p0,q0,pars=pars,indices=1,scale="z") # compute L-regions

plot(0,xlim=c(0,8),ylim=c(0,8),type="n",xaxs="i",yaxs="i",
     xlab=expression("Z"[P]),ylab=expression("Z"[Q]),
     main=expression(widehat("cfdr")^2))
lines(v1$x[1,],v1$y)

for (i in 1:6) {
set.seed(i)

zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(N-n1p-n1pq,sd=1))
zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(N-n1p-n1pq-n1q,sd=1))

p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

pars2=fit.4g(cbind(p,q))$pars

v2=vlxl(c(p0,p),c(q0,q),pars=pars2,indices=1,scale="z")
lines(v2$x[1,],v2$y,col="red")
}

rlim=2
wr=which(zp^2 + zq^2 > rlim^2)
points(abs(zp[wr]),abs(zq[wr]),cex=0.5)
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black") # don't plot all points where dense, to cut down on plot size

legend("topright",c("True","Est."),col=c("black","red"),lty=1)

if (save_pdf) dev.off()



#### PLOT FOR STANDARD cFDR ####

if (save_pdf) pdf(paste0(output_dir,"local_vs_cdf2.pdf"),width=4,height=4)



# N=10000; n1p=200; n1q=200; n1pq=200; sp=3; sq=3 # as for previous plot
# pars=c((N-n1p-n1q-n1pq)/N,n1p/N,n1q/N,sp,sp,sq,sq)

p0=2*pnorm(-4); q0=1 # force curve through point (p0,q0) but don't include it in points used to fit parameters

v1=vlx(p0,q0,pars=pars,indices=1,scale="z",adj=0) # compute L regions


plot(0,xlim=c(0,8),ylim=c(0,8),type="n",xaxs="i",yaxs="i",
     xlab=expression("Z"[P]),ylab=expression("Z"[Q]),
     main=expression(widehat("cFDR")^2))
lines(v1$x[1,],v1$y)

for (i in 1:6) {
set.seed(i)

zp=c(rnorm(n1pq,sd=sp),rnorm(n1p,sd=sp),rnorm(N-n1p-n1pq,sd=1))
zq=c(rnorm(n1pq,sd=sq),rnorm(n1p,sd=1),rnorm(n1q,sd=sq),rnorm(N-n1p-n1pq-n1q,sd=1))

p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

pars2=fit.4g(cbind(p,q))$pars

v2=vlx(c(p0,p),c(q0,q),pars=pars2,indices=1,scale="z",adj=0)
lines(v2$x[1,],v2$y,col="red")
}

rlim=2
wr=which(zp^2 + zq^2 > rlim^2)
points(abs(zp[wr]),abs(zq[wr]),cex=0.5)
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")

legend("topright",c("True","Est."),col=c("black","red"),lty=1)

dev.off()



###########################################################################
## TWAS analysis ##########################################################
###########################################################################
#
# The computation of L-regions for the cFDR data is fairly slow. This code
#  will generate the output file ./data/twas_summary.RData only if the file
#  does not already exist.


###########################################################################
## Generate summary statistics, L regions, and v-values ###################
###########################################################################

outfile=paste0(data_dir,"twas_summary.RData")

if (!file.exists(outfile)) {


## Read raw data - BRCA
rxbrca0=read.table(paste0(data_dir,"BCAC.dat"),header=T)
l1=levels(rxbrca0$PANEL); l1x=l1[which(grepl("GTEx",l1))]; rxbrca=rxbrca0[which(rxbrca0$PANEL %in% l1x),]
cname1=paste0(rxbrca$PANEL,rxbrca$ID); rxbrca=rxbrca[match(unique(cname1),cname1),]; cname1=unique(cname1); rownames(rxbrca)=cname1

## Read raw data - OCA
rxoca0=read.table(paste0(home_dir,"OCAC.dat"),header=T)
l1=levels(rxoca0$PANEL); l1x=l1[which(grepl("GTEx",l1))]; rxoca=rxoca0[which(rxoca0$PANEL %in% l1x),]
cname1=paste0(rxoca$PANEL,rxoca$ID); rxoca=rxoca[match(unique(cname1),cname1),]; cname1=unique(cname1); rownames(rxoca)=cname1

# Extract p-values and assign each gene-tissue pair a name.
p_brca=rxbrca$TWAS.P; names(p_brca)=paste0(rxbrca$PANEL,rxbrca$ID)
p_oca=rxoca$TWAS.P; names(p_oca)=paste0(rxoca$PANEL,rxoca$ID)



# Merge BRCA and OCA datasets
inx=intersect(names(p_brca),names(p_oca))
p_brca=p_brca[inx]; p_oca=p_oca[inx]

# remove NAs
w=which(is.finite(-qnorm(p_brca/2)+ -qnorm(p_oca/2)))
p_brca=p_brca[w]; p_oca=p_oca[w]

# Z-scores
z_brca=-qnorm(p_brca/2); z_oca=-qnorm(p_oca/2)


# Folds for CV: same gene is always in the same fold
fold=as.numeric(as.factor(id))
nfold=max(fold)

# Other data
chr=as.character(rxbrca[names(p_brca),]$chr)
pos=as.numeric(rxbrca[names(p_brca),]$pos)
gene=as.character(rxbrca[names(p_brca),]$ID)
tissue=as.character(rxbrca[names(p_brca),]$PANEL)



### BRCA

sub_brca=which(z_oca^2 + z_brca^2 > 4^2 ) # only bother with potentially interesting variables

xbrca1=matrix(-1,length(sub_brca),2004);  # L-regions for cFDR1s
xbrca2=xbrca1 # L-regions for cFDR1

r_brca=rank(p_brca)

totalx=0
for (xfold in 1:max(fold)) {
inx=which(fold==xfold)
sub=intersect(sub_brca,inx)
if (length(sub)>0) {
  v1=vl(p_brca,p_ov,adj=T,indices=sub,fold=inx,mode=2,nv=2000)
  v2=vl(p_brca,p_ov,adj=F,indices=sub,fold=inx,mode=2,nv=2000)
  xbrca1[match(sub,sub_brca),]=v1$x; ybrca1=v1$y
  xbrca2[match(sub,sub_brca),]=v2$x; ybrca2=v2$y
} else v1=NULL
totalx=totalx+length(sub)
print(c(xfold,totalx))
}


#### OCA 

sub_oca=which(z_oca^2 + z_brca^2 > 4^2 )

xoca1=matrix(-1,length(sub_oca),2004);  # L-regions for cFDR1s
xoca2=xoca1 # L-regions for cFDR1

r_oca=rank(p_oca)

totalx=0
for (xfold in 1:max(fold)) {
inx=which(fold==xfold)
sub=intersect(sub_oca,inx)
if (length(sub)>0) {
  v1=vl(p_oca,p_brca,adj=T,indices=sub,fold=inx,mode=2,nv=2000)
  v2=vl(p_oca,p_brca,adj=F,indices=sub,fold=inx,mode=2,nv=2000)
  xoca1[match(sub,sub_oca),]=v1$x; yoca1=v1$y
  xoca2[match(sub,sub_oca),]=v2$x; yoca2=v2$y
} else v1=NULL
totalx=totalx+length(sub)
print(c(xfold,totalx))
}

### Integrate over L-regions 

# Parametrisation of Q|H0: BRCA|OCA
pars_brca=fit.2g(p_oca[which(p_brca> 0.5)])$pars
pi0_brca=pars_brca[1]
sigma_brca=pars_brca[2]

# Parametrisation of Q|H0: OCA|BRCA
pars_oca=fit.2g(p_brca[which(p_oca> 0.5)])$pars
pi0_oca=pars_oca[1]
sigma_oca=pars_oca[2]

# Integrate over L: BRCA|OCA
i_brca1=il(xbrca1,ybrca1,pi0_null=pi0_brca,sigma_null=sigma_brca,dist="norm") # adjusted cFDR
i_brca2=il(xbrca2,ybrca2,pi0_null=pi0_brca,sigma_null=sigma_brca,dist="norm") # unadjusted cFDR

# Integrate over L: OCA|BRCA
i_oca1=il(xoca1,yoca1,pi0_null=pi0_oca,sigma_null=sigma_oca,dist="norm") # adjusted cFDR
i_oca2=il(xoca2,yoca2,pi0_null=pi0_oca,sigma_null=sigma_oca,dist="norm") # unadjusted cFDR


### Save everything 

save(p_brca,p_oca,fold,sub_brca,sub_oca,
     chr,pos,gene,tissue,
     xbrca1,xbrca2,xoca1,xoca2,ybrca1,ybrca2,yoca1,yoca2,
     pi0_brca,sigma_brca,pi0_oca,sigma_oca,
     i_brca1,i_brca2,i_oca1,i_oca2,
     file=outfile)

}




###########################################################################
## Analyse cFDR results ###################################################
###########################################################################

## Load cFDR data
load(outfile)

for (brca in c(FALSE,TRUE)) { # run BRCA|OCA, then OCA|BRCA

# shorthands
if (brca) {
 xv1=xbrca1; xv2=xbrca2; yv1=ybrca1; yv2=ybrca2; sub=sub_brca
 p=p_brca; q=p_oca; v1=i_brca1; v2=i_brca2; pi0=pi0_brca; sigma=sigma_brca; 
} else {
 xv1=xoca1; xv2=xoca2; yv1=yoca1; yv2=yoca2; sub=sub_oca
 p=p_oca; q=p_brca; v1=i_oca1; v2=i_oca2; pi0=pi0_oca; sigma=sigma_oca; 
}

k=length(sub); # number of values for which L-regions were computed
xsub=3:2002 # set of indices of finite-valued points of L-curves

# FDR control level
alpha=1e-6

# Benjamini-Hochberg method on p-values
r=rank(p)
hp=which(r <= max(r[which(p/r < alpha/length(p))]))

# Benjamini-Hochberg method on v-values, adjusted
vbc1=rep(1,length(p)); vbc1[sub]=v1; rvbc1=rank(vbc1)
mvx1=max(rvbc1[which(vbc1/rvbc1 < alpha/length(vbc1))]); mvxl1=which(sub==match(mvx1,rvbc1))
h1=which(rvbc1 <= mvx1)

# Benjamini-Hochberg method on v-values, unadjusted
vbc2=rep(1,length(p)); vbc2[sub]=v2; rvbc2=rank(vbc2)
mvx2=max(rvbc2[which(vbc2/rvbc2 < alpha/length(vbc2))]); mvxl2=which(sub==match(mvx2,rvbc2))
h2=which(rvbc2 <= mvx2)


zp=-qnorm(p/2); zq=-qnorm(q/2)


## Draw plot of Z-scores and rejection regions

if (brca) header="BRCA|OCA" else header="OCA|BRCA"
if (brca) plot_name=paste0("twas_regions_brca.pdf") else plot_name=paste0("twas_regions_oca.pdf")

# save as pdf if switch is set
if (save_pdf) pdf(paste0(output_dir,plot_name),width=5,height=5)


rlim=3; osub=which(zp^2 + zq^2 > rlim^2)
plot(zp[osub],zq[osub],cex=0.3,xlim=c(0,8),ylim=c(0,8),xaxs="i",yaxs="i",
     main=header,xlab=c("Z (BRCA)","Z (OCA)")[2-brca],ylab=c("Z (BRCA)","Z (OCA)")[1+brca])
draw.ellipse(0, 0, rlim+0.1, rlim+0.1, col = "black", border = "black")

points(zp[h1],zq[h1],col="red",pch=3,cex=0.5)
points(zp[hp],zq[hp],col="blue",pch=5,cex=0.5)

border_check=rep(1000,length(v1))
for (i in (mvx1-20):(mvx1+20)) 
   border_check[i]=sum(abs((1:length(p) %in% h1) - in.out(cbind(xv1[order(v1)[i],],yv1),cbind(p,q))))
bmin=order(v1)[which.min(border_check)]
lines(-qnorm(xv1[bmin,xsub]/2),-qnorm(yv1[xsub]/2),col="red",lwd=2,lty=1)

border_check=rep(1000,length(v2))
for (i in (mvx2-20):(mvx2+20)) 
   border_check[i]=sum(abs((1:length(p) %in% h2) - in.out(cbind(xv2[order(v2)[i],],yv2),cbind(p,q))))
bmin=order(v2)[which.min(border_check)]
lines(-qnorm(xv2[bmin,xsub]/2),-qnorm(yv2[xsub]/2),col="red",lwd=2,lty=2)



abline(v=-qnorm(max(p[hp])/2),col="blue",lwd=2,lty=3)

legend("bottomleft",c("Obs. Z scores", expression(paste("H"[0]," rej. by p-val")),
       expression(paste("H"[0]," rej. by cFDR"))),bg="white",
       col=c("black","blue","red"),pch=c(1,5,3),pt.cex=c(0.3,0.75,0.75))

if (save_pdf) dev.off()






#### Summary table 
if (brca) table_name="tab_brca.txt" else table_name="tab_oca.txt"

tab=data.frame(gene=gene[h1],tissue=tissue[h1],p=p[h1],q=q[h1],v=vbc1[h1])
write.table(tab,file=paste0(output_dir,table_name),quote=F,row.names=F,col.names=T)






#### Number of discoveries at various thresholds

k=dim(xv1)[1]

hpng=c(); hpnt=c() # gene discoveries, total discoveries, p-value
hcng=c(); hcnt=c() # gene discoveries, total discoveries, cFDR (adjusted)

an=10^-seq(3,10,length.out=100) # values of alpha

for (i in 1:length(an)) {

alpha=an[i]

# Benjamini-Hochberg method on p-values
r=rank(p)
hp=which(r <= max(r[which(p/r < alpha/length(p))]))

# Benjamini-Hochberg method on v-values, adjusted
vbc1=rep(1,length(p)); vbc1[sub]=v1; rvbc1=rank(vbc1)
mvx1=max(rvbc1[which(vbc1/rvbc1 < alpha/length(vbc1))]); mvxl1=which(sub==match(mvx1,rvbc1))
h1=which(rvbc1 <= mvx1)

hpnt=c(hpnt,length(hp))
hpng=c(hpng,length(unique(gene[hp])))

hcnt=c(hcnt,length(h1))
hcng=c(hcng,length(unique(gene[h1])))

}


if (brca) header="BRCA|OCA" else header="OCA|BRCA"
if (brca) plot_name=paste0("twas_fdrcut_brca.pdf") else plot_name=paste0("twas_fdrcut_oca.pdf")

if (save_pdf) pdf(paste0(output_dir,plot_name))

mm=c(hcng/hpng,hcnt/hpnt)
mh=range(mm[which(is.finite(mm))])
plot(0,xlim=range(-log10(an)),type="n",
     ylim=mh,
     xlab=expression(paste("-log"[10],"(FDR) (",alpha,")")),ylab="Ratio of N discoveries",
     main=header);
legend("topleft",c("Genes","Gene-tissue pairs"),lty=c(2,1),
       lwd=2,col=c("red","blue"))
abline(h=1,lwd=2,lty=1)
lines(-log10(an),hcng/hpng,col="red",lwd=2,lty=2)
lines(-log10(an),hcnt/hpnt,col="blue",lwd=2,lty=3)

if (save_pdf) dev.off()

}