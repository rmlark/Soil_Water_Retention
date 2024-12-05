############################################################################################
#
#
#
GALS<-function(nl,vG,Dz,thet_i,dmax,h,plot="F"){
#
#
#
#########################################################################################
#
# h	rainfall matrix : column 1 time (min), depth (cm)
# nl  number of layers
# vG  matrix (nrow=nl,ncol=4), contains thr, ths, alp, nscal  for each layer  
#					alp in kPa
# Dz<-c(5,10)  vector, thickness of layers /cm
# thet_i[nl]  initial vwc for each layer
# dmax<-1000  maximum profile depth (cm)
# plot  whether or not to generate a plot of the model output
#
#
###############################################################################
#
#  References
#
# Guarracino, L. 2007. Estimation of saturated hydraulic conductivity Ks from
# the van Genuchten shape parameter α. Water Resources Research, 43,
# W11502.
#
# Liu, J., Zhang, J., Feng, J. 2008. Green-Ampt model for layered soils with
# nonuniform initial water content under unsteady infiltration. Soil Science
# Society of America Journal, 72, 1041–1047.
#
# Morel-Seytoux, H.J., Meyer, P.D., Nachabe, M., Tourna, J., van Genuchten,
# M.T., Lenhard, R.J. 1996. Parameter equivalence of the Brooks-Corey
# and van Genuchten soil characteristics: preserving the effective capillary
# drive. Water Resources Research, 32, 1251–1258.
#
#########################################################################################

nt<-nrow(h)  # number of rows in rainfall data

H<-matrix(0,nrow=nt,ncol=3) #Time /min, rf /cm

H[,1:2]<-h

for(i in 2:nt){
H[i,3]<-H[(i-1),3]+H[i,2]
}

colnames(H)<-c("Time_min","rainfall","cumulative_rf")

del.h<-vector("numeric",nt)

for(i in 2:nt){
del.h[i]<-(h[i,1]-h[(i-1),1])
}

if(min(del.h[2:nt])==max(del.h[2:nt])){
t.int<-min(del.h[2:nt])	
}else{
stop("Irregular time steps")
}

# t.int is time step in rainfall data (minutes)
#
#############################################################
#
# Soil profile, derived


vG[,3]<-vG[,3]/10			#/10 to turn kPa to cm 

# 1.  Capillary drive, following Morel-Seytoux et al. (1996)

S<-vector("numeric",nl)
m<-1-1/vG[,4]
m2<-m^2
m3<-m^3
S<-(0.046*m+2.07*m2+19.5*m3)/((1+4.7*m+16*m2)*vG[,4])

# 2.  Increase in vwc on saturation from initial

M<-vG[,2]-thet_i

# 3.  Layers

PD<-matrix(nrow=nl,ncol=2)	#prof desc
PD[1,1]<-0
PD[1,2]<-Dz[1]
if(nl>1){
for (i in 2:nl){
PD[i,1]<-PD[(i-1),2]
PD[i,2]<-PD[(i-1),2]+Dz[i]
}
}
PD[nl,2]<-dmax

colnames(PD)<-c("min.Depth","max.Depth")

r.names<-vector("character",nl)
for (i in 1:nl){
r.names[i]<-paste("Layer",i,sep=".")
}
rownames(PD)<-r.names

# 4.  K by layer, following Guarracino (2007) but in units cm/ unit time

time.steps<-60/t.int #time steps per hour

K<-(4.65e4*vG[,2]*vG[,3]*vG[,3])/(24*time.steps)

Dz.b.K<-Dz/K  # ${\rm d}z_i/K_i$ as in term in summation in denominator
		  # on RHS of Liu et al. (2008) Eq 11	 


######################################################################

# State vectors

f<-rep(0,nt) 	#Potential infiltration rate
F<-rep(0,nt) 	#Infiltration
FC<-rep(0,nt) 	#Cumulative infiltration
R<-rep(0,nt)      #Runoff
RC<-rep(0,nt)     #Cumulative runoff
d<-rep(0,nt)	#Depth of wetting front


#####################################################################
r.dep<-min(Dz)/50

for (i in 2:nt){     
depth<-max(r.dep,d[i-1])
layer<-which.layer(depth,PD)
if(layer==1){l.layer<-depth}else{l.layer<-depth-PD[layer,1]}  
		# l.layer is l_{i+1} in Liu et al (2008) Eq[1]
if(layer==1){sum.Dz<-0}else{sum.Dz<-PD[(layer-1),2]}
		# sum.Dz is $\sum_{l=1}^i {\rm d}z_i$ in Liu et al (2008) Eq[10]
if(layer==1){sumDz.b.K<-0}else{sumDz.b.K<-sum(Dz.b.K[1:(layer-1)])}
		# ${\rm d}z_i/K_i$ in Liu et al (2008) Eq[10]
fp<-(sum.Dz+S[layer]+l.layer)/(sumDz.b.K+(l.layer/K[layer])) #cm per 15 min
		#fp computed from Liu et al (2008) Eq[10]
f[i]<-fp
F[i]<-min(H[i,2],fp)		#infiltration
FC[i]<-FC[(i-1)]+F[i]		#cumulative infiltration
R[i]<-max(0,(H[i,2]-fp))	#runoff
RC[i]<-RC[(i-1)]+R[i]		#cumulative runoff
d[i]<-d[(i-1)]+F[i]/(M[layer]) #current depth of wetting front
}

if(plot==T){
dev.new()

par(mfrow=c(2,2))
plot(H[,1],H[,2],type="b",pch=16,col="blue",xlab="Time /min",ylab="Rainfall /cm")

plot(H[,1],H[,2],type="n",xlab="Time /min",ylab="Depth /cm")
points(H[,1],R,type="b",pch=16,col="chocolate4")
points(H[,1],F,type="b",pch=16,col="cyan4")
legend("topright",legend=c("Runoff","Infiltration"),pch=16,col=c("chocolate4","cyan4"),lty=1)


plot(H[,1],H[,3],type="b",pch=16,xlab="Time /min",ylab="Cumulative depth /cm",col="blue",ylim=c(0,(max(H[,3])*1.5)))
points(H[,1],RC,type="b",pch=16,col="chocolate4")
points(H[,1],FC,type="b",pch=16,col="cyan4")
legend("topleft",legend=c("Rainfall","Runoff","Infiltration"),pch=16,col=c("blue","chocolate4","cyan4"),lty=1)

plot(H[,1],-d,type="n",pch=16,xlab="Time /min",ylab="Depth of wetting front /cm",ylim=c(-50,0))
for(i in 1:nt){
lines(c(H[i,1],H[i,1]),c(0,-d[i]),lwd=17,col="cadetblue4",lend=1)
#lines(c(H[i,1],H[i,1]),c(-d[i],min(-d)),lwd=17,col="burlywood3",lend=1)
lines(c(H[i,1],H[i,1]),c(-d[i],-50),lwd=17,col="coral4",lend=1)
}
}

OP<-cbind(H,F,FC,R,RC,d)
return(OP)
}


############################################################################################


which.layer<-function(x,PD){
l<-which((x>=PD[,1])&(x<PD[,2]))
return(l)
}


############################################################################################


initialize.wc<-function(vG,init,nl,smd=0.4){
#
# function to return a vector of nl initial VWCs
#
# init: init rule
#	init = 1 : residual water content
#	init = 2 : field capacity
#	init = 3 : 40% smd
#
#  smd: soil moisture deficit, expressed as a proportion of AWC
#
#


thet_i<-vector("numeric",nl)
# Set initial water content for all layers

if(init==1){
thet_i<-vG[,1]
}else{
for (j in 1:nl){
#
thr<-vG[j,1]
ths<-vG[j,2]
alp<-vG[j,3]
nscal<-vG[j,4]
#
fcap<-swcc(33,thr,ths,alp,nscal)
pwp<-swcc(1500,thr,ths,alp,nscal)
AW<-fcap-pwp
#
if(init==2){
thet_i[j]<-fcap
}else{
thet_i[j]<-pwp+(smd*AW)
}
}
}
return(thet_i)
}

###############################################################################

GALS.summary<-function(OP){

summ<-matrix(ncol=1,nrow=6)
rownames(summ)<-c("Total rainfall /cm","Total runoff /cm",
	"Maximum rainfall rate /cm per 15 min","Maximum runoff rate /cm per 15 min",
	"Runoff /%","Wetted depth /cm")
colnames(summ)<-"Value"
summ[1,1]<-max(OP[,3])
summ[2,1]<-max(OP[,7])
summ[3,1]<-max(OP[,2])
summ[4,1]<-max(OP[,6])
summ[5,1]<-100*round((max(OP[,7])/max(OP[,3])),2)
summ[6,1]<-max(OP[,8])

return(summ)
}
#
#############################################################################
#

GALS.plot<-function(OP,maxdep=50,lpo1="topright",lpo2="topleft"){

H<-OP[,1:3]
F<-OP[,4]
FC<-OP[,5]
R<-OP[,6]
RC<-OP[,7]
d<-OP[,8]
nt<-nrow(OP)

dev.new(height=50,width=60)

par(mfrow=c(2,2))
plot(H[,1],H[,2],type="b",pch=16,col="darkblue",xlab="Time /min",ylab="Rainfall /cm")

plot(H[,1],H[,2],type="n",xlab="Time /min",ylab="Depth /cm")
points(H[,1],R,type="b",pch=16,col="burlywood4")
points(H[,1],F,type="b",pch=16,col="cyan4")
legend(lpo1,legend=c("Runoff","Infiltration"),pch=16,col=c("burlywood4","cyan4"),lty=1)


plot(H[,1],H[,3],type="b",pch=16,xlab="Time /min",ylab="Cumulative depth /cm",col="blue",ylim=c(0,(max(H[,3])*1.5)))
points(H[,1],RC,type="b",pch=16,col="burlywood4")
points(H[,1],FC,type="b",pch=16,col="cyan4")
legend(lpo2,legend=c("Rainfall","Runoff","Infiltration"),pch=16,col=c("darkblue","burlywood4","cyan4"),lty=1)

plot(H[,1],-d,type="n",pch=16,xlab="Time /min",ylab="Depth of wetting front /cm",ylim=c(-maxdep,0))
for(i in 1:nt){
lines(c(H[i,1],H[i,1]),c(0,-d[i]),lwd=17,col="cadetblue4",lend=1)
lines(c(H[i,1],H[i,1]),c(-d[i],-50),lwd=17,col="coral4",lend=1)
}
}


############################################################################################





############################################################################
#
mFDR.alpha.investment<-function(p.values,a_fdr){

#
# Given a sequence of p.values this function tests them with
# marginal false discovery rate control at level a_fdr
# and with alpha investment following the strategy set out
# by Lark (2017). 
#
#	Lark, R.M. 2017. Controlling the marginal false discovery 
#	rate in inferences from a soil data set with α-investment. 
#	European Journal of Soil Science, 68, 221–234. 
#
# The output is a list containing Wealth, a vector with the initial 
# alpha wealth and the amount after successive tests and threshold.pval
# a matrix with the threshold values for mFDR and reported p-values

Nt<-length(p.values)

#########################################################################
#
# Set alpha value and payout (omega)

alp<-a_fdr
omeg<-alp*(1-alp)
hstar<-0
W<-alp

#  Vectors to hold outputs

Wealth<-c(alp,rep(0,Nt))
Level<-rep(0,Nt)
pval<-rep(0,Nt)

################################################
#
# Loop over all tests


for (j in 1:Nt){

p<-p.values[j]
if(j==Nt){final<-"T"}else{final<-"F"}

invest<-test_and_update(p,j,hstar,omeg,alp,W,final)

#invest[1] #accept (0) or reject (1) Ho
#invest[2] #alphaj level of reported test  
W<-invest[3] #Wj  current alpha wealth
hstar<-invest[4] #hstar, index of last rejected hypothesis

Wealth[j+1]<-W
Level[j]<-invest[2]
pval[j]<-p
}

op<-list(Wealth=Wealth,threshold.pval=cbind(Level,p.values))
return(op)
}

############################################################################
#
SWRC.mFDR.plot<-function(op,par.order,pos="bottomright"){

#
# Plot output of mFDR control on a sequence of tests
# for adding successive SWRC parameters allowed to 
# differ between groups.
#
# op is the output list from mFDR.alpha.investment
# par.order is a vector giving the index for the 
# successive parameters in the standard order 
# thr,ths,alp,nscal
#
# so if the sequence is alp, nscal, ths, thr then 
# order is c(3,4,2,1).  This ensures correct symbols
# on x axis
#
# pos is position for the legend, following the base R
# legend function (default value is "bottom.right","none"
# means no legend to be drawn.
#

Level<-op$threshold.pval[,1]
pval<-op$threshold.pval[,2]

par(mar=c(4,5,3,3))
Nt<-length(Level)
Test<-seq(0,Nt)
Testl<-seq(1,Nt)
plot(Testl,Level,xlim=c(1,Nt),ylim=c(0,0.08),pch=16,ylab="Probability",xlab="",
yaxt="n",xaxt="n",cex.lab=1.1)
mtext("Test",1,0.5,line=2)

alpha<-expression(paste(alpha))
theta_s<-expression(paste(theta[s]))
theta_r<-expression(paste(theta[r]))

test.lab<-c(theta_r,theta_s,alpha,"n")

axis(1,seq(1,4),test.lab[par.order])
axis(2,c(0,0.01,0.03,0.05,0.08),label=c("0","0.01","0.03","0.05",">0.05"),
cex.axis=0.95,las=2,mgp=c(3,0.75,0))
lines(Testl,Level)
axis.break(2,0.065,style="slash")

plotpval<-pval
for (i in 1:Nt){
if(pval[i]>0.05)plotpval[i]<-0.08
}
points(Testl,plotpval)
lines(Testl,plotpval,lty=3)
if(pos!="none"){
legend(pos,legend=c("FDR threshold","p-value"),pch=c(16,1),lty=c(1,3))
}
}


###########################################################################

test_and_update<-function(p,j,hstar,omeg,alp,W,final){

#
#  Based on Foster & Stine (2008).  JRSS(B) 70, 429--444.
#
# p is p-value for current (jth) test
# hstar is index of last rejected hypothesis
# om is payout
# W is alpha wealth after (j-1)th test (on function call)


#first, compute alpha_j

# This uses the rule in Equation (15) of Foster and Stine (2008)

if(final=="F"){
alphaj<-W/(1+j-hstar)
} else {
alphaj<-W
}

# conservative hopeful thrifty strategy

alphaj<-min(alphaj,alp)


# Now, test p-value against alphaj

if(p<=alphaj){

# null rejected

Rj<-1
hstar<-j
Wj<-W+omeg

} else {

Rj<-0
Wj<-W-(alphaj/(1-alphaj))

}

# return outputs: value of Rj, alphaj, 
# current alpha wealth, hstar (possibly updated) 

# Rj 1 if null is rejected, 0 if accepted
# alphaj level of reported test
# Wj  current alpha wealth
# hstar, index of last rejected hypothesis

return(c(Rj,alphaj,Wj,hstar))
}


############################################################################################
############################################################################################
#
#  Soil physics:  methods to plot WRC data, fit and plot Van Genuchten models and to 
#   interpret the parameters
#
############################################################################################
############################################################################################



swcc<-function(h,thr,ths,alp,nscal)thr+(ths-thr)/(1+(alp*h)^nscal)^(1-(1/nscal))


VanGenuchten.fit.single<- function(data.df,init.vals,roundoff=4){
  
  theta<-data.df$theta
  h<-data.df$h

  swcc<-function(h,thr,ths,alp,nscal)thr+(ths-thr)/(1+(alp*h)^nscal)^(1-(1/nscal))
  
  nlmvwc<-nls(theta~swcc(h,thr,ths,alp,nscal), data=data.df, start=c(thr=init.vals[1],
  ths=init.vals[2],alp=init.vals[3],nscal=init.vals[4]))
  op<-list("Coefficient.estimates"=as.matrix(round(coef(nlmvwc),roundoff)))
  return(op)
  
}



VanGenuchten.fit.compare<-function(data.df,init.vals,g.name,par.var.null,
par.var.full,cov.name,i.seed){

#Function to fit VanGenuchten WRC model to data where one or more
#treatments, or other groups, can be compared with respect to one or
#more of the four parameters.  Two models are fitted, a null model with
#some number of parameters common to all groups, and an alternative model
#in which a larger subset of the parameters differ between the groups.

#The coefficients for each group in the alternative model (i.e. with most 
#parameters differing between the groups) are exported as are the 
#log-likelihood ratio for the comparison between the alternative and the null
#models.

#data.df data frame with tensions (kPa) in variable h, VWC (proportion) in
#  variable theta, a factor (g.name) which groups all observations from the
#  same soil sample (or the mean from a single plot), and a factor (cov.name)
#  which specifies the groups to be compares (e.g. treatments)
#
#init.vals a vector with starting point for parameters, in order "thr","ths","alp","nscal"
#par.var.null a vector coding which parameters (same order as init.vals) are constant (0)
# and which vary between the levels of the covariate in the null model.
#par.var.full  a vector coding which parameters (same order as init.vals) are constant (0)
# and which vary between the levels of the covariate in the alternative model.

#g.name  factor which groups observations from a single specimen (or plot average)
#cov.name factor specifying the treatments or other groups to be compared
#


library(saemix)
dfname<-deparse(substitute(data.df))

g.col<-which(names(data.df)==g.name) #which column of data frame contains groupings
ng<-nlevels(data.df[,g.col])

c.col<-which(names(data.df)==cov.name) #which column of data frame contains covariate
Ng<-nlevels(data.df[,c.col])  #nlevels of covariate



N=nrow(data.df)
Group=vector("numeric",N)
Group[1]<-1

for (i in 2:N){
if(data.df[i,g.col]==data.df[(i-1),g.col]){
Group[i]<-Group[(i-1)]}else{
Group[i]<-Group[(i-1)]+1}
}
data.df[,g.col]<-factor(Group)

saemix.options<-saemixControl(nbiter.saemix = c(400, 200),
  nbiter.sa = NA,
  nb.chains = 1,
  fix.seed = TRUE,
  seed = i.seed,
  nmc.is = 5000,
  nu.is = 4,
  nbdisplay = 100,
  displayProgress=F,print=F,save=F,save.graphs=F,print.is=F) 


wrc.data<-saemixData(data.df,name.group=g.name,
name.predictors="h",name.response="theta",name.covariates=cov.name,verbose=F)

vanGenuchten_null<-saemixModel(model=vanG,
description="Van Genuchten water release curve",
psi0=matrix(init.vals,ncol=4,byrow=T,dimnames=
list(NULL,c("thr","ths","alp","nscal"))),name.response="theta",
covariate.model=matrix(par.var.null,ncol=4,byrow=T),
transform.par=c(3,3,1,1),verbose=F)

pdf(file = NULL)
wrc_null.fit<-saemix(vanGenuchten_null,wrc.data,control=saemix.options)
dev.off()

vanGenuchten_groups<-saemixModel(model=vanG,
description="Van Genuchten water release curve",
psi0=matrix(init.vals,ncol=4,byrow=T,dimnames=
list(NULL,c("thr","ths","alp","nscal"))),name.response="theta",
covariate.model=matrix(par.var.full,ncol=4,byrow=T),
transform.par=c(3,3,1,1),verbose=F)

pdf(file = NULL)
wrc_groups.fit<-saemix(vanGenuchten_groups,wrc.data,control=saemix.options)
dev.off()

ldf<-sum(par.var.full)-sum(par.var.null)

coeffs<-matrix(nrow=4,ncol=Ng)
rownames(coeffs)<-colnames((coef(wrc_groups.fit))$population$psi)
colnames(coeffs)<-levels(data.df[,which(names(data.df)==cov.name)])

reftab<-table(data.df[,g.col],data.df[,c.col])

for (i in 1:Ng){ #over levels of factor
row.L<-min(which(reftab[,i]!=0))
#index.factor.level<-min(which(data.df[,c.col]==(levels(data.df[,c.col]))[i]))
#rep.level<-as.numeric(data.df[index.factor.level,g.col])
coeffs[,i]<-((coef(wrc_groups.fit))$population$psi)[row.L,]
}

par.index<-which(summary(wrc_groups.fit@results,print=F)$fixed.effects$p!="-")
pars<-(summary(wrc_groups.fit@results,print=F)$fixed.effects$Parameter[par.index])
par.est<-as.numeric(summary(wrc_groups.fit@results,print=F)$fixed.effects$Estimate[par.index])
par.se<-as.numeric(summary(wrc_groups.fit@results,print=F)$fixed.effects$SE[par.index])
p.val<-as.numeric(summary(wrc_groups.fit@results,print=F)$fixed.effects$p[par.index])

inference<-list("Parameter effects"=pars,"Estimates"=par.est,"Standard error"=par.se,"P.value"=p.val)
COMP<-cbind(par.var.null,par.var.full)
rownames(COMP)<-c("thr","ths","alp","nscal")
colnames(COMP)<-c("Group.dependent.parameters.null","Group.dependent.parameters.full")
op<-list("Coefficient.estimates"=coeffs,"Comparison"=COMP,"Inference"=inference)

return(op)
}

############################################################################################
VanGenuchten.fit.group<-function(data.df,init.vals,g.name){


#Function to fit VanGenuchten WRC model to a replicated set of observations
# from a single group or treatment.

#The coefficients for each group in the alternative model (i.e. with most 
#parameters differing between the groups) are exported as are the 
#log-likelihood ratio for the comparison between the alternative and the null
#models.

#data.df data frame with tensions (kPa) in variable h, VWC (proportion) in
#  variable theta, a factor (g.name) which groups all observations from the
#  same soil sample (or the mean from a single plot).

#init.vals a vector with starting point for parameters, in order "thr","ths","alp","nscal"
#par.var.null a vector coding which parameters (same order as init.vals) are constant (0)
# and which vary between the levels of the covariate in the null model.
#par.var.full  a vector coding which parameters (same order as init.vals) are constant (0)
# and which vary between the levels of the covariate in the null model.

#g.name  factor which groups observations from a single specimen (or plot average)
#cov.name factor specifying the treatments or other groups to be compared
#


library(saemix)
dfname<-deparse(substitute(data.df))

g.col<-which(names(data.df)==g.name) #which column of data frame contains groupings
ng<-nlevels(data.df[,g.col])

saemix.options<-saemixControl(displayProgress=F,print=F,save=F,save.graphs=F,print.is=F) 


wrc.data<-saemixData(data.df,name.group=g.name,
name.predictors="h",name.response="theta",verbose=F)

vanGenuchten<-saemixModel(model=vanG,
description="Van Genuchten water release curve",
psi0=matrix(init.vals,ncol=4,byrow=T,dimnames=
list(NULL,c("thr","ths","alp","nscal"))),name.response="theta",
transform.par=c(3,3,1,1),verbose=F)

pdf(file = NULL)
wrc.fit<-saemix(vanGenuchten,wrc.data,control=saemix.options)
dev.off()


coeffs<-matrix(nrow=4,ncol=1)
sterr<-matrix(nrow=4,ncol=1)
rownames(coeffs)<-wrc.fit@results@name.fixed
rownames(sterr)<-wrc.fit@results@name.fixed
colnames(coeffs)<-c("Estimate")
colnames(sterr)<-c("Standard error")
coeffs[,1]<-wrc.fit@results@fixed.effects
sterr[,1]<-wrc.fit@results@se.fixed

op<-list("Coefficient.estimates"=coeffs,"Standard.errors"=sterr)

return(op)
}

############################################################################################


VanGenuchten.fit.compare<-function(data.df,init.vals,g.name,par.var.null,
par.var.full,cov.name){

#Function to fit VanGenuchten WRC model to data where one or more
#treatments, or other groups, can be compared with respect to one or
#more of the four parameters.  Two models are fitted, a null model with
#some number of parameters common to all groups, and an alternative model
#in which a larger subset of the parameters differ between the groups.

#The coefficients for each group in the alternative model (i.e. with most 
#parameters differing between the groups) are exported as are the 
#log-likelihood ratio for the comparison between the alternative and the null
#models.

#data.df data frame with tensions (kPa) in variable h, VWC (proportion) in
#  variable theta, a factor (g.name) which groups all observations from the
#  same soil sample (or the mean from a single plot), and a factor (cov.name)
#  which specifies the groups to be compares (e.g. treatments)
#
#init.vals a vector with starting point for parameters, in order "thr","ths","alp","nscal"
#par.var.null a vector coding which parameters (same order as init.vals) are constant (0)
# and which vary between the levels of the covariate in the null model.
#par.var.full  a vector coding which parameters (same order as init.vals) are constant (0)
# and which vary between the levels of the covariate in the alternative model.

#g.name  factor which groups observations from a single specimen (or plot average)
#cov.name factor specifying the treatments or other groups to be compared
#


library(saemix)
dfname<-deparse(substitute(data.df))

g.col<-which(names(data.df)==g.name) #which column of data frame contains groupings
ng<-nlevels(data.df[,g.col])

c.col<-which(names(data.df)==cov.name) #which column of data frame contains covariate
Ng<-nlevels(data.df[,c.col])  #nlevels of covariate

i.seed<-sample(10000000:100000000,1)

N=nrow(data.df)
Group=vector("numeric",N)
Group[1]<-1

for (i in 2:N){
if(data.df[i,g.col]==data.df[(i-1),g.col]){
Group[i]<-Group[(i-1)]}else{
Group[i]<-Group[(i-1)]+1}
}
data.df[,g.col]<-factor(Group)

saemix.options<-saemixControl(nbiter.saemix = c(400, 200),
  nbiter.sa = NA,
  nb.chains = 1,
  fix.seed = TRUE,
  seed = i.seed,
  nmc.is = 5000,
  nu.is = 4,
  nbdisplay = 100,
  displayProgress=F,print=F,save=F,save.graphs=F,print.is=F) 


wrc.data<-saemixData(data.df,name.group=g.name,
name.predictors="h",name.response="theta",name.covariates=cov.name,verbose=F)

vanGenuchten_null<-saemixModel(model=vanG,
description="Van Genuchten water release curve",
psi0=matrix(init.vals,ncol=4,byrow=T,dimnames=
list(NULL,c("thr","ths","alp","nscal"))),name.response="theta",
covariate.model=matrix(par.var.null,ncol=4,byrow=T),
transform.par=c(3,3,1,1),verbose=F)

pdf(file = NULL)
wrc_null.fit<-saemix(vanGenuchten_null,wrc.data,control=saemix.options)
dev.off()

vanGenuchten_groups<-saemixModel(model=vanG,
description="Van Genuchten water release curve",
psi0=matrix(init.vals,ncol=4,byrow=T,dimnames=
list(NULL,c("thr","ths","alp","nscal"))),name.response="theta",
covariate.model=matrix(par.var.full,ncol=4,byrow=T),
transform.par=c(3,3,1,1),verbose=F)

pdf(file = NULL)
wrc_groups.fit<-saemix(vanGenuchten_groups,wrc.data,control=saemix.options)
dev.off()

ldf<-sum(par.var.full)-sum(par.var.null)

coeffs<-matrix(nrow=4,ncol=Ng)
rownames(coeffs)<-colnames((coef(wrc_groups.fit))$population$psi)
colnames(coeffs)<-levels(data.df[,which(names(data.df)==cov.name)])

reftab<-table(data.df[,g.col],data.df[,c.col])

for (i in 1:Ng){ #over levels of factor
row.L<-min(which(reftab[,i]!=0))
coeffs[,i]<-((coef(wrc_groups.fit))$population$psi)[row.L,]
}




L<-2*(wrc_groups.fit@results@ll.is-wrc_null.fit@results@ll.is)
p.val<-1-pchisq(L,ldf)

inference<-list("L"=L,"P.value"=p.val)
COMP<-cbind(par.var.null,par.var.full)
rownames(COMP)<-c("thr","ths","alp","nscal")
colnames(COMP)<-c("Group.dependent.parameters.null","Group.dependent.parameters.full")
op<-list("Coefficient.estimates"=coeffs,"Comparison"=COMP,"Inference"=inference)

return(op)
}

############################################################################################


trform<-function(x,tr,direc){
#
#  transform parameters based on saemix conventions
#  tr = 0 (none), 1 (log), 2 (probit), 3 (logit)
#  direc fw (forward transform) or bw (backtransform)
#

if(tr==0){y=x}
if(tr==1){y=tr.log(x,direc)}
if(tr==3){y=tr.logit(x,direc)}

return(y)
}


tr.log<-function(x,direc){
if (direc=="fw"){y<-log(x)}else{y<-exp(x)}
return(y)
}

tr.logit<-function(x,direc){
if (direc=="fw"){
y<-log(x/(1-x))
}else{
e.x<-exp(x)
y<-e.x/(1+e.x)
}
return(y)
}




############################################################################################

vanG<-function(ps1,id,xidep){
h<-xidep[,1]
thr<-ps1[id,1]
ths<-ps1[id,2]
alp<-ps1[id,3]
nscal<-ps1[id,4]
mscal<-mscal<-1-1/nscal
den<-(1+((abs(alp*h))^nscal))^mscal
return(thr+((ths-thr)/den))
}


############################################################################################

vanGdraw<-function(psi,h){
thr<-psi[1]
ths<-psi[2]
alp<-psi[3]
nscal<-psi[4]
mscal<-1-(1/nscal)
den<-(1+((abs(alp*h))^nscal))^mscal
vwc<-thr+((ths-thr)/den)
return(as.numeric(vwc))
}

############################################################################################

plot.wrc.data<-function(data.df,coeffs,hvar="h",thetavar="theta",
xlab="Tension /kPa",maxth=-1,groups="none",main=""){
if(missing(coeffs)){models=FALSE}else{models=TRUE}

#Function to plot WRC data from dataframe data.df
#hvar optional name for tension variable in data.df, default "h"
#thetavar optional name for vwc in data.df, default "theta"
#xlab optional label for tension axis, default "Tension /kPa",
#maxth optional maximum vwc for plot 
#groups optional name for groups to be distinguished on plot
#main  optional main title for plot

library(RColorBrewer)


mypalette<-brewer.pal(9,"Set1")
n<-nrow(data.df)

hcol=which(names(data.df)==hvar)
thcol=which(names(data.df)==thetavar)
if(groups!="none"){
gcol=which(names(data.df)==groups)
gnames<-as.character(levels(data.df[,gcol]))
ngr<-length(gnames)
cols<-mypalette[as.numeric((data.df[,gcol]))]
}else{
cols<-rep("black",n)
}
if(maxth<0){
ylim=c(0,(max(data.df[,thcol])*1.01))}else{
ylim=c(0,maxth)
}

options("scipen"=100, "digits"=4) 
plot(data.df[,hcol],data.df[,thcol],ylim=ylim,xlab=xlab,ylab=expression(paste(theta)),
pch=16,col=cols,log="x",main=main)
if(groups!="none"){
legend("topright", legend = gnames,
pch = 16, col = mypalette)
}

if(models==TRUE){
H<-c(seq(0.05,1,0.001),seq(1,1500,0.1))
if(groups=="none"){
  thdr<-vanGdraw(coeffs[,1],H)
  lines(H,thdr,lwd=2)
  }else{
for (i in 1:ngr){
  thdr<-vanGdraw(coeffs[,i],H)
  lines(H,thdr,lwd=2,col=mypalette[i])
}
}
}

}

############################################################################################

SQ.indices<-function(coeffs,bulk_density,t_m,t_FC, t_PWP){

# coeffs: coefficients for m soils in wrc coefficient format
# i.e. a 4 x m matrix in which the rows correspond respectively
# to thr, ths, alp, nscal
#
# bulk_density an optional argument.  A m-vector of bulk density values
# (g cm^-3).  If these are not provided then the value is inferred from
# ths, assuming a particle density of 2.65 g cm^-3 (this may underestimate
# particle density for very ferruginous soils, and overestimate it for soils
# with a large organic content.  (Landon, J.R. (ed), 1984.  Booker Tropical
# Soil Manual (1st Edn) Longman, Harlow.  page 97, section 6.5).  
# Required to compute Dexter's S value.
#
# t_m tension (kPa) at which VWC ~ matrix porosity.  
#  defaulting to 4.9 kPa following Reynolds et al. (2007) (page 318, 
#  max. equiv. pore diameter of 0.06 mm.
#
# t_FC tension (kPa) for field capacity.  An optional argument. Default value
# is 9.8 following Reynolds
# t_PWP tension (kPa) for permanent wilting point.  An optional argument. 
# Default 1471 following Reynolds



ncl<-ncol(coeffs)
if (missing(bulk_density)){
pd<-2.65
bulk_density<-pd*(1-coeffs[2,])
}

if(missing(t_m)){
t_m=4.9
}

if(missing(t_FC)){
t_FC=9.8
}

if(missing(t_PWP)){
t_PWP=1471
}



#if(t_FC==9.8){
#cat(paste("\n\tField Capacity is set at 9.8 kPa by default.\n"))
#cat(paste("\tThis follows Reynolds et al (2007), the source\n"))
#cat(paste("\tof these indices, but 33 kPa is more usually used.\n\n"))
#cat(paste("\tTo set FC to 33 kPa, include the optional argument:\n"))
#cat(paste("\tt_FC=33\n\n\n"))
#}

noutp<-matrix(nrow=ncl,ncol=6)
rownames(noutp)<-colnames(coeffs)
colnames(noutp)<-c("Dexter.S","Total_Porosity","Macropores","Relative_Water_Capacity",
            "Plant_Available_Water_Capacity","Air_Capacity")

interp<-matrix(" ",nrow=ncl,ncol=6)
rownames(interp)<-colnames(coeffs)
colnames(interp)<-c("Dexter.S","Total_Porosity","Macropores","Relative_Water_Capacity",
            "Plant_Available_Water_Capacity","Air_Capacity")


for(i in 1:ncl){


thr<-coeffs[1,i]
ths<-coeffs[2,i]
alp<-coeffs[3,i]
nscal<-coeffs[4,i]
bd<-bulk_density[i]

# i.  Dexter S

thr_g<-thr/bd #convert volumetric WC to gravimetric (see Dexter (2004) 
ths_g<-ths/bd # section 6 paragraph 2

noutp[i,1]<-round((nscal)*(ths_g-thr_g)*((((2*nscal)-1)/(nscal-1))^((1/nscal)-2)),3)

Inter<-"Poor microstructural quality"
if(noutp[i,1]>0.035){Inter<-"Good microstructural quality"}
if(noutp[i,1]<=0.02){Inter<-"Very poor microstructural quality"}

interp[i,1]<-Inter
# ii.  Total porosity

noutp[i,2]<-ths
interp[i,2]<-"~"
# iii. Macroporosity following Reynolds et al, 2007 

thm<-vanGdraw(coeffs[,i],t_m) #matrix porosity
noutp[i,3]<-round((ths-thm),2)

if(noutp[i,3]<=0.04){Inter<-"Macroporosity degraded"
}else{Inter<-"Undegraded (medium to fine textured)"}

interp[i,3]<-Inter

# iv. Relative water capacity Reynolds et al, 2007 

thFC<-vanGdraw(coeffs[,i],t_FC)
noutp[i,4]<-round((thFC/ths),1)

Inter<-"RWC in optimal range for microbial activity"
if(noutp[i,4]<=0.6){Inter<-"RWC suboptimal for microbial activity"}
if(noutp[i,4]>0.7){Inter<-"RWC exceeds range for optimal microbial activity"}


interp[i,4]<-Inter

# v. Plant available water capacity 

thPWP<-vanGdraw(coeffs[,i],t_PWP)
noutp[i,5]<-round(thFC-thPWP,2)
if(noutp[i,5]>0.2){Inter<-"PAW ideal"}
if(noutp[i,5]<=0.2){Inter<-"PAW good"}
if(noutp[i,5]<=0.15){Inter<-"PAW limited"}
if(noutp[i,5]<=0.10){Inter<-"PAW poor"}


interp[i,5]<-Inter


#print(c(i,thFC,thPWP))
#vi.  Air capacity

noutp[i,6]<-round(ths-thFC,2)

Inter<-"AC adequate, but marginal in fine-textured soil"
if(noutp[i,6]>0.15){Inter<-"AC good for all soils"}
if(noutp[i,6]<=0.1){Inter<-"Aeration deficit likely in rooting zone"}

interp[i,6]<-Inter

} #end of loop over classes


return(list("Indices"=noutp,"Interpretation"=interp))
}

############################################################################################

print.indices<-function(indices,ret=F){
roundit<-c(3,2,2,2,2,2)
indices$Indices<-round(indices$Indices,roundit)
nclass<-nrow(indices$Indices)
if(nclass>1){
outp<-data.frame(matrix(nrow=(nclass*6),ncol=2))
colnames(outp)<-c("Index value","Interpretation")
ronames<-c("Dexter's S","Total porosity","Macroporosity","Relative water capacity",
"Plant-available water capacity","Air capacity")
long.names<-c("Dexter's S","Total porosity","Macroporosity","Relative\n water capacity",
"Plant-available\n water capacity","Air capacity")
class.names<-rownames(indices$Indices)

dfro<-0
for (icl in 1:nclass){
for (i in 1:6){
dfro<-dfro+1
rownames(outp)[dfro]<-paste(class.names[icl],ronames[i],sep=", ")
cat(paste(long.names[i],class.names[icl],indices$Indices[icl,i],indices$Interpretation[icl,i],sep='\t','\n'))
outp[dfro,1]<-indices$Indices[icl,i]
outp[dfro,2]<-indices$Interpretation[icl,i]
}
}
}else{
outp<-data.frame(matrix(nrow=(nclass*6),ncol=2))
colnames(outp)<-c("Index value","Interpretation")
ronames<-c("Dexter's S","Total porosity","Macroporosity","Relative water capacity",
"Plant-available water capacity","Air capacity")
long.names<-c("Dexter's S","Total porosity","Macroporosity","Relative\n water capacity",
"Plant-available\n water capacity","Air capacity")


for (i in 1:6){
rownames(outp)[i]<-ronames[i]
cat(paste(long.names[i],indices$Indices[i],indices$Interpretation[i],sep='\t\t','\n'))
outp[i,1]<-indices$Indices[i]
outp[i,2]<-indices$Interpretation[i]
}
}



if(ret){return(outp)}
}

############################################################################################

SQ.BSt.indices<-function(coeffs,bulk_density,t_m,t_FC, t_PWP){

n.sam<-nrow(coeffs)
n.gr<-ncol(coeffs)/4

index.names<-c("S","Total_por","Macropores","RWC","PAWC","AC")
cnames<-colnames(coeffs)
names<-vector("character",n.gr)

ico<-0
for(i in seq(1,(1+(n.gr-1)*4),4)){
ico<-ico+1
names[ico]<-strsplit(cnames[i],"thr.")[[1]][2]
}

output<-matrix(nrow=n.sam,ncol=(n.gr*6))
op.names<-vector("character",1)

for(i in 1:n.gr){
op.names<-c(op.names,paste(names[i],index.names,sep="."))
}
colnames(output)<-op.names[-1]

for (it in 1:n.sam){
in.put<-matrix(coeffs[it,],nrow=4)
colnames(in.put)<-names
index.vals<-SQ.indices(in.put)$Indices
ind.v<-vector("numeric",1)
for (j in 1:n.gr)
ind.v<-c(ind.v,index.vals[j,])
output[it,]<-ind.v[-1]
}
return(output)
}

############################################################################################


########################################################################
plot.pvd<-function(coeffs,Title="",Plot=T){

H<-c(seq(0.05,1,0.001),seq(1,1500,0.1))
EPD<-300/H
H.shind<-intersect(which(EPD>0.5),which(EPD<1000))
H.sh<-H[H.shind]

psd<-unlist(lapply(H.sh,Sh,coeffs))

if(Plot){
options("scipen"=100, "digits"=4) 
ylab="Pore volume density"
xlab=expression(paste("Equivalent pore diameter /",mu,"m"))

plot(H.sh,psd,type="l",log="x",ylab=ylab,xlab=xlab,
main=Title)
}
op<-cbind(H.sh,psd)
return(op)
}

#
#

Sh<-function(h,psi){
tr<-psi[1]
ts<-psi[2]
a<-psi[3]
n<-psi[4]
m<-1-(1/n)
SH<-(m*n*(ts-tr)*(a^n)*(h^n)*(1+(a*h)^n)^(-(m+1)))#change of sign as h here stands for matric potential
return(SH)
}

#

############################################################################################


############################################################################################

VanGenuchten.fit.compare.bs<-function(data.df,init.vals,g.name,
par.var.full,cov.name,i.seed,nboot,bsmeth="parametric"){

#Function to fit VanGenuchten WRC model to data where one or more
#treatments, or other groups, can be compared with respect to one or
#more of the four parameters.  The function returns bootstrap samples of
#the model parameters

#data.df data frame with tensions (kPa) in variable h, VWC (proportion) in
#  variable theta, a factor (g.name) which groups all observations from the
#  same soil sample (or the mean from a single plot), and a factor (cov.name)
#  which specifies the groups to be compares (e.g. treatments)
#
#init.vals a vector with starting point for parameters
#par.var.full a vector coding which parameters, in order "thr","ths","alp","nscal",
#are constant (0) and which vary between the levels of the covariate.

#g.name  factor which groups observations from a single specimen (or plot average)
#cov.name factor specifying the treatments or other groups to be compared
#i.seed a seed value for the random number generator
#nboot  number of bootstrap samples to draw for each parameter
#bsmeth  method for bootstrapping.  Default is "parametric".  There are other options
#for the saemix.bootstrap function, "conditional" is one.  Others less useful here.


library(saemix)
dfname<-deparse(substitute(data.df))

g.col<-which(names(data.df)==g.name) #which column of data frame contains groupings
ng<-nlevels(data.df[,g.col])

c.col<-which(names(data.df)==cov.name) #which column of data frame contains covariate
Ng<-nlevels(data.df[,c.col])  #nlevels of covariate
lev.names<-levels(data.df[,c.col]) #level names


N=nrow(data.df)
Group=vector("numeric",N)
Group[1]<-1

for (i in 2:N){
if(data.df[i,g.col]==data.df[(i-1),g.col]){
Group[i]<-Group[(i-1)]}else{
Group[i]<-Group[(i-1)]+1}
}
data.df[,g.col]<-factor(Group)

saemix.options<-saemixControl(nbiter.saemix = c(400, 200),
  nbiter.sa = NA,
  nb.chains = 1,
  fix.seed = TRUE,
  seed = i.seed,
  nmc.is = 5000,
  nu.is = 4,
  nbdisplay = 100,
  displayProgress=F,print=F,save=F,save.graphs=F,print.is=F) 


wrc.data<-saemixData(data.df,name.group=g.name,
name.predictors="h",name.response="theta",name.covariates=cov.name,verbose=F)

vanGenuchten_groups<-saemixModel(model=vanG,
description="Van Genuchten water release curve",
psi0=matrix(init.vals,ncol=4,byrow=T,dimnames=
list(NULL,c("thr","ths","alp","nscal"))),name.response="theta",
covariate.model=matrix(par.var.full,ncol=4,byrow=T),
transform.par=c(3,3,1,1),verbose=F)

pdf(file = NULL)
wrc_groups.fit<-saemix(vanGenuchten_groups,wrc.data,control=saemix.options)
dev.off()

bs<-saemix.bootstrap(wrc_groups.fit,method=bsmeth,nboot=nboot)


Ng #number of levels in categorical covariate

#Remove first column of bs 
bs<-bs[,-1]
#Remove omega values
npar<-4+sum(par.var.full)
bs<-bs[,1:npar]
par.sets<-matrix(nrow=nboot,ncol=(4*Ng))
colnames(par.sets)<-1:(4*Ng)
o.co<-0
for(i in 1:Ng){
for(j in 1:4){
o.co<-o.co+1
colnames(par.sets)[o.co]<-paste(c("thr","ths","alp","nscal")[j],lev.names[i],sep=".")
}
}

cpvf<-vector("numeric",4)
cpvf[1]<-par.var.full[1]
for(i in 2:4){
cpvf[i]<-cpvf[(i-1)]+par.var.full[i]
}

cols.lev1<-cpvf+1:4-par.var.full
cols.lev2<-vector("numeric",4)
cols.lev2<-cols.lev1*(par.var.full-1)^2

for(i in 1:4){
par.sets[,i]<-bs[,cols.lev1[i]]
}

for(i in 1:4){
if(par.var.full[i]==0){
par.sets[,(i+4)]<-bs[,cols.lev2[i]]
}else{
k<-cols.lev1[i]
tr<-c(3,3,1,1)[i]
par.sets[,(i+4)]<-
trform((trform(bs[,k],tr,"fw")+bs[,(k+1)]),tr,"bw")
}
}

return(par.sets)
}

###############################





