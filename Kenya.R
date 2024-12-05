#########################################################################################
#
#
#  Demonstration script for use of R codes in the CEPHaSoilPhys_1.R file available from
#  LINK HERE
#
#  Fitting a van Genuchten SWRC model to data on volumetric water content of six soil
#  samples from Kenya, which belong to two World Reference Base (1998) classes: 
#  Nitisols and Solonetz. Inference about differences with respect to van Genuchten
#  parameters, computation of soil quality indices, plotting SWRC and pose size density
#  plots.  Parametric bootstrap sampling of the model parameters. Green Ampt model runs
#  using the rainfall data as presented in the paper.
#  
#
#  Data from WOSIS (Batjes et al., 2017) are in the file Kenya_water2.csv
#  All under licence CC-BY-NC https://creativecommons.org/licenses/by-nc/4.0/deed.en
#
#
#  Method Undisturbed soil sample in metal/PVC-ring, saturated from field water content
#  then desorption of water to target tension either (0.1 to 50 kPa) by a hanging water
#  column or (250, 1500 kPa) on a pressure plate
#
#
#  R.M. Lark  29th November 2024
#
#
#  References
#
# 	Batjes, N.H., Ribeiro, E., van Oostrum, A., Leenaars, J., de Mendes, J.
#	2016. Standardised soil profile data for the world (WoSIS, July 2016
#	snapshot),http://dx.doi.org/10.17027/isric-wdcsoils.20160003
#
#  	World Reference Base 1998. World Reference Base on Soils. Rome, Italy:
#	FAO-ISRIC. 	
#
#########################################################################################
#


source("CEPHaSoilPhys_1.R")


library(saemix)
library(plotrix)
library(TeachingDemos)


data.df<-read.csv("Kenya_water2.csv",header=T,stringsAsFactors=T)
head(data.df)

#  data.df data frame with tensions (kPa) in variable h, VWC (proportion) in
#  variable theta, a factor (g.name) which groups all observations from the
#  same soil sample (or the mean from a single plot), and a factor (cov.name)
#  which specifies the groups to be compares (e.g. treatments)

plot.wrc.data(data.df,hvar="h",thetavar="theta",
xlab="Tension /kPa",maxth=-1,groups="WRB",main="Kenyan soils")

###############################################################################		
#													#
#   All van Genuchten parameters are given in the order thr,ths,alp,nscal.	#
#													#
###############################################################################

init.vals<-c(0.1,0.5,0.2,0.4)  # parameter values to initiate the algorithm

# First model run

par.var.null<-c(0,0,0,0)  #Which parameters differ between groups in null model? (None here)
par.var.full<-c(0,1,0,0)  #Which parameters differ between groups in alternative model? (ths Here) 

op.1<-VanGenuchten.fit.compare(data.df,init.vals,g.name="Profile_ID",par.var.null,
par.var.full,cov.name="WRB")

print(op.1$Inference) #Evidence to reject a null hypothesis of common values of ths
p.1<-op.1$Inference$P.value

# Second model run

par.var.null<-c(0,1,0,0)  #Common and fixed parameters as in model op.1
par.var.full<-c(1,1,0,0)  # thr and ths differ between soil classes

op.2<-VanGenuchten.fit.compare(data.df,init.vals,g.name="Profile_ID",par.var.null,
par.var.full,cov.name="WRB")

print(op.2$Inference) #Evidence to reject a null hypothesis of common values of ths
p.2<-op.2$Inference$P.value

# Third model run

par.var.null<-c(1,1,0,0)  #Common and fixed parameters as in model op.2
par.var.full<-c(1,1,1,0)  # thr, ths and alp differ between soil classes

op.3<-VanGenuchten.fit.compare(data.df,init.vals,g.name="Profile_ID",par.var.null,
par.var.full,cov.name="WRB")

print(op.3$Inference) #No evidence to reject a null hypothesis of common values of alp
p.3<-op.3$Inference$P.value


# Fourth model run

par.var.null<-c(1,1,0,0)  #Common and fixed parameters as in model op.2
par.var.full<-c(1,1,0,1)  # thr, ths and nscal differ between soil classes

op.4<-VanGenuchten.fit.compare(data.df,init.vals,g.name="Profile_ID",par.var.null,
par.var.full,cov.name="WRB")

print(op.4$Inference) #No evidence to reject a null hypothesis of common values of nscal
p.4<-op.4$Inference$P.value


#  Plot the p-values and thresholds for marginal false discovery rate control with alpha investment

p.values<-c(p.1,p.2,p.3,p.4)

FDR.op<-mFDR.alpha.investment(p.values,0.05)

 # to plot, the testing order for parameters (given the standard sequence thr, ths,alp,nscal)
 # is 2,1,3,4

SWRC.mFDR.plot(FDR.op,c(2,1,3,4))

#
#  copy the selected model to selected.model
#

selected.model<-op.2

#
# Plot the data and fitted models
#

plot.wrc.data(data.df,selected.model$Coefficient.estimates,groups="WRB",hvar="h",
thetavar="theta",main="Kenyan soils")

#
# Plot the corresponding pore volume density functions
#

op.nitosol<-plot.pvd(selected.model$Coefficient.estimates[,1],Plot=F)
op.solonetz<-plot.pvd(selected.model$Coefficient.estimates[,2])
lines(op.nitosol,col="red")

#
# Compute and print out the soil quality indices
#

index.values<-SQ.indices(selected.model$Coefficient.estimates)

print.indices(index.values,ret=F)

#
#  Generate a parametric bootstrap set of model parameters
#

i.seed<-char2seed("KenyaBootstrap")

bs.set<-VanGenuchten.fit.compare.bs(data.df,init.vals,g.name="Profile_ID",par.var.full=c(1,1,0,0),
cov.name="WRB",i.seed=i.seed,10,"parametric")
#
print(bs.set)
#

# Compute soil quality index values for each boostrap sample
#
SQ.BSt.indices(bs.set)
#

#################################################################################################
#
#  Run the Green Ampt model for a rainfall data set in Liempe_29_Jan_22_10am_24hrs.dat
#  This contains successive 15-min rainfall depths (cm) over a 24-hour period.

h<-as.matrix(read.table("Liempe_29_Jan_22_10am_24hrs.dat",header=T),ncol=2) #rainfall data
nl<-1  #Number of layers
Dz<-30 #Layer thickness
dmax<-30 #Maximum profile depth


vG<-matrix(nrow=nl,ncol=4) #Matrix for SWRC parameters

#use van Genuchten parameters for first group (Nitisols) 
# or second (Solonetz)

soil.group<-"Solonetz"
c.col<-which(colnames(selected.model$Coefficient.estimates)==soil.group)

vG[1,]<-selected.model$Coefficient.estimates[,c.col]

# init: init rule
#	init = 1 : residual water content
#	init = 2 : field capacity
#	init = 3 : 40% soil moisture deficit

init=3

theta_i<-initialize.wc(vG,init,nl) 

OP<-GALS(nl,vG,Dz,theta_i,dmax,h)

GALS.summary(OP)

GALS.plot(OP,30)

#
#
##################################################################################################

