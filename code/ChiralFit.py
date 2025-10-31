# Chiral Fit 

import os
import numpy as np
import ChiralData as dat
import ChiralFunc as fun
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit


#############
# Read Data #
#############

x=dat.ChiralData("../data/")
x.ReadData()

Tcenter=144
Tdelta=14

####################
# Prepare Data:    #
# 1. (T,m) = xdata #
# 2. M             #
# 3. err           #
####################
T=np.array([])
m=np.array([])
M=np.array([])
err=np.array([])
# choose here the data sets to be included
# for dkey in ["l568m160", "l568m80", "l408m40", "l328m27", "l328m20"]:
for dkey in ["l568m160", "l568m80"]:
   index=x.GetIndexFromDkey(dkey)
   # include only temperatures in the range [Tcenter-Tdelta,Tcenter+Tdelta]
   for i in range(len(x.data[index]["sdata"]["T"])):
      if np.abs(x.data[index]["sdata"]["T"][i]-Tcenter)<Tdelta :
         m=np.append(m,1.0/x.data[index]["ml_by_ms"])
         T=np.append(T,x.data[index]["sdata"]["T"][i])
         M=np.append(M,x.data[index]["sdata"]["M"][i])
         err=np.append(err,x.data[index]["sdata"]["M_err"][i])

# stack data for practical purpos 
xdata=np.vstack((T,m))
   
################
# Fit Function #
################
def residual(params, xdata, M, err):
  Tc=params['Tc']
  z0=params['z0']
  h0=params['h0']
  T,m=xdata
  # calculate t
  # calculate h
  # calculate z
  # calculate model prediciton 
    
  return (M-model)/err

def _M(params, xdata):
  Tc=params['Tc']
  z0=params['z0']
  h0=params['h0']
  T,m=xdata
  # calculate t
  # calculate h
  # calculate z
  # calculate model prediciton 
    
  return model

###############
# perform fit #
###############
params=Parameters()
params.add('Tc', value=144)
params.add('z0', value=1.5)
params.add('h0', value=1.0e-6)
a_t=0.0 ## you might include a regular term here

out=minimize(residual, params, args=(xdata, M, err))
print(fit_report(out))

##########
# plot 1 #
##########
colorid=0
for dkey in ["l568m160", "l568m80", "l408m40", "l328m27", "l328m20"]:
   index=x.GetIndexFromDkey(dkey)
   Tdata=x.data[index]["sdata"]["T"]
   Mdata=x.data[index]["sdata"]["M"]
   edata=x.data[index]["sdata"]["M_err"]
   colorstr="C{:d}".format(colorid)
   # plot data with error bars
   # calculate model predictions
   # plot 
   plt.errorbar(Tdata,Mdata,yerr=edata, fmt='+', color=colorstr, ecolor=colorstr, capsize=3)
   colorid+=1
   Tfit=np.linspace(Tcenter-Tdelta,Tcenter+Tdelta,100)
   mfit=np.array([x.data[index]["ml_by_ms"]]*len(Tfit))
   mfit=np.divide(1.0,mfit)
   xfit=np.vstack((Tfit,mfit))
   Mfit=_M(out.params,xfit)
   plt.plot(Tfit, Mfit, colorstr)
plt.show()

##########
# plot 2 #
##########

# create the plot with the data colapse here
