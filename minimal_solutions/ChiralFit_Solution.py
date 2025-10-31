# Chiral Fit 

import os
import numpy as np
import ChiralData_Solution as dat
import ChiralFunc_Solution as fun
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

# stack (T, m) for practical purpos
xdata=np.vstack((T,m))
   
################
# Fit Function #
################
def residual(params, xdata, M, err):
  Tc=params['Tc']
  z0=params['z0']
  h0=params['h0']
  T,m=xdata
  t=np.divide(np.add(T,-Tc),Tc)
  h=np.abs(np.divide(m,h0))
  z=z0*np.divide(t,np.power(m,1.0/fun.Delta))
  model=np.zeros(len(t))
  # do not allow for large z values
  # which might give numerical problems with inverting z(theta)
  for i in range(len(t)):
    if(np.abs(z[i])>20):
      z[i]=20.0*np.sign(z[i])
    model[i]=fun.fG(fun.theta_fun(z[i]))*np.power(h[i],1.0/fun.delta)+b1*h[i]*t[i]
    
  return (M-model)/err

def _M(params, xdata):
  Tc=params['Tc']
  z0=params['z0']
  h0=params['h0']
  T,m=xdata
  t=np.divide(np.add(T,-Tc),Tc)
  h=np.abs(np.divide(m,h0))
  z=z0*np.divide(t,np.power(m,1.0/fun.Delta))
  model=np.zeros(len(t))
  for i in range(len(t)):
    if(np.abs(z[i])>20):
      z[i]=20.0*np.sign(z[i])
    model[i]=fun.fG(fun.theta_fun(z[i]))*np.power(h[i],1.0/fun.delta)+b1*h[i]*t[i]
    
  return model

###############
# Perform Fit #
###############
params=Parameters()
params.add('Tc', value=144)
params.add('z0', value=1.5)
params.add('h0', value=1.0e-6)
b1=0.0 ## you might include a regular term here

out=minimize(residual, params, args=(xdata, M, err))
print(fit_report(out))

##########
# Plot 1 #
##########

#the chiral condensat as function of T for different masses
colorid=0
for dkey in ["l568m160", "l568m80", "l408m40", "l328m27", "l328m20"]:
   index=x.GetIndexFromDkey(dkey)
   Tdata=x.data[index]["sdata"]["T"]
   Mdata=x.data[index]["sdata"]["M"]
   edata=x.data[index]["sdata"]["M_err"]
   colorstr="C{:d}".format(colorid)
   plt.errorbar(Tdata,Mdata,yerr=edata, fmt='+', color=colorstr, ecolor=colorstr, capsize=3)
   Tfit=np.linspace(Tcenter-Tdelta,Tcenter+Tdelta,100)
   mfit=np.array([x.data[index]["ml_by_ms"]]*len(Tfit))
   mfit=np.divide(1.0,mfit)
   xfit=np.vstack((Tfit,mfit))
   Mfit=_M(out.params,xfit)
   plt.plot(Tfit, Mfit, colorstr)
   colorid+=1
plt.show()

#########
# Plot2 #
#########

#the chiral condensat as function of T for different masses
colorid=0
for dkey in ["l568m160", "l568m80"]:
   index=x.GetIndexFromDkey(dkey)
   Tdata=x.data[index]["sdata"]["T"]
   Mdata=x.data[index]["sdata"]["M"]
   edata=x.data[index]["sdata"]["M_err"]
   tdata=np.divide(np.add(Tdata,-out.params['Tc']),out.params['Tc'])
   mdata=np.divide(1.0,np.array([x.data[index]["ml_by_ms"]]*len(tdata)))
   zdata=out.params['z0']*np.divide(tdata,np.power(mdata,1/fun.Delta))
   hdata=mdata/out.params['h0']
   Mdata=np.divide(Mdata,np.power(hdata,1/fun.delta))
   edata=np.divide(edata,np.power(hdata,1/fun.delta))
   colorstr="C{:d}".format(colorid)
   plt.errorbar(zdata,Mdata,yerr=edata, fmt='+', color=colorstr, ecolor=colorstr, capsize=3)
   Tfit=np.linspace(Tcenter-Tdelta,Tcenter+Tdelta,100)
   tfit=np.divide(np.add(Tfit,-out.params['Tc']),out.params['Tc'])
   mfit=np.array([x.data[index]["ml_by_ms"]]*len(Tfit))
   mfit=np.divide(1.0,mfit)
   zfit=out.params['z0']*np.divide(tfit,np.power(mfit,1/fun.Delta))
   xfit=np.vstack((Tfit,mfit))
   Mfit=_M(out.params,xfit)
   fGfit=np.array([])
   for i in range(len(tfit)):
      fGfit=np.append(fGfit,fun.fG(fun.theta_fun(zfit[i])))
   plt.plot(zfit, fGfit, colorstr)
   colorid+=1
plt.show()
