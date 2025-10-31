# Chiral data module

import os
import numpy as np
import scipy.interpolate

class ChiralData():

   def __init__(self,mydir):
      self.datadir=mydir
      self.data=[]
      self.Nt=[6,8,12]
      self.V =[16,24,32,40,42,48,56,60,72]
      self.ml=[20,27,40,60,80,160]

      id=0
      for i in self.Nt:
         for j in self.V:
            for k in self.ml: 
               dkey="l"+str(j)+str(i)+"m"+str(k)
               file=self.datadir+"/"+self.GetFileName(dkey)
               if (os.path.isfile(file)):
                  id=id+1
                  print("found data set {:2d}:  Nt = {:2d}  V = {:2d}  ml/ms = {:3d}  key = {:s}".format(id,i,j,k,dkey))
                  dset={"dkey":dkey,"Nt":i,"V":j,"ml_by_ms":k,"sdata":{}}
                  self.data.append(dset)

   def ReadData(self):
        for index in range(len(self.data)):
           self.ReadDataSet(index)

   def ReadDataSet(self,index):
      dkey=self.data[index]["dkey"]
      file=self.datadir+"/"+self.GetFileName(dkey)
      print(file)
      (beta,T,mass, M, M_err, chi, chi_err)=np.loadtxt(file,usecols=(0,1,2,3,4,5,6), unpack=True)
      self.data[index]["sdata"]={"beta":beta, "T":T, "M":M, "M_err":M_err, "chi":chi, "chi_err":chi_err}
      
   def GetFileName(self,dkey):
        return "M_chiM_"+dkey+".dat"
     
   def GetIndexFromDkey(self,dkey):
      for index in range(len(self.data)):
         if (self.data[index]["dkey"]==dkey):
            return index

      
if (__name__ == "__main__"):
    import matplotlib.pyplot as plt

    #############
    # Read Data #
    #############
    x=ChiralData("../data/")
    x.ReadData()


    #############################
    # Plot 1: Chrial Condensate #
    #############################

    TInter=np.linspace(138,168,100)
    Tcut=150
    colorid=0
    for dkey in ["l568m160", "l568m80", "l408m40", "l328m27", "l328m20"]:
       colorstr="C{:d}".format(colorid)
       index=x.GetIndexFromDkey(dkey)
       T=x.data[index]["sdata"]["T"]
       M=x.data[index]["sdata"]["M"]
       plt.plot(T,M,colorstr+"o")
       MInter=scipy.interpolate.interp1d(T,M,kind="cubic",fill_value="extrapolate")
       plt.plot(TInter,MInter(TInter),colorstr+"-")
       colorid=colorid+1
    plt.show()

    #################################
    # Plot 2: Chrial Susceptibility #
    #################################

    # add code here!
