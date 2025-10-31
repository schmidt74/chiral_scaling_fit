import numpy as np
import scipy.optimize

beta=0.380
delta=4.824
Delta=beta*delta

h3=0.306
h5=-0.00338
theta0=1.359

def h(theta):
    return (theta+h3*np.power(theta,3)+h5*np.power(theta,5))*np.power(1-np.power(theta,2)/np.power(theta0,2),2)

def z(theta):
    return (1-np.power(theta,2))/(np.power(theta0,2)-1)*np.power(theta0,1/beta)*np.power(h(theta)/h(1),-1/Delta)

def theta_fun(x):
    root=scipy.optimize.bisect(lambda y: x-z(y), 1e-5,0.999*theta0, args=(), xtol=2e-12, rtol=8.881784197001252e-16, maxiter=100, full_output=False, disp=True)
    return root
    
def fG(theta):
    return theta*np.power(h(theta)/h(1),-1/delta)

      
if (__name__ == "__main__"):
    import matplotlib.pyplot as plt

    #generate equidistant points in theta
    THETA=np.linspace(0,0.99*theta0,200)
    #calculate z(theta)
    x=z(THETA)
    #calculate fG(theta)
    y=fG(THETA)

    #generate equidistant points in z
    X=np.linspace(-30,30,100)
    #calculate theta(z)
    THETA2=np.array([])
    for i in range(len(X)):
        #print(i, X[i], t(X[i]), fG(theta_fun(X[i])))
        THETA2=np.append(THETA2,theta_fun(X[i]))
    #calculte fG(theta)
    Y=fG(THETA2)

    plt.plot(THETA,x,"x")
    plt.plot(THETA2,X,"+")
    plt.show()

    plt.plot(x,y,"x")
    plt.plot(X,Y,"+")
    plt.xlim(-30,30)
    plt.ylim(0,4)
    plt.show()

