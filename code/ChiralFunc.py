import numpy as np
import scipy.optimize

beta=0.380
delta=4.824
Delta=beta*delta

h3=0.306
h5=-0.00338
theta0=1.359

# helper function, used inside the Schofield parametrization
def h(theta):
    # add code here
    return result

def z(theta):
    # add code here
    return result

# inverse of the above.
# use numerical method to solve z(theta)=z0
# e.g. besection from scipy.optimie
def t(x):
    # add code here
    return root

#scaling function in the Schofield parametriziation
def fG(theta):
    # add code here
    return result

      
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
        #print(i, X[i], t(X[i]), fG(t(X[i])))
        THETA2=np.append(THETA2,t(X[i]))
    #calculte fG(theta)
    Y=fG(THETA2)

    ####################
    # Plot 1: z(theta) #
    ####################
    plt.plot(THETA,x,"x")
    plt.plot(THETA2,X,"+")
    plt.show()

    #################
    # Plot 2: fG(z) #
    #################
    plt.plot(x,y,"x")
    plt.plot(X,Y,"+")
    plt.xlim(-30,30)
    plt.ylim(0,4)
    plt.show()

