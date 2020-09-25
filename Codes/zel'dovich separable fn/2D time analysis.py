import numpy as np
import math, time
import matplotlib,scipy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.fftpack import fft, ifft,fft2, ifft2, fftfreq,rfftfreq, fftshift
print('''Different eras:
      1. Radiation Dominated era
      2. Matter Dominated era
      3. Dark-energy Dominated era''')
ch=int(input("Enter the era:")) #set to 3
H=0.7088 #km/s/Mpc
t=-10
def av(t):
    if ch==1:
        ax = t**0.5
    elif ch==2:
        ax = t**(2/3)
    elif ch==3:
        ax = math.exp(H*t)
    return ax
Om=float(input("Enter the value of Omega_m:"))  #set to 1
Ol=float(input("Enter the value of Omega_lambda:")) #set to 0
def D1p(a):
    if Om==1 and Ol==0:
        D1p = a
    elif Om<1 and Ol==0:
        x=((1/Om)-1)
        D1p = (1+(3/x)+3*math.pow((1+x)/x**3,0.5)*math.log(((1+x)**0.5)-(x**0.5)))
    else:
        D1p = ((2.5*a*Om)/((Om**(4/7))-Ol+(1+0.5*Om)*(1+(Ol/70))))
    return D1p
#Zeldovich 2D time
N=[]
TIME=[]
nvalues=[]
for i in range(1,11):
    nvalues.append(math.pow(2,i))
print(nvalues)
for n in nvalues:
    start = time.time()
    for i in range(100):
        xlim=n/5
        ylim=n/5
        A=1
        B=1
        C=1
        initx = np.linspace(0,xlim,n)
        inity = np.linspace(0,ylim,n)
        posx,posy = np.meshgrid(initx,inity)
        #Fourier Transform
        x, y = np.meshgrid(initx, inity, sparse=False, indexing='ij')
        #FUNCTION
        sigma=0.05
        #delta = 1/(2*np.pi*sigma**2) * np.exp(-0.5 * (x ** 2 + y ** 2)/sigma**2)
        delta = np.sin(x)*np.sin(y)
        F = fft2(delta)
        T = initx[-1]-initx[0]
        xaf = (fftfreq(F.shape[0],d=2))*2*np.pi
        xaf+=1e-16
        yaf = (fftfreq(F.shape[1],d=2))*2*np.pi
        yaf+=1e-16
        kx, ky = np.meshgrid(xaf, yaf, sparse=False, indexing='ij')
        #x-coordinate
        psitildeix=(((1j)*F*kx)/(kx**2+ky**2))
        psiix=np.real((ifft2(psitildeix)))
        #y-coordinate
        psitildeiy=(((1j)*F*ky)/(kx**2+ky**2))
        psiiy=np.real((ifft2(psitildeiy)))
        i+=1
    end = time.time()
    N.append(n*n)
    print(n*n)
    TIME.append((end-start)/100)
plt.plot(N,TIME)
plt.xlabel("Number of particles")
plt.ylabel("Time taken")
plt.title("2D time analysis")
plt.show()
