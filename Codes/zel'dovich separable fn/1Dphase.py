import numpy as np
import math
import matplotlib,scipy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.fftpack import fft, ifft,fftfreq,rfftfreq, fftshift
print('''Different eras:
      1. Radiation Dominated era
      2. Matter Dominated era
      3. Dark-energy Dominated era''')
ch=int(input("Enter the era:")) #set to 3
H=0.7088 #km/s/Mpc
t=-3
def av(t):
    if ch==1:
        ax = t**0.5
    elif ch==2:
        ax = t**(2/3)
    elif ch==3:
        ax = math.exp(H*t)
    return ax
n=int(input("Enter the number of particles:")) #set to 50
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
#Zeldovich
A=1
B=1
xlim=n/5
init = np.linspace(0,xlim,n)
#Fourier Transform
T = init[-1]-init[0]
N = n
x=init
delta = A*np.sin(B*x)
xf = (fftfreq(N,T/N))*2*np.pi
xf+=1e-16
deltilde=(fft(delta))
psitildei=(((1j)*deltilde)/(xf))
psii=np.real((ifft(psitildei)))
delt=0.03
time=200
for i in range(0,time):
    a1=av(t+i*delt)
    a2=av(t+(i+1)*delt)
    D1=D1p(a1)
    dD1=((D1p(a2)-D1p(a1))/delt)
    psi=D1*psii
    dpsi=dD1*psii
    plt.xlim(0,xlim)
    plt.ylim(-n/20,n/20)
    plt.ion()
    plt.suptitle("Sinusoidal density contrast")
    plt.title("Frame No.: %d"%i+" out of %d"%time)
    plt.scatter(np.mod(init+psi,xlim),dpsi,s=10)
    #plt.savefig('Images/%d'%i)
    plt.show()
    plt.pause(0.001)
    if i!=time-1:
        plt.clf()
