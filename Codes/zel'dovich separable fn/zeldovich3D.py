import numpy as np
import math
import matplotlib,scipy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.fftpack import fft, ifft,fftn, ifftn,fftfreq,rfftfreq, fftshift
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
n=int(input("Enter the number of particles:")) #set to 20
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
xlim=n/5
ylim=n/5
zlim=n/5
A=1
B=1
C=1
initx = np.linspace(0,xlim,n)
inity = np.linspace(0,ylim,n)
initz = np.linspace(0,zlim,n)
posx,posy,posz = np.meshgrid(initx,inity,initz)
#Fourier Transform
x, y, z= np.meshgrid(initx, inity, initz, sparse=False, indexing='ij')
#FUNCTION
sigma=0.05
delta = 1/(2*np.pi*sigma**2) * np.exp(-0.5 * (x ** 2 + y ** 2 + z ** 2)/sigma**2)
#delta = np.sin(x)*np.sin(y)*np.sin(z)
F = fftn(delta)
T = initx[-1]-initx[0]
xaf = (fftfreq(F.shape[0],d=3))*2*np.pi
xaf+=1e-16
yaf = (fftfreq(F.shape[1],d=3))*2*np.pi
yaf+=1e-16
zaf = (fftfreq(F.shape[2],d=3))*2*np.pi
zaf+=1e-16
kx, ky, kz= np.meshgrid(xaf, yaf, zaf, sparse=False, indexing='ij')
#x-coordinate
psitildeix=(((1j)*F*kx)/(kx**2+ky**2+kz**2))
psiix=np.real((ifftn(psitildeix)))
#y-coordinate
psitildeiy=(((1j)*F*ky)/(kx**2+ky**2+kz**2))
psiiy=np.real((ifftn(psitildeiy)))
#z-coordinate
psitildeiz=(((1j)*F*kz)/(kx**2+ky**2+kz**2))
psiiz=np.real((ifftn(psitildeiz)))
#psimx,psimy = np.meshgrid(psiix,psiiy)
delt=0.05
time=250
fig = plt.figure()
plt.ion()
for i in range(0,time):
    ax = fig.add_subplot(111, projection='3d')
    a=av(t+i*delt)
    D1=D1p(a)
    psix=D1*psiix
    psiy=D1*psiiy
    psiz=D1*psiiz
    plt.suptitle("Gaussian density contrast")
    plt.title("Frame No.: %d"%i+" out of %d"%time)
    ax.set_xlim3d(0,xlim)
    ax.set_ylim3d(0,ylim)
    ax.set_zlim3d(0,zlim)
    ax.scatter3D(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),np.mod(posz+psiz,zlim),s=1)
    plt.savefig('Images/%d'%i)
    plt.show()
    plt.pause(0.01)
    if i!=time-1:
        plt.clf()
