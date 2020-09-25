import numpy as np
import math
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
xlim=n/5
ylim=n/5
zlim=n/5
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
#psimx,psimy = np.meshgrid(psiix,psiiy)
delt=0.05
time=200
fig = plt.figure(figsize=[16,9])
plt.ion()
cc=40
for i in range(0,time):
    a1=av(t+i*delt)
    a2=av(t+(i+1)*delt)
    D1=D1p(a1)
    dD1=((D1p(a2)-D1p(a1))/delt)
    dpsix=dD1*psiix
    psix=D1*psiix
    dpsiy=dD1*psiiy
    psiy=D1*psiiy
    plt.ion()
    plt.suptitle("Gaussian density contrast")
    #plots
    #Subplot 1
    ax = fig.add_subplot(3,2,1, projection='3d',aspect = 'auto')
    plt.title("Frame No.: %d"%i+" out of %d"%time)
    ax.set_xlim3d(0,xlim)
    ax.set_ylim3d(0,ylim)
    ax.set_zlim3d(-zlim/2,zlim/2)
    ax.plot_wireframe(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),dpsix,rcount=cc,ccount=cc)
    ax.set_zlabel("Vx")
    #Subplot 2
    ax = fig.add_subplot(3,2,2, projection='3d',aspect = 'auto')
    ax.set_xlim3d(0,xlim)
    ax.set_ylim3d(0,ylim)
    ax.set_zlim3d(-zlim/2,zlim/2)
    ax.plot_wireframe(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),dpsiy,rcount=cc,ccount=cc)
    ax.set_zlabel("Vy")
    #Subplot 3
    ax = fig.add_subplot(3,2,3, projection='3d',aspect = 'auto')
    ax.set_xlim3d(0,xlim)
    ax.set_ylim3d(0,ylim)
    ax.set_zlim3d(-zlim/2,zlim/2)
    ax.contourf(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),dpsix)
    ax.set_zlabel("Vx")
    #Subplot 4
    ax = fig.add_subplot(3,2,4, projection='3d',aspect = 'auto')
    ax.set_xlim3d(0,xlim)
    ax.set_ylim3d(0,ylim)
    ax.set_zlim3d(-zlim/2,zlim/2)
    ax.contourf(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),dpsiy)
    ax.set_zlabel("Vy") 
    #Subplot5
    plt.subplot(3,2,5,aspect = 'equal')
    plt.xlim(0,xlim)
    plt.ylim(0,ylim)
    plt.quiver(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),dpsix,dpsiy)
    #Subplot6
    plt.subplot(3,2,6,aspect = 'equal')
    plt.xlim(0,xlim)
    plt.ylim(0,ylim)
    plt.scatter(np.mod(posx+psix,xlim),np.mod(posy+psiy,ylim),s=0.8)
    fig.tight_layout()
    plt.savefig('Images/%d'%i)
    plt.show()
    plt.pause(0.01)
    if i!=time-1:
        plt.clf()
