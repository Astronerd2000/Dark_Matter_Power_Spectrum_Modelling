import math, cmath,timeit
import matplotlib
from matplotlib import pyplot as plt
#initialization
T = int(input("Enter the Time Period (even): "))
N = int(input("Enter the number equally spaced samples (even): "))
freal = []
fimag = []
w = []
delt= T/N
const=-2*math.pi*(1j)
#range of t: [-T/2,T/2)
#range of n: [-N/2, N/2-1]
def DFT(m,N,delt):
    fm = 0
    for n in range(-N//2,N//2):
        f=f2(n*delt)
        fm+=(f*cmath.exp((const*m*n)/N))
    return fm
def f1(t):
    if t!=0:
        f1 = (math.sin(2*math.pi*t)/t)
    else:
        f1=1
    return f1
def f2(t):
    f2 = (math.exp(-0.5*t*t)/math.sqrt(2*math.pi))
    return f2
start = timeit.timeit()
for m in range(-N//2,N//2):
    freal.append((DFT(m,N,delt).real*delt))
    fimag.append((DFT(m,N,delt).imag*delt))
    w.append(2*math.pi*m/T)
    #print(DFT(m,N,delt))
end = timeit.timeit()
print(end - start)
plt.subplot(1,2,1)
plt.plot(w,freal)
plt.title("Real part Vs w(=2$\pi$m/T)")
plt.subplot(1,2,2)
plt.plot(w, fimag)
plt.title("Imaginary part Vs w(=2$\pi$m/T)")
plt.suptitle("Gaussian Function")
plt.show()
