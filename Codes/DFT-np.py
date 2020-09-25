import math, cmath,timeit
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
#initialization
T = int(input("Enter the Time Period (even): "))
N = int(input("Enter the number equally spaced samples (even): "))
freal = []
fimag = []
w = []
delt= T/N
const=(-2*math.pi*(1j))/N
wl = np.empty([N,N],dtype =complex)
fn = np.empty([N],dtype =complex)
#range of t: [-T/2,T/2)
#range of n: [-N/2, N/2-1]
def DFT(N,delt):
    for n in range(0,N):
        fn[n] = f1((n-N//2)*delt)
        for m in range(0,N):
            wl[m][n]= cmath.exp(const*(m-N//2)*(n-N//2))
    return np.dot(wl,fn)
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
F=DFT(N,delt)
for m in range(0,N):
    freal.append((F[m].real*delt))
    fimag.append((F[m].imag*delt))
    w.append(2*math.pi*(m-N//2)/T)
    #print(DFT(N,delt)[m])
end = timeit.timeit()
print(end - start)

plt.plot(w,freal)
plt.title("Real part Vs w")
plt.show()
plt.plot(w, fimag)
plt.title("Imaginary part Vs w")
plt.show()
