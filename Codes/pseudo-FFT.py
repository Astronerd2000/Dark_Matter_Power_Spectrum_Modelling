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
const=(-2*math.pi*(1j))/N
#range of t: [-T/2,T/2)
#range of n: [-N/2,N/2-1]
def pfft(x):
    fm = 0
    for n in range(-N//2,0):
        fe=f1(2*n*delt)
        fo=f1((2*n+1)*delt)
        fm+=(cmath.exp(const*m*2*n)*fe + cmath.exp(const*(m)*(2*n+1))*fo)
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
    freal.append((pFFT(m,N,delt).real*delt))
    fimag.append((pFFT(m,N,delt).imag*delt))
    w.append(2*math.pi*m/T)
    #print(pFFT(m,N,delt))
end = timeit.timeit()
print(end - start)
plt.plot(w,freal)
plt.title("Real part Vs w")
plt.show()
plt.plot(w, fimag)
plt.title("Imaginary part Vs w")
plt.show()

