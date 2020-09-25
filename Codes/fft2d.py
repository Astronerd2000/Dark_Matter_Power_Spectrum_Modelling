import math, cmath,timeit
import matplotlib, scipy
from matplotlib import pyplot as plt
import numpy as np
from scipy.fftpack import fft, ifft,fftfreq,rfftfreq, fftshift
#initialization
T = 40#int(input("Enter the Time Period (even): "))
N = 1600#int(input("Enter the number equally spaced samples (even): "))
delt= T/N
const=-2*math.pi*(1j)
t = np.linspace(-T/2, T/2, N)
xf = rfftfreq(N,1/T)[1:]
f1 = np.piecewise(t, [t!=0 , t==0], [np.sin(2*np.pi*t)/t,1])
f2 = (np.exp(-0.5*t*t)/np.sqrt(2*math.pi))
start = timeit.timeit()
y= fftshift(fft(f2))[1:]
iy=(ifft(y))
yf = np.abs(y)
iyf = np.abs(iy)
end = timeit.timeit()
print(xf)
print(end - start)
plt.plot(xf,iyf)
plt.grid()
plt.show()

