import math, cmath
import matplotlib
from matplotlib import pyplot as plt
#initialization
T = int(input("Enter the Time Period (even): "))
N = int(input("Enter the number equally spaced samples (even): "))
xr = []
xi = []
fr = []
fi = []
w = []
delt= T/N
const=-2*math.pi*(1j)
wn=cmath.exp(const/N)
f=[]
#range of t: [0,T)
#range of n: [0,N-1]
pi = math.pi
def complex_dft(xr, xi, n):
	rex = [0] * n
	imx = [0] * n
	for k in range(0, n):  # exclude n
		rex[k] = 0
		imx[k] = 0
	for k in range(0, n):  # for each value in freq domain
		for i in range(0, n):  # correlate with the complex sinusoid
			sr =  math.cos(2 * pi * k * i / n)
			si = -math.sin(2 * pi * k * i / n)
			rex[k] += xr[i] * sr - xi[i] * si
			imx[k] += xr[i] * si + xi[i] * sr
	return rex, imx

# FFT version based on the original BASIC program
def fft_basic(rex, imx, n):
	pi = 3.141592653589793
	m = int(math.log(n, 2))  # float to int
	j = n / 2

	# bit reversal sorting
	for i in range(1, n - 1):  # [1,n-2]
		if i >= j:
			# swap i with j
			rex[i], rex[j] = rex[j], rex[i]
			imx[i], imx[j] = imx[j], imx[i]
		k = n / 2
		while (1):
			if k > j:
				break
			j -= k
			k /= 2
		j += k

	for l in range(1, m + 1):  # each stage
		le = int(math.pow(2, l))  # 2^l
		le2 = le / 2
		ur = 1
		ui = 0
		sr =  math.cos(pi / le2)
		si = -math.sin(pi / le2)
		for j in range(1, le2 + 1):  # [1, le2] sub DFT
			for i in xrange(j - 1, n - 1, le):  #  for butterfly
				ip = i + le2
				tr = rex[ip] * ur - imx[ip] * ui
				ti = rex[ip] * ui + imx[ip] * ur
				rex[ip] = rex[i] - tr
				imx[ip] = imx[i] - ti
				rex[i] += tr
				imx[i] += ti
			tr = ur
			ur = tr * sr - ui * si
			ui = tr * si + ui * sr
def f1(t):
    if t!=0:
        f1 = (math.sin(2*math.pi*t)/t)
    else:
        f1=1
    return f1
def f2(t):
    f2 = (math.exp(-0.5*t*t)/math.sqrt(2*math.pi))
    return f2
for n in range(0,N):
    f.append(f1(n*delt))
    fr.append(f[n].real)
    fi.append(f[n].imag)
    p = cmath.sin(2 * pi * (1j) / n)
    xr.append(p.real)
    xi.append(p.imag)
rex, imx = complex_dft(xr, xi, n)
fft_basic(fr, fi, n)
'''plt.plot(w,freal)
plt.title("Real part Vs w")
plt.show()
plt.plot(w, fimag)
plt.title("Imaginary part Vs w")
plt.show()'''

