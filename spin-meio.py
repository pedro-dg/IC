import numpy as np
import matplotlib.pyplot as plt

B = np.linspace(-10,10,1000)
T = np.linspace(0.1,100,1000)
m = np.array([-1/2,1/2])
expoente = m*(B[:,np.newaxis]/T)[...,np.newaxis]
Z = np.sum(np.exp(expoente),axis=2)

F = -T[np.newaxis,:]*np.log(Z)
dB = B[1] - B[0]
M = -np.gradient(F,dB,axis=0)

plt.plot(T,M[-1,:],color="red",label="B=10T")
plt.plot(T,M[500,:],color="blue",label="B=5T")
plt.plot(T,M[0,:],color="green",label="B=0T")
plt.grid(True)
plt.xlim(0,100)
plt.xlabel(r'$k_B$T/$\mu_B$')
plt.ylabel("Magnetização (A/m)")
plt.title("Magnetização x Temperatura")
plt.legend()
plt.show()

plt.plot(B,M[:,0], color="red",label="T=0.1K")
plt.plot(B,M[:,20], color="blue",label=f"T={T[20]:.2f}K")
plt.plot(B,M[:,40], color="green",label=f"T={T[40]:.2f}K")
plt.grid(True)
plt.xlim(-10,10)
plt.xlabel(r'$\mu_B$B/$k_B$')
plt.ylabel("Magnetização (A/m)")
plt.title("Magnetização x Campo Magnético")
plt.legend()
plt.show()