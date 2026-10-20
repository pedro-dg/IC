import numpy as np
import matplotlib.pyplot as plt

#Definindo parâmetros
B = np.linspace(0.1,10,1000)
T = np.linspace(0.1,10,1000)
x = B[:,np.newaxis]/T
J = 7/2
m = np.arange(J,J+1)
#Cálculo de Z
expoente = -m*x[...,np.newaxis]
z = np.sum(np.exp(expoente),axis=2)

#Cálculo da Magnetização
F = -T[:, np.newaxis]*np.log(z)
dB = B[1] - B[0]
M = -np.gradient(F,dB,axis=0)

#Customização
plt.grid(True)
plt.xlabel("T(K)")
plt.ylabel("M/M_0")
plt.legend()
plt.title("Magnetização x Temperatura")

#Plot
plt.plot(T,M[-1,:],color="blue",label="B=1000T")
plt.plot(T,M[0,:],color="red",label="B=0.1T")
plt.plot(T,M[500,:],color="green",label="B=500T")



plt.show()

