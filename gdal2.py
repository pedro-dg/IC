import numpy as np
import matplotlib.pyplot as plt

#Definindo os parâmetros

J = 7/2
T = np.linspace(0.1,10,1000)
B = np.linspace(0,10,1000)
x = B[:,np.newaxis]/T
m = np.arange(-J,J+1)

#Efetuando o cálculo da função partição

def calculo_z(m,x,T):
    expoente = m*x[...,np.newaxis]
    z = np.sum(np.exp(expoente),axis=2)
    F = -T[np.newaxis,:]*np.log(z)
    return z,F

#Efetuando o cálculo de S
def calculo_s(T,F):
    dT = T[1] - T[0]
    s = -np.gradient(F,dT,axis=1)
    return s

#Efetuando o cálculo de M
def calculo_M(F):
    dB = B[1] - B[0]
    M = -np.gradient(F,dB,axis=0)
    return M

def tratar_dados(x,s):
    x_flat = x.ravel()
    s_flat = s.ravel()
    sort_index = np.argsort(x_flat)
    x_sorted = x_flat[sort_index]
    s_sorted = s_flat[sort_index]
    return x_sorted,s_sorted

def plot_s(x_sorted,s_sorted):
    plt.grid(True)
    plt.xlabel("x = B/T")
    plt.ylabel("S (J/mol.K)")
    plt.plot(x_sorted,s_sorted)
    plt.show()
    plt.close()
    return
def plot_M(T,M):
    plt.plot(T,M[-1,:],color="blue",label="B = 10T")
    plt.plot(T,M[0,:],color="green",label="B = 0T")
    plt.plot(T,M[500,:],color="red",label=f"B = {B[500]:.2f}T")
    plt.legend()
    plt.xlabel("T(K)")
    plt.ylabel("M/M_0")
    plt.xlim(0.1,10)
    plt.grid(True)
    plt.title("Magnetização x Temperatura a campo constante")
    plt.show()
#Plot de S a B constante
def plot_sB(s,T):
    plt.grid(True)
    plt.xlabel(r"$k_B$T/$\mu_B$")
    plt.ylabel("S (J/mol.K)")
    plt.plot(T,s[0,:],color="red",label="B=0T")
    plt.plot(T,s[500,:],color="green",label="B=5T")
    plt.plot(T,s[-1,:],color="blue",label="B=10T")
    plt.legend()
    plt.xlim(0,10)
    plt.title("Entropia x Temperatura a campo constante")
    plt.show()
    plt.close()
    return


#Plot da função analítica
def plot_analitico():
    x_an = np.linspace(0,10,1000)
    s_an = np.log(np.sinh(8*x_an)/np.sinh(x_an)) -8*x_an*(1/np.tanh(8*x_an)) + x_an*(1/np.tanh(x_an))
    plt.plot(x_an,s_an)
    plt.grid(True)
    plt.show()
    return

z,F = calculo_z(m,x,T)
s = calculo_s(T,F)
x_s,s_s = tratar_dados(x,s)
M = calculo_M(F)


plot_s(x_s,s_s)
plot_sB(s,T)
plot_M(T,M)
plot_analitico()
