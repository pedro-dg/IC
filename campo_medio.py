import numpy as np
import matplotlib.pyplot as plt

#Definindo a função de Brillouin
def Bj(x,J):
    B_J = (2*J + 1)/(2*J)*(1/np.tanh((2*J+1)/(2*J)*x)) - 1/(2*J)*(1/np.tanh(1/(2*J)*x))
    return B_J

C = np.linspace(-10,-1,50)
D = np.linspace(1,10,50)
x = np.linspace(0.1,10,100)


#Encontrando as raízes de uma função qualquer

def inversa(f,a,b,tol=1e-6,iter=1000):

    if f(a)*f(b) > 0:
        raise ValueError("O intervalo não possui mudança de sinal!")
    if f(a) == 0:
        return a
    if f(b) == 0:
        return b
    
    for n in range(iter):

        c = (a+b)/2

        if (b-a)/2 < tol:
            return c

        if f(a)*f(c) < 0:
            b = c
        else:
            a = c

    raise RuntimeError("O intervalo não convergiu!")

#Encontrando as raízes para C fixo

def MxT(C,D):
    x0 = np.array([])
    for d in D:
        intercepto = inversa(lambda x: Bj(x,1/2) - (C + d*x),0.1,10)
        x0 = np.append(x0,intercepto)
    return x0

#Encontrando as raízes para D fixo

def MxH(C,D):
    x0 = np.array([])
    for c in C:
        intercepto = inversa(lambda x: Bj(x,1/2) - (c + D*x),0.1,10)
        x0 = np.append(x0,intercepto)
    return x0

x0_T = MxT(C[-1],D)
y = Bj(x0_T,1/2)
plt.plot(D,y)
plt.xlabel('kT/2zJS²')
plt.ylabel("M")
plt.grid(True)
plt.show()

x0_H = MxH(C,D[-1])
y = Bj(x0_H,1/2)
plt.plot(C,y)
plt.xlabel(r'$g\mu_0 H_0/2zS$')
plt.ylabel('M')
plt.grid(True)
plt.show()
