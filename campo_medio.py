import numpy as np
import matplotlib.pyplot as plt

def Bj(x,J):
    B_J = (2*J + 1)/(2*J)*(1/np.tanh((2*J+1)/(2*J)*x)) - 1/(2*J)*(1/np.tanh(1/(2*J)*x))
    return B_J

def inversa(f,y,a,b,tol=1e-6,iter=10000):

    
    a = np.full(y.shape,a)
    b = np.full(y.shape,b)

    for n in range(iter):
        
        c = (a+b)/2

        esquerda = f(a) - y
        meio = f(c) - y

        cond = esquerda*meio <= 0

        a = np.where(cond,a,c)
        b = np.where(cond,c,b)

        if np.all(b - a < tol):
            break

    return (a+b)/2

        
y = np.linspace(0,1.5,100)
B = inversa(lambda x: Bj(x,7/2), y, -1.1,10)

plt.plot(y,B)
plt.grid(True)
plt.show()