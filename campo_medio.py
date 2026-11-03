import numpy as np
import matplotlib.pyplot as plt

C = np.linspace(-0.5,0,100)
D = np.linspace(0.1,0.99,100)
x = np.linspace(0,2,1000)


#Definindo a função de Brillouin

def Bj(x,J):

    x = np.asarray(x,dtype=float)
    B_J = np.zeros_like(x)

    mask = x != 0

    B_J[mask] = (2*J + 1)/(2*J)*(1/np.tanh((2*J+1)/(2*J)*x[mask])) - 1/(2*J)*(1/np.tanh(1/(2*J)*x[mask]))
    return B_J

#Encontrar intervalos onde há mudança de sinal em f
def intervalo(f,xmin,xmax,n=1000):
    
    pontos = np.linspace(xmin,xmax,n)
    
    #Calculando a função f em cada ponto
    
    valores = f(pontos)

    #Índices onde a condição f(a)*f(b) < 0 é verdadeira

    indices = np.where(valores[:-1]*valores[1:] < 0)[0]

    return [(pontos[i], pontos[i+1]) for i in indices]
    


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
    
        f = lambda x: Bj(x,1/2) - (C + d*x)

        a = 1e-6
        b = 1/d + 1

        intercepto = inversa(f,a,b)
        x0 = np.append(x0,intercepto)

    return x0

#Encontrando as raízes para D fixo

def MxH(C,D):

    resultados = []
    for c in C:
        f = lambda x: Bj(x,1/2) - (c + D*x)

        intervalos = intervalo(f,-30,30,10000)

        raizes = [inversa(f,a,b) for a,b in intervalos]

        resultados.append((c,raizes))
        
    return resultados

#Calculando os intervalos onde há raízes para cada C com D constante
def raizes_C(C,D):

    ab = []

    for c in C:

        valores = intervalo(
            lambda x: Bj(x,1/2) - (c + D*x),
            -10,10
            )
        ab.append((c,valores))

    for c,interval in ab:
        for a,b in interval:
            print(f"{c}: \t ({a},{b})")

#Calculando os intervalos onde há raízes para cada D com C constante

def raizes_D(C,D):

    ab = []

    for d in D:

        valores = intervalo(
            lambda x: Bj(x,1/2) - (C + d*x),
            -10,10
            )
        ab.append((d,valores))

    for d,interval in ab:
        for a,b in interval:
            print(f"{d}: \t ({a},{b})")

#raizes_D(C[-1],D)
x0_T = MxT(C[-1],D)
y = Bj(x0_T,1/2)

plt.plot(D,y)
plt.title("Magnetização x Temperatura à campo nulo")
plt.xlabel('kT/2zJS²')
plt.ylabel(r'$M/M_0$')
plt.grid(True)
plt.savefig("./plots/Campo Médio/MxT.png",dpi=300)
plt.show()

#raizes_C(C,D[0])

raizes = MxH(C,D[0])
r1 = []
for c,x0 in raizes:
    r1.append(x0[2])

y = Bj(r1,1/2)
plt.plot(C,y)
plt.show()