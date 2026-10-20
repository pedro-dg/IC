import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,10,1000)
B = np.linspace(0,10,1000)
T = np.linspace(0,10,1000)
z = np.sinh(8*x)/np.sinh(x)
F = -T*np.log(z)
dT = T[1] - T[0]
s = -np.gradient(F,dT)
plt.plot(x,z)
plt.show()