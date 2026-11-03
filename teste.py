import numpy as np
import matplotlib.pyplot as plt

B = np.linspace(0.1,10,1000)[:, np.newaxis]
T = np.linspace(0.1,10,1000)
z = np.sinh(8*B/T)/np.sinh(B/T)
F = -T*np.log(z)
dT = T[1] - T[0]
s = -np.gradient(F,dT)
x = B/T
plt.plot(x,z)
plt.show()