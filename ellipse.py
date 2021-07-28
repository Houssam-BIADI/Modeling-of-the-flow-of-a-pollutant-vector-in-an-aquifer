import numpy as np
import matplotlib.pyplot as plt

# from matrice import *

Lx = 1000
Ly = 1200
dx = 20
dy = 20
xc = 800
yc = 1000
a = 100
b = 80
nx = int(Lx // dx + 1)
ny = int(Ly // dy + 1)

C = np.zeros(nx * ny)

for i in range(nx):
    for j in range(ny):
        k = i + nx * j

        x = i * dx
        y = j * dy

        if (x - xc) ** 2 / a ** 2 + (y - yc) ** 2 / b ** 2 <= 1:
            C[k] = 150
        else:
            C[k] = 5

Cplt = np.reshape(C, (ny, nx))
plt.imshow(Cplt, origin='lower',cmap='jet')
plt.show()
