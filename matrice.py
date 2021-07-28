import ellipse
import numpy as np
import matplotlib.pyplot as plt
from vitesse import *
dt = 10
T_etude = 2 * 365
dm = 8.64 * 10 ** -5
Tau = 0.66
alpha = 20

"""
CORRECTIONS

- Prise en compte de la porosité (vitesse des pores vs vitesse de Darcy)
- Utilisation de abs pour éviter que la dispersivité soit négative
- Correction du signe devant Ad (par rapport à votre définition de Ad)
- Correction (mais qui n'a pas d'impact dx vs dy dans l'advection)
- Correction Dx vs Dy (dans les noeuds du coeur + 1 des bords de la matrice Ad)

"""

Dx = np.zeros((nx) * (ny))
Dy = np.zeros((nx) * (ny))
Aa = np.zeros((nx * ny, nx * ny))
Ad = np.zeros((nx * ny, nx * ny))
Bd=np.zeros((nx) * (ny))

#u = -0.00001*86400*np.ones(nx*ny)
#v = -0.00001*86400*np.ones(nx*ny)
u = u/0.3
v = v/0.3


for i in range((nx) * (ny)):
    Dx[i] = dm * Tau + alpha * abs(u[i])

for i in range((nx) * (ny)):
    Dy[i] = dm * Tau + alpha * abs(v[i])

for i in range(nx - 1):
    for j in range(ny - 1):
        k = i + nx * j
        Aa[k, k + 1] = -u[k] * (dt / dx)
        Aa[k, k + nx] = -v[k] * (dt / dy)
        Aa[k, k] =  u[k] * (dt / dx) + v[k] * (dt / dy)



for i in range(1, nx - 1):
    for j in range(1, ny - 1):
        k = i + nx * j
        Ad[k, k + 1] = - Dx[k] * dt / dx ** 2
        Ad[k, k - 1] = - Dx[k] * dt / dx ** 2
        Ad[k, k] = Dx[k] * dt / dx ** 2 + Dx[k + 1] * dt / dx ** 2 + Dy[k] * dt / dy ** 2 + Dy[k + nx] * dt / dy ** 2
        Ad[k, k - nx] = - Dy[k] * dt / dy ** 2
        Ad[k, k + nx] = - Dy[k + nx] * dt / dy ** 2

for j in range(1,ny-1):
    i=0
    k=i+nx*j
    Ad[k, k] = + Dy[k] * dt / dy ** 2 + Dy[k + nx] * dt / dy **2 + Dx[k] * dt / dx ** 2
    Ad[k, k + 1] = - Dx[k] * dt / dx ** 2
    Ad[k, k - nx] = - Dy[k] * dt / dy ** 2
    Ad[k, k + nx] = - Dy[k + nx] * dt / dy ** 2

for i in range(1,nx-1):
    j=0
    k=i+nx*j
    Ad[k, k + 1] = - Dx[k] * dt / dx ** 2
    Ad[k, k - 1] = - Dx[k] * dt / dx ** 2
    Ad[k, k] = Dx[k] * dt / dx ** 2 + Dx[k + 1] * dt / dx ** 2 + Dy[k + nx] * dt / dy ** 2
    Ad[k, k + nx] = - Dy[k + nx] * dt / dy ** 2

i = 0
j = 0
k = i + nx * j
Ad[k, k + 1] = - Dx[k] * dt / dx ** 2
Ad[k, k] = Dx[k + 1] * dt / dx ** 2 + Dy[k + nx] * dt / dy ** 2
Ad[k, k + nx] = - Dy[k + nx] * dt / dy ** 2

Ct = ellipse.C
for i in range(T_etude // dt):
    Ct = np.dot((np.eye(nx * ny) + Aa),Ct)
    Ct = np.linalg.solve(np.eye(nx * ny) + Ad,Ct)


    Cplt2 = np.reshape(Ct, (ny, nx))
    plt.imshow(Cplt2, origin='lower',cmap='jet')
    plt.colorbar()
    plt.show()


