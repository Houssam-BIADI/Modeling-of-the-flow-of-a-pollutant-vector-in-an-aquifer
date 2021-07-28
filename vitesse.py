# -*- coding:utf-8 -*-

__projet__ = "modl"
__nom_fichier__ = "vitesse"
__author__ = "Prénom Nom"
__date__ = "octobre 2019"
from charge_hydr import *
u = np.zeros((nx) * (ny))
v = np.zeros((nx) * (ny))
for i in range(nx-1):
    for j in range(ny-1):
        k=i+nx*j
        u[k]=-k1*(h[k+1]-h[k])*(1/dx)
        v[k] = -k1 * (h[k + nx] - h[k]) * (1 / dy)

for j in range(ny-1):
    i=nx-1
    k = i + nx * j
    u[k]=-q/rhau
    v[k] = -k1 * (h[k + nx] - h[k]) * (1 / dy)

for i in range(nx-1):
      j=ny-1
      k = i + nx * j
      v[k]=-q/rhau
      u[k] = -k1 * (h[k + 1] - h[k]) * (1 / dx)

j = ny-1
i = nx-1
k = i + nx * j
u[k] = -q / rhau
v[k]=-q/rhau

u1= np.reshape(u, (ny, nx))
plt.imshow(u1, origin='lower',cmap='jet')
plt.colorbar()
plt.show()
v1= np.reshape(v, (ny, nx))
plt.imshow(v1, origin='lower',cmap='jet')
plt.colorbar()
plt.show()

