import numpy as np
from Epaisseur import*
import matplotlib.pyplot as plt

Lx = 1000
Ly = 1200
dx = 20
dy = 20
xc = 800
yc = 1000
a = 100
b = 80
k1=25.92
rhau=10**3
q=216
nx = int(Lx // dx + 1)
ny = int(Ly // dy + 1)
#T = np.zeros((nx) * (ny))
Ae=np.eye((nx * ny))
Be=np.zeros((nx) * (ny))
ep=calcule_epaisseur(nx,ny)

T = k1 * ep

"""
for i in range(1,nx - 1):
    for j in range(1,ny - 1):
        k=i+nx*j
        Ae[k,k+1]=T[k+1]*(1/dx**2)
        Ae[k,k]=-(1/dx**2*T[k+1]+T[k+nx]*1/dy**2+T[k]*(1/dx**2+1/dy**2))
        Ae[k,k-1]=T[k]*(1/dx**2)
        Ae[k,k-nx]=T[k]*(1/dy**2)
        Ae[k,k+nx]=T[k+nx]*(1/dy**2)


for j in range(ny-1):
    i=nx-1
    k=i+nx*j
#    Ae[k,k]=1
#    Ae[k][k-1]=-1

    Ae[k, k] = -(T[k + nx] * 1 / dy ** 2 + T[k] * (1 / dx ** 2 + 1 / dy ** 2))
    Ae[k,k-1]=T[k]*(1/dx**2)
    Ae[k, k - nx] = T[k] * (1 / dy ** 2)
    Ae[k, k + nx] = T[k + nx] * (1 / dy ** 2)
    Be[k]=-q*T[k]*1/(rhau*k1*dx)


for i in range(0,nx-1):
    j=ny-1
    k=i+nx*j
#    Ae[k][k]=1
#    Ae[k][k-nx]=-1

    Ae[k, k + 1] = T[k + 1] * (1 / dx ** 2)
    Ae[k, k] = -(1 / dx ** 2) * T[k + 1] + T[k] * (1 / dx ** 2  + 1 / dy ** 2)
    Ae[k,k-1]=T[k]*(1/dx**2)
    Ae[k, k - nx] = T[k] * (1 / dy ** 2)
    Be[k]=-q*T[k]*1/(rhau*k1*dx)
i=nx-1
j=ny-1
k=nx*j+i

Ae[k, k] = -( T[k] * (1 / dx ** 2 + 1 / dy ** 2))
Ae[k, k - 1] = T[k] * (1 / dx ** 2)
Ae[k, k - nx] = T[k] * (1 / dy ** 2)
Be[k] = -q * T[k] * 1 / (rhau * k1 * dx) - q * T[k] * 1 / (rhau * k1 * dx)
"""

for i in range(1,nx - 1):
    for j in range(1,ny - 1):
        k=i+nx*j
        Ae[k,k+1]=T[k+1]*1/dx**2
        Ae[k,k]=-(T[k+1]*1/dx**2+T[k+nx]*1/dy**2+T[k]*(1/dx**2+1/dy**2))
        Ae[k,k-1]=T[k]*1/dx**2
        Ae[k,k-nx]=T[k]*1/dy**2
        Ae[k,k+nx]=T[k+nx]*1/dy**2


for j in range(1,ny-1):
    i=nx-1
    k=i+nx*j
#    Ae[k,k]=1
#    Ae[k][k-1]=-1

    Ae[k, k] = -(T[k + nx]*1/dy**2  + T[k]*(1/dx**2+1/dy**2))
    Ae[k, k - nx] = T[k]*1/dy**2
    Ae[k, k - 1] = T[k]*1/dx**2
    Ae[k, k + nx] = T[k + nx]*1/dy**2
    Be[k]=-q*T[k]/(rhau*k1*dx)


for i in range(1,nx-1):
    j=ny-1
    k=i+nx*j
#    Ae[k][k]=1
#    Ae[k][k-nx]=-1

    Ae[k, k + 1] = T[k + 1]*1/dx**2
    Ae[k, k] = - T[k + 1]*1/dx**2 - T[k]*(1/dx**2+1/dy**2)
    Ae[k,k-1]=T[k]*1/dx**2
    Ae[k, k - nx] = T[k]*1/dy**2
    Be[k]=-q*T[k]/(rhau*k1*dy)
i=nx-1
j=ny-1
k=nx*j+i

Ae[k, k] = -( 2*T[k] )
Ae[k, k - 1] = T[k]
Ae[k, k - nx] = T[k]
Be[k] = -q * T[k] / (rhau * k1*dx ) - q * T[k]  / (rhau * k1*dy )

h=np.linalg.solve(Ae,Be)
h1= np.reshape(h, (ny, nx))
plt.imshow(h1, origin='lower',cmap='jet')
plt.colorbar()
plt.show()












