# -*- coding:utf-8 -*-

__projet__ = "modl"
__nom_fichier__ = "Epaisseur"
__author__ = "Prénom Nom"
__date__ = "octobre 2019"

import numpy as np
import math as mt
dt = 364
dx = 10
dy = 10
nx = int(1000 // dx + 1)
ny = int(1200 // dy + 1)
vx = -0.864
T_etude = 2 * 365
vy = -0.864
dm = 8.64 * 10 ** -5
Tau = 0.66
alpha = 20

# Définition de l'épaisseur (coordonnées globales)

def calcule_epaisseur(nx,ny):
    """Décrit le profil des épaisseurs."""

    ep = np.zeros(nx*ny)
    for i in range(nx):
        for j in range(ny):
            k = i+j*nx
            r = 1-mt.sqrt((2*i - nx)**2 + (2*j - nx)**2)/mt.sqrt(nx**2 + ny**2)
            ep[k] = -25*r**2 + 30
    return ep


