"""
Created on Wed Aug 17 13:06:57 2022

@author: Konstantinos Dalampekis
"""

import numpy as np
from numpy import genfromtxt

import matplotlib.pyplot as plt

L = genfromtxt('/Volumes/Mac Os - Δεδομένα/Python_Projects/Fibre_Laser/Kostas/L.txt', delimiter=',')
print(np.shape(L))
wavelength = genfromtxt('/Volumes/Mac Os - Δεδομένα/Python_Projects/Fibre_Laser/Kostas/wavelength.csv', delimiter=',')
print(np.shape(wavelength))
Pa = genfromtxt('/Volumes/Mac Os - Δεδομένα/Python_Projects/Fibre_Laser/Kostas/Pa.txt', delimiter=',')
print(np.shape(Pa))

Pa = genfromtxt('/Volumes/Mac Os - Δεδομένα/Python_Projects/Fibre_Laser/Kostas/Pa.txt', delimiter=',').reshape((len(wavelength),len(L)))
Pn = genfromtxt('/Volumes/Mac Os - Δεδομένα/Python_Projects/Fibre_Laser/Kostas/Pn.txt', delimiter=',').reshape((len(wavelength),len(L)))
print(np.shape(Pn))

PaLog=10*np.log(Pa/Pa[0,0])
PnLog=10*np.log(Pn/Pn[0,0])
PTLog=10*np.log((Pa+Pn)/(Pa[0,0]+Pn[0,0]))

# plt.figure()
X, Y = np.meshgrid(L, wavelength)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, PaLog)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Pump Pa Absorption/Gain in dB')
ax.set_xlabel('Length (m)')
ax.set_ylabel('Wavelength (nm)')
plt.show()

# plt.figure()
X, Y = np.meshgrid(L, wavelength)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, PnLog)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Pump Pn Absorption/Gain in dB')
ax.set_xlabel('Length (m)')
ax.set_ylabel('Wavelength (nm)')
plt.show()

# plt.figure()
X, Y = np.meshgrid(L, wavelength)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, PTLog)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Pump PTotal Absorption/Gain in dB')
ax.set_xlabel('Length (m)')
ax.set_ylabel('Wavelength (nm)')
plt.show()