import pandas as pd
from pandas import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import glob

Bf = np.loadtxt('B_rz_f.txt')
Bi = np.loadtxt('B_rz_i.txt')
print(Bf.shape)

def grafica(x):
    fig = plt.figure(figsize=(13,13))
    ax = fig.add_subplot(111, projection="3d")
    X, Y = np.mgrid[0:100:1, 0:129:1]
    Z = x
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel('Puntos de malla x',fontsize=20)
    ax.set_ylabel('Puntos de malla y ',fontsize=20)
    ax.set_zlabel('Campo magnetico normalizado ',fontsize=20)
    ax.plot_surface(X, Y, Z, cmap="gnuplot", lw=0.5, rstride=1, cstride=1)
    plt.savefig('campo_{0}.jpg'.format(i))

i=0
grafica(Bi)
i=1
grafica(Bf)


Bif = np.loadtxt('B_z_i_f.txt')

plt.figure(figsize=(15,10))
plt.plot(Bif[:,0], Bif[:,1],Bif[:,2],color="blue", linewidth=1.5,)

plt.tick_params(labelsize = 20)
plt.xlabel('Puntos de malla',fontsize=20)
plt.ylabel('Campo magnetico normalizado',fontsize=20)
plt.title('Corte longitudinal del campo inicial y final',fontsize=20)
plt.grid(True)
plt.savefig('Campo_i_f.jpg')


e = np.loadtxt('Electron_test.txt')


plt.figure(figsize=(15,10))
plt.plot(e[:,0], e[:,3],color="blue", linewidth=1.5,)

plt.tick_params(labelsize = 20)
plt.xlabel('tiempo',fontsize=20)
plt.ylabel('posicion z',fontsize=14)
plt.title('z constante',fontsize=20)
plt.grid(True)
plt.savefig('test_ele_z.jpg')

plt.figure(figsize=(15,10))
plt.plot(e[:,0], e[:,4],color="blue", linewidth=1.5,)

plt.tick_params(labelsize = 20)
plt.xlabel('tiempo',fontsize=20)
plt.ylabel('r',fontsize=20)
plt.grid(True)
plt.savefig('test_ele_gamma.jpg')

plt.figure(figsize=(15,10))
plt.plot(e[:,1], e[:,2],color="blue", linewidth=1.5,)

plt.tick_params(labelsize = 20)
plt.xlabel('Posicion x',fontsize=20)
plt.ylabel('Posicion y ',fontsize=20)
plt.title('Electron en centro de la cavidad',fontsize=20)
plt.grid(True)
plt.savefig('test_ele_2d.jpg')




