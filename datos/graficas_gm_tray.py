import pandas as pd
from pandas import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import glob

files = glob.glob('*.txt') # IPython magic
print(files)
ele = ([pd.read_csv(f, sep=" ", names=['tiempo','Posicion x','Posicion y','Posicion z', 'gamma']) for f in files])

def gm(x):
    
    plt.figure(figsize=(15,10))
    
    plt.tick_params(labelsize = 20)
    plt.xlabel('tiempo (μs)',fontsize=20)
    plt.ylabel('ɣ',fontsize=20)
    plt.title('Evolucion del factor ɣ para el electron {0}'.format(x+1),fontsize=20)
    plt.grid(True)

    plt.plot(ele[x]['tiempo'], ele[x]['gamma'],color="blue", linewidth=1.5,)
    plt.savefig('graficas_gamma/gamma_{0}.jpg'.format(x))


def tray(x):
    
    fig = plt.figure(figsize=(13,13))
    ax = fig.add_subplot(111, projection="3d")

    zline = ele[x]['Posicion x']
    xline = ele[x]['Posicion y']
    yline = ele[x]['Posicion z']

    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel('Posicion x',fontsize=20)
    ax.set_ylabel('Posicion y ',fontsize=20)
    ax.set_zlabel('longitud de la cavidad',fontsize=20)
    ax.plot3D(xline, yline, zline, 'r')

    plt.savefig('graficas_trayectoria/trayectoria_{0}.jpg'.format(x))
print('hola-0')
for i in range(len(files)):
	tray(i)
	print('hola')
    #gm(i)
		
    
