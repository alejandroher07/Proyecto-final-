import pandas as pd
from pandas import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import glob

files = glob.glob('*.txt') # IPython magic

tiempo = ([pd.read_csv(f, sep=" ", names=['Posicion x','Posicion y','Posicion z', 'gamma','radio']) for f in files])

def grafica(x):

    plt.figure(figsize=(15,10))

    plt.tick_params(labelsize = 20)
    plt.xlabel('Longitud de la cavidad [cm]',fontsize=20)
    plt.ylabel('Radio [cm]',fontsize=20)
    plt.grid(True)

    colors = tiempo[x]['gamma']
    R = np.sqrt(tiempo[x]['Posicion x']*tiempo[x]['Posicion x']+tiempo[x]['Posicion y']*tiempo[x]['Posicion y'])

    plt.scatter(tiempo[x]['Posicion z'],R,c=colors)
    plt.savefig('tiempo{0}.jpg'.format(x))

for i in range(len(files)):
    grafica(i)

