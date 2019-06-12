'''================================================================================
  ================================================================================
    ESTE PRGRAMA LEE LOS ARCHIVOS DE LA TRAYECTORIA PARA 100 ELECTRONES Y LOS 
	UTILIZA COMO DATAFRAME.

    GRAFICA LOS ARCHIVOS LEIDOS, SU EVOLUCION ESPACIAL Y ENERGETICO. PARA TODA LA 
	ESCALA DE TIEMPO.
  ================================================================================
  ================================================================================'''

##***********************************************
#     SE IMPORTAN LAS LIBRERIAS NECESARIAS 
#***********************************************

import pandas as pd
from pandas import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import glob #LIBRERIA PARA LEER TODOS LOS .txt


##***********************************************
#     LECTRUTA DE LOS ARCHIVOS ALMACENANDOLOS 
#	EN LA VARIABLE ele, UTILIZANDO PANDAS
#***********************************************

files = glob.glob('*9.txt') 
print(files)
ele = ([pd.read_csv(f, sep=" ", names=['tiempo','Posicion x','Posicion y','Posicion z', 'gamma']) for f in files])


##***********************************************
#     FUNCIONES 
#***********************************************

def gm(x):
    '''ESTA FUNCION REALIZA UN GRAFICO 2D DE LA EVOLUCION 
	ENEGETICA  DEL ELECTRON SELECCIONADO CON LA VARIAIBLE x'''
    plt.figure(figsize=(15,10))
    
    plt.tick_params(labelsize = 20)
    plt.xlabel('tiempo (μs)',fontsize=20)
    plt.ylabel('ɣ',fontsize=20)
    plt.title('Evolucion del factor ɣ para el electron {0}'.format(x+1),fontsize=20)
    plt.grid(True)

    plt.plot(ele[x]['tiempo'], ele[x]['gamma'],color="blue", linewidth=1.5,)
    plt.savefig('graficas_gamma/gamma_{0}.jpg'.format(x))


def tray(x):
     '''ESTA FUNCION REALIZA UN GRAFICO 2D DE LA EVOLUCION 
	ESPACIAL DEL ELECTRON SELECCIONADO CON LA VARIAIBLE x'''   
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

##***********************************************
#     CICLO PARA PASAR POR TODOS LOS ARCHIVOS  
#***********************************************

for i in range(len(files)):
	tray(i)
    	gm(i)
		
    
