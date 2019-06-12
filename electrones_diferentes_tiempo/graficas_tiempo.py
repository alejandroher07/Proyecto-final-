'''================================================================================
  ================================================================================
    ESTE PRGRAMA LEE LOS ARCHIVOS DE POSICION DE 1000 ELECTRONES EN UN DETERMINADO  
	TIEMPO Y LOS ALMACENA COMO DATAFRAME.

    GRAFICA LOS ARCHIVOS LEIDOS, SU POSICION EN ESE INSTANTE ADICIONAL DIFERENCIADO 
	POR LA ENERGIA QUE SE TIENE.
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
#	EN LA VARIABLE tiempo, UTILIZANDO PANDAS
#***********************************************

files = glob.glob('*.txt') # IPython magic
tiempo = ([pd.read_csv(f, sep=" ", names=['Posicion x','Posicion y','Posicion z', 'gamma','radio']) for f in files])


##***********************************************
#     FUNCIONES 
#***********************************************

def grafica(x):
    '''ESTA FUNCION REALIZA UN GRAFICO 2D DE LA POSICION 
 	 DEL LOS ELECTRONES EN UN DETERMINADO TIEMPO
	EL INSTANTE CAMBIO DEPENDIENDO DE x'''
    plt.figure(figsize=(15,10))

    plt.tick_params(labelsize = 20)
    plt.xlabel('Longitud de la cavidad [cm]',fontsize=20)
    plt.ylabel('Radio [cm]',fontsize=20)
    plt.grid(True)

    colors = tiempo[x]['gamma']
    R = np.sqrt(tiempo[x]['Posicion x']*tiempo[x]['Posicion x']+tiempo[x]['Posicion y']*tiempo[x]['Posicion y'])

    plt.scatter(tiempo[x]['Posicion z'],R,c=colors)
    plt.savefig('tiempo{0}.jpg'.format(x))

##***********************************************
#     CICLO PARA PASAR POR TODOS LOS ARCHIVOS  
#***********************************************
for i in range(len(files)):
    grafica(i)

