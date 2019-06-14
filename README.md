El objetivo de este trabajo es demostrar  cómo se puso en práctica las herramientas trabajadas en la asignatura. Por lo cual se utilizo un archivo Single_Particle_Bessel_Circ_Polarization-campo.cpp para obtener los datos de la trayectoria de electrones los cuales son acelerados en un dispositivo Gyrac, como primera etapa se utilizo un ordenar.sh para organizar los datos y ubicarlos en sus receptivas carpetas. Posteriormente se realizaron diferentes .py para hacer de una forma sencilla gráficas para estos archivos utilizando las herramientas que pandas tiene a su disposición.  

Finalmente con los datos obtenidos se realizo un análisis resaltando los resultados interesante.

# Titulo del trabajo:

## ESTUDIO COMPUTACIONAL DE LA ACELERACIÓN AUTO-RESONANTE DE ELECTRONES EN CAMPOS MAGNÉTICOS VARIABLES EN EL TIEMPO


### Teoria: 
K. S. Golovanisky. En 1982 propuso este esquema. El problema de este mecanismo es el control que se debe tener del campo magnético, ya que debe crecer en una proporción adecuada para mantener la  condición ECR. Debido a esta dificultad y al hecho que este esquema fue diseñado para la aceleración 2D de electrones, su implementación experimental se logró apenas en el año 2017, utilizando un campo magnético no homogéneo pero creciente en el tiempo; donde un plasma confinado y sometido a pulsos de campo magnético de duración del orden de los microsegundos, se logran acelerar los electrones hasta un valor energético suficiente para la producción de rayos X.
# Mecanismo Gyrac

Un electrón en el mecanismo Gyrac es acelerado utilizando un campo magnético que crece suavemente
en el tiempo para mantener la condición de resonancia a pesar del crecimiento del factor relativista a
medida que la partı́cula gana energı́a.


![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/Screenshot%20from%202019-06-12%2014-43-39.png) 



# Campo de microondas
Las componentes de los campos eléctrico y magnético de un modo cilı́ndrico $TE_{111}$ polarizado circularmente excitado en la cavidad

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/campo.png) 

# Campo magnético no homogéneo

Para crear un campo magnético variable en el tiempo se ubican 2 bobinas en los extremos de una
cavidad resonante cilı́ndrica con corriente constante en el mismo sentido, y dos boninas centrales con
corriente variable en sentido opuesto a las dos boninas ya mencionadas, creando un pozo en la región
central. Las bobinas centrales funcionan como bobinas de control para aumentar el campo magnético
correctamente en el tiempo. Los electrones son acelerados en la parte central, donde se tiene mayor
control sobre las condiciones para el sostenimiento de la aceleración auto-resonante. Este
diseño fue propuesto por andreev y colaboradores en el 2012 y es utilizado en la actualidad para la
generación del campo magnético no homogéneo en un dispositivo de Gyrac de forma experimental.


![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/cavidad_corriente_1.jpg) 

# Organización de la información



![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/Screenshot%20from%202019-06-12%2014-50-55.png) 

# Testing

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/campo_0.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/campo_1.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/test_ele_gamma.jpg) 

# Electrones en diferentes instantes de tiempo

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/tiempo0.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/tiempo2.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/tiempo1.jpg) 

# Resultados de interés

Se leen los datos como dataframe utilizando pandas 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/Screenshot%20from%202019-06-12%2015-13-46.png) 

# Resultados de interés

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/gamma_0.jpg) 
![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/trayectoria_0.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/gamma_3.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/trayectoria_8.jpg) 

![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/gamma_7.jpg)


![alt text](https://github.com/alejandroher07/Proyecto-final-/blob/master/imagenes/trayectoria_9.jpg) 
















