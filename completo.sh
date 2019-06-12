#!/bin/bash
	'''ME COMPILA Y EJECUTA EL ARCHIVO EN Single_Particle_Bessel_Circ_Polarization-campo.cpp  Y LLAMA EL HA BASH PARA ORDENAR ARCHIVOS DE SALIDA  '''
	
	c++ Single_Particle_Bessel_Circ_Polarization-campo.cpp -o ejecutable

	./ejecutable

	bash ordenar.sh
