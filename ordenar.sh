#!/bin/bash
	# Ordena los archivos creados por el codigo .cpp
	#mkdir datos	
	#mkdir testing
	#mkdir Posiciones
	
	for i in resul*.txt;
	do
	mv $i ./datos/;
	done

	for i in Posi*.txt;
	do
	mv $i ./Posiciones/;
	done

	for i in B_*.txt;
	do
	mv $i ./testing/;
	done


