#!/bin/bash
	'''ESTE BASH ME ORDENA LOS ARCHIVOS EN LA CARPTA CORRESPONDIENTE '''

	#mkdir datos	
	#mkdir testing
	#mkdir Posiciones
	
	for i in resul*.txt;
	do
	mv $i ./datos/;
	done

	for i in Posi*.txt;
	do
	mv $i ./electrones_diferentes_tiempo/;
	done

	for i in B_*.txt;
	do
	mv $i ./testing/;
	done

	mv Electron_test.txt ./testing/;
