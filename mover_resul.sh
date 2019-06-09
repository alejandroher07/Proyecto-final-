#!/bin/bash
	# copia los archivos resultados y los pone en la carpeta datos
	for i in resul*.txt;
	do
	mv $i ./datos/;
	done
	

