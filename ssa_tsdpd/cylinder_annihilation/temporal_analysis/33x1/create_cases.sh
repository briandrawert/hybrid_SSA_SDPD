#!/bin/bash

counter=1
while [ $counter -le 100 ]
do
	echo Create directory $counter
	mkdir $counter
	cd $counter
	cp ../cylinder_annihilation.lmp .
	mpirun -np 1  ../../../../../../src/lmp_mpi -in cylinder_annihilation.lmp
	cd ..
	((counter++))	

done
echo All done

