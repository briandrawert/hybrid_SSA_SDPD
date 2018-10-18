#!/bin/bash

counter=1
while [ $counter -le 100 ]
do
	echo Create directory $counter
	mkdir $counter
	cd $counter
	cp ../diffusion1d_ssa.lmp .
	mpirun -np 1  ../../../../../src/lmp_mpi -in diffusion1d_ssa.lmp
	cd ..
	((counter++))	

done
echo All done

