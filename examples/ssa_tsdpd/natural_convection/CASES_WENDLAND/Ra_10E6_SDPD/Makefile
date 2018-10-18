0OBJ_DIR = ./

UNAME=$(shell uname)

SRCS =  natural_convection_heated_cylinder.lmp

NP = 16
BIN = ../../../../src/lmp_mpi
RUNCMD = mpirun -np $(NP) -mca btl ^openib $(BIN) -in $(SRCS)

run :	
	clear
	$(RUNCMD)

clean :
	rm -rf *.vtk
	rm -rf log.*
	rm -rf *.lammpstrj
	rm -rf *.profile
	rm -rf *.dat
