UNAME=$(shell uname)


OBJ_DIR = ./

SRCS = 	cylinder_annihilation.lmp


ifeq ($(UNAME), Darwin)
BIN = ../../../src/lmp_mac
RUNCMD = DYLD_LIBRARY_PATH=/Users/brian/Desktop/ssa_tDPD/local/lib $(BIN) -in $(SRCS)
else
NP= 1 #$(shell nproc)
BIN = ../../../src/lmp_mpi
RUNCMD = mpirun -np $(NP) $(BIN) -in $(SRCS)
endif
	
run :	
	clear
	$(RUNCMD)

clean :
	rm -rf *.vtk
	rm -rf log.*
	rm -rf *.lammpstrj
	rm -rf *.profile
	rm -rf *.dat
