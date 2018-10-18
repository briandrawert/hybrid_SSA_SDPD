OBJ_DIR = ./

UNAME=$(shell uname)

SRCS =  diffusion2d_ssa.lmp

	

ifeq ($(UNAME), Darwin)
BIN = ../../../src/lmp_mac
RUNCMD = DYLD_LIBRARY_PATH=/Users/brian/Desktop/ssa_tDPD/local/lib $(BIN) -in $(SRCS)
else
#NP=$(shell nproc)
NP = 1
BIN = ../../../../src/lmp_mpi
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
