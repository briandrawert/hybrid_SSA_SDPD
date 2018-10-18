OBJ_DIR = ./

UNAME=$(shell uname)

SRCS =  cell_model.lmp

	

ifeq ($(UNAME), Darwin)
BIN = ../../../src/lmp_mac
RUNCMD = DYLD_LIBRARY_PATH=/Users/brian/Desktop/ssa_tDPD/local/lib $(BIN) -in $(SRCS)
else
#NP=$(shell nproc)
NP = 1
BIN = ../../../../../src/lmp_serial
#RUNCMD = xterm -hold -e gdb -ex run $(BIN) < $(SRCS)
RUNCMD = $(BIN) < $(SRCS)
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
