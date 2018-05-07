# hybrid_SSA_SDPD
### A hybrid smoothed dissipative particle dynamics (SDPD) spatial stochastic simulation algorithm (sSSA) for advection-diffusion-reaction problems

by:  Brian Drawert, Bruno Jacob, Zhen Li, Tau-Mu Yi, Linda Petzold


We have developed a new algorithm which merges discrete stochastic simulation, using the spatial stochastic simulation algorithm (sSSA), with the particle based fluid dynamics simulation framework of smoothed dissipative particle dynamics (SDPD). This hybrid algorithm enables discrete stochastic simulation of spatially resolved chemically reacting systems on a mesh-free dynamic domain with a Lagrangian frame of reference. SDPD combines two popular mesoscopic techniques: smoothed particle hydrodynamics and dissipative particle dynamics (DPD), linking the macroscopic and mesoscopic hydrodynamics effects of these two methods. We have implemented discrete stochastic simulation using the reaction-diffusion master equations (RDME) formalism, and deterministic reaction-diffusion equations based on the SDPD method. We validate the new method by comparing our results to four canonical models, and demonstrate the versatility of our method by simulating a flow containing a chemical gradient past a yeast cell in a microfluidics chamber.

Keywords: Particle Based Fluid Dynamics, Reaction-Diffusion Master Equation, Discrete Stochastic Simulation



To run the examples found in the paper, you can use our pre-configued docker containers.  First make sure you have [Docker installed](https://www.docker.com/) and running.  Then type the following commands into the terminal to run each of the examples.  The result will be to create a number of VTK files in the current directory.  You can view these files with a program such as [ParaView](https://www.paraview.org/).

From the main text: 

### Validation for one-dimensional reaction-diffusion:
```
docker run -it -v "`pwd`":/work  briandrawert/hybrid_ssa_sdpd /bin/bash -c "/run_1d_rxn_diffusion.sh"
```

### Validation for natural convection in a cylinder inside a square enclosure:
```
docker run -it -v "`pwd`":/work  briandrawert/hybrid_ssa_sdpd /bin/bash -c "/run_natural_convection.sh"
```

### Application: micro-channel reactive flow past a cell:
```
docker run -it -v "`pwd`":/work  briandrawert/hybrid_ssa_sdpd /bin/bash -c "/run_cell_example.sh"
```

From the SI:

### Validation for one-dimensional diffusion:
```
docker run -it -v "`pwd`":/work  briandrawert/hybrid_ssa_sdpd /bin/bash -c "/run_1d_diffusion.sh"
```

### Validation for two-dimensional diffusion:
```
docker run -it -v "`pwd`":/work  briandrawert/hybrid_ssa_sdpd /bin/bash -c "/run_2d_diffusion.sh"
```

