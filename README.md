# hybrid_SSA_SDPD
### A hybrid smoothed dissipative particle dynamics (SDPD) spatial stochastic simulation algorithm (sSSA) for advection-diffusion-reaction problems

by:  Brian Drawert, Bruno Jacob, Zhen Li, Tau-Mu Yi, Linda Petzold


We have developed a new algorithm which merges discrete stochastic simulation, using the spatial stochastic simulation algorithm (sSSA), with the particle based fluid dynamics simulation framework of smoothed dissipative particle dynamics (SDPD). This hybrid algorithm enables discrete stochastic simulation of spatially resolved chemically reacting systems on a mesh-free dynamic domain with a Lagrangian frame of reference. SDPD combines two popular mesoscopic techniques: smoothed particle hydrodynamics and dissipative particle dynamics (DPD), linking the macroscopic and mesoscopic hydrodynamics effects of these two methods. We have implemented discrete stochastic simulation using the reaction-diffusion master equations (RDME) formalism, and deterministic reaction-diffusion equations based on the SDPD method. We validate the new method by comparing our results to four canonical models, and demonstrate the versatility of our method by simulating a flow containing a chemical gradient past a yeast cell in a microfluidics chamber.

Keywords: Particle Based Fluid Dynamics, Reaction-Diffusion Master Equation, Discrete Stochastic Simulation
