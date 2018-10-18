/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Zhen Li (zhen_li@brown.edu)
   the CRUNCH group, Division of Applied Mathematics, Brown University

   This fix was designed to reset the concentration in tDPD
------------------------------------------------------------------------- */
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ssa_tdpd_resetc.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "iostream"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSsaTdpdResetC::FixSsaTdpdResetC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"cc_resetc") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix SsaTdpdResetC command, first error.");

  int arg_index = 3;
  step = atoi(arg[arg_index++]);
  ctype = atoi(arg[arg_index++]);

  if(ctype > atom->num_tdpd_species) //added (must add num_ssa_species later) 
      error->all(FLERR,"Illegal fix cc_resetc command -- type option wrong.\n");   

  if (strcmp(arg[arg_index],"circle") == 0)       	index = 0;
  else if (strcmp(arg[arg_index],"rectangle") == 0)     index = 1;
  else error->all(FLERR,"Illegal fix cc_resetc command, symbol error.");
  arg_index++;

  if(index == 0){
        if (narg != 10 ) error->all(FLERR,"Illegal fix cc_resetc command, index0");
	
	center[0] = atof( arg[arg_index++] );
	center[1] = atof( arg[arg_index++] );
	radius = atof( arg[arg_index++] );
        value  = atof( arg[arg_index++] );
  }
  else if(index == 1){
        if (narg != 11 ) error->all(FLERR,"Illegal fix cc_resetc command, index1");
	
	center[0] = atof( arg[arg_index++] );
	center[1] = atof( arg[arg_index++] );
	length = atof( arg[arg_index++] );
	width  = atof( arg[arg_index++] );
        value  = atof( arg[arg_index++] );
  }
  else error->all(FLERR,"Illegal fix cc_resetc command, index error.");


  printf("caution:fix_cc_resetc, cpu is %d\n", comm->me);
  MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

FixSsaTdpdResetC::~FixSsaTdpdResetC()
{
}

/* ---------------------------------------------------------------------- */

int FixSsaTdpdResetC::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSsaTdpdResetC::init()
{

}

/* ---------------------------------------------------------------------- */


void FixSsaTdpdResetC::post_integrate()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **C = atom->C;
  double **Q = atom->Q;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double drx, dry, rsq;
  double radius_sq = radius*radius;

  if(update->ntimestep==step) 
  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
		C[i][ctype] = 0.0;
		if(index == 0){
		    drx = x[i][0] - center[0];
		    dry = x[i][1] - center[1];
		    rsq = drx*drx + dry*dry;
		    if(rsq < radius_sq)			
			C[i][ctype]= value;
		}
		else if(index == 1){
		    drx = x[i][0] - center[0];
		    dry = x[i][1] - center[1];
		    if(abs(drx) < length && abs(dry) < width)			
			C[i][ctype]= value;
		}
    	}

  }
}

/* ---------------------------------------------------------------------- */

