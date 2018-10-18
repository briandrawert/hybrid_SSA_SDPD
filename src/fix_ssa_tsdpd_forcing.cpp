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

   This fix was designed to act as a forcing zone in tDPD
------------------------------------------------------------------------- */
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ssa_tsdpd_forcing.h"
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

FixSsaTsdpdForcing::FixSsaTsdpdForcing(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"ssa_tsdpd_forcing") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix SsaTsdpdForcing command, first error.");


  ssa_case = 0;
  tsdpd_case = 0;
  velocity_case = 0;

  int arg_index = 3;
 
  if (strcmp(arg[arg_index],"ssa") == 0)  ssa_case = 1;
  else if (strcmp(arg[arg_index],"tsdpd") == 0)  tsdpd_case = 1;
  else if (strcmp(arg[arg_index],"velocity") == 0) velocity_case = 1;
  else error->all(FLERR,"Illegal argument[3]. Choose <ssa>, <tsdpd> or <velocity>");
  
  arg_index++;
  step = atoi(arg[arg_index++]);

  if (ssa_case==1) {
    ctype = atoi(arg[arg_index++]);
    if(ctype > atom->num_ssa_species) 
        error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command: species id > num_ssa_species.\n");   
  }
  else if (tsdpd_case==1) {
     ctype = atoi(arg[arg_index++]);
     if(ctype > atom->num_tdpd_species)  
        error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command: species id > num_tsdpd_species.\n");   
  }
  else if (velocity_case==1) { 
     vtype = atoi(arg[arg_index++]);
     if(vtype > 3)  
        error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command: velocity id must be either 1 (vx), 2 (vy)  or 3 (vz).\n");   
  }

    if (strcmp(arg[arg_index],"circle") == 0)       	index = 0;
    else if (strcmp(arg[arg_index],"rectangle") == 0)     index = 1;
    else error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command, symbol error.");
    arg_index++;

    if(index == 0){
        if (narg != 11 ) error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command, index0");
	
	center[0] = atof( arg[arg_index++] );
	center[1] = atof( arg[arg_index++] );
	radius = atof( arg[arg_index++] );
	if (tsdpd_case == 1) value  = atof( arg[arg_index++] );
	else if (ssa_case == 1) value_int = int(atof( arg[arg_index++] ));
	else if (velocity_case == 1) value = atof( arg[arg_index++] );
    }
    else if(index == 1){
        if (narg != 12 ) error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command, index1");
	
	center[0] = atof( arg[arg_index++] );
	center[1] = atof( arg[arg_index++] );
	length = atof( arg[arg_index++] );
	width  = atof( arg[arg_index++] );
	if (tsdpd_case == 1) value  = atof( arg[arg_index++] );
	else if (ssa_case == 1) value_int = int(atof( arg[arg_index++] ));
	else if (velocity_case == 1) value = atof( arg[arg_index++] );
    }
    else error->all(FLERR,"Illegal fix ssa_tsdpd_forcing command, index error.");


  MPI_Barrier(world);


}

/* ---------------------------------------------------------------------- */

FixSsaTsdpdForcing::~FixSsaTsdpdForcing()
{
}

/* ---------------------------------------------------------------------- */

int FixSsaTsdpdForcing::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpdForcing::init()
{

}

/* ---------------------------------------------------------------------- */


void FixSsaTsdpdForcing::post_integrate()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double **C = atom->C;
  double **Q = atom->Q;
  int **Cd = atom->Cd;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double drx, dry, rsq;
  double radius_sq = radius*radius;

  if(update->ntimestep > step) 
  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
		if(index == 0){
		    drx = x[i][0] - center[0];
		    dry = x[i][1] - center[1];
		    rsq = drx*drx + dry*dry;
		    if(rsq < radius_sq) {
			if (tsdpd_case==1) C[i][ctype] = value;
			if (ssa_case==1) Cd[i][ctype] = value_int;
			if (velocity_case==1) v[i][vtype] = value;
                        
		    }
		}
		else if(index == 1){
		    drx = x[i][0] - center[0];
		    dry = x[i][1] - center[1];
		    if(fabs(drx) < length && fabs(dry) < width){
			if (tsdpd_case==1) C[i][ctype] = value; 
			if (ssa_case==1) Cd[i][ctype] = value_int;
			if (velocity_case==1) v[i][vtype] = value;
		    }
		}
    	}

  }
}

/* ---------------------------------------------------------------------- */

