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

   This fix was designed to add source in tDPD
------------------------------------------------------------------------- */
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ssa_tdpd_source.h"
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

FixSsaTdpdSource::FixSsaTdpdSource(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"ssa_tdpd_source") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix SsaTdpdSource command");

  int arg_index = 3;
  step = atoi(arg[arg_index++]);
  typec = atoi(arg[arg_index++]);
  wdim  = atoi(arg[arg_index++]);

  if(typec > atom->num_tdpd_species-1) 
      error->all(FLERR,"Illegal fix ssa_tdpd_source command -- type option wrong.\n");   

  if (strcmp(arg[arg_index],"circle") == 0)       	index = 0;   // circle
  else if (strcmp(arg[arg_index],"rectangle") == 0)     index = 1;   //rectangle
  else if (strcmp(arg[arg_index],"rectangle_set") == 0) index = 2;   //rectangle
  else error->all(FLERR,"Illegal fix ssa_tdpd_source command, choice error");

  arg_index++; //skips the geometry definition (i.e., skips the strings circle or rectangle

  //if circle
  if(index == 0){
    if (narg != 11 ) error->all(FLERR,"Illegal fix ssa_tdpd_source command");
	
	center[0] = atof( arg[arg_index++] );
	center[1] = atof( arg[arg_index++] );
	radius = atof( arg[arg_index++] );
    value  = atof( arg[arg_index++] );
  }
  //if rectangle
  else if(index == 1 || index == 2)  {
    if (narg != 12 ) error->all(FLERR,"Illegal fix ssa_tdpd_source command, ");
	
    center[0] = atof( arg[arg_index++] );
    center[1] = atof( arg[arg_index++] );
    length = atof( arg[arg_index++] );
    width  = atof( arg[arg_index++] );
    value  = atof( arg[arg_index++] );
  }
  else error->all(FLERR,"Illegal fix ssa_tdpd_source command");


  printf("caution:fix_ssa_tdpd_source, cpu is %d\n", comm->me);
  MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

FixSsaTdpdSource::~FixSsaTdpdSource()
{
}

int FixSsaTdpdSource::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSsaTdpdSource::init()
{

}

/* ---------------------------------------------------------------------- */


void FixSsaTdpdSource::post_force(int vflag)
{
  double **x = atom->x;
  double **C = atom->C;
  double **Q = atom->Q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double drx, dry, rsq;
  double radius_sq = radius*radius;


  if(update->ntimestep > step) {
    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            if(index == 0){
                drx = x[i][(wdim+1)%3] - center[0];
                dry = x[i][(wdim+2)%3] - center[1];
                rsq = drx*drx + dry*dry;
                if(rsq < radius_sq)			
                Q[i][typec] += value;
            }
            else if(index == 1){
                drx = x[i][(wdim+1)%3] - center[0];
                dry = x[i][(wdim+2)%3] - center[1];
                if(abs(drx) < length && abs(dry) < width)			
                Q[i][typec] += value;
            }
            else if(index == 2){
                drx = x[i][(wdim+1)%3] - center[0];
                dry = x[i][(wdim+2)%3] - center[1];
                //printf("FixSsaTdpdSource::post_force, index=2 wdim=%i center=%f %f typec=%i length=%f width=%f xy=%f %f dr=%f %f\n",wdim, center[0], center[1], typec, length, width, x[i][(wdim+1)%3], x[i][(wdim+2)%3], drx, dry);
                if(abs(drx) < length && abs(dry) < width)
                C[i][typec] = value;
            }
        }
    }
}
}

/* ---------------------------------------------------------------------- */

