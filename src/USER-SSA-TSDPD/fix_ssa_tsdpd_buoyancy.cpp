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

#include "string.h"
#include "stdlib.h"
#include "fix_ssa_tsdpd_buoyancy.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "error.h"
#include "update.h"
#include "comm.h"

using namespace LAMMPS_NS; 
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSsaTsdpdBuoyancy::FixSsaTsdpdBuoyancy(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix ssa_tsdpd/buoyancy command");
  
  boussinesq_ssa_flag = boussinesq_sdpd_flag = gravity_flag = 0;

  if (strcmp(arg[3],"boussinesq/ssa") == 0)  boussinesq_ssa_flag = 1;
  else if (strcmp(arg[3],"boussinesq/sdpd") == 0)  boussinesq_sdpd_flag = 1;
  else if (strcmp(arg[3],"gravity") == 0)  gravity_flag = 1;
  else error->all(FLERR,"Illegal type of force in fix ssa_tsdpd/buoyancy command. Valid options: <boussinesq/ssa>, <boussinesq/sdpd> or <gravity>");


  acceleration    = atof(arg[4]);
  rank_coordinate = atoi(arg[5]);
  rank_buoyancy   = atoi(arg[6]);
  C_ref           = atof(arg[7]);


  if (boussinesq_ssa_flag == 1  && atom->num_ssa_species == 0) error->all(FLERR,"Buoyancy force is being used based on SSA species, but number of SSA species is zero");
  if (boussinesq_ssa_flag == 1  && atom->num_ssa_species > 0  ) {
    if(rank_buoyancy >= atom->num_ssa_species) error->all(FLERR,"Index out of bounds: rank of buoyancy concentration term > number SSA species");
    if(atom->concentration_conversion == 1) error->warning(FLERR,"Concentration conversion number not set or being considered unitary. Please specify the conversion factor (molecules -> volumetric concentration) in the atom_style ssa_tsdpd command.");
  }

  xflag = yflag = zflag = 0;
  if (rank_coordinate == 0) xflag = 1;
  else if (rank_coordinate == 1) yflag = 1;
  else if (rank_coordinate == 2) zflag = 1;


  else error->all(FLERR,"Illegal fix ssa_tsdpd/buoyancy command");

  if (xflag  && domain->xperiodic )
    error->all(FLERR,"Cannot use buoyancy force in a periodic dimension");
  if (yflag  && domain->yperiodic)
    error->all(FLERR,"Cannot use buoyancy force in a periodic dimension");
  if (zflag && domain->zperiodic)
    error->all(FLERR,"Cannot use buoyancy force in a periodic dimension");

  MPI_Barrier(world);

}

/* ---------------------------------------------------------------------- */

FixSsaTsdpdBuoyancy::~FixSsaTsdpdBuoyancy()
{

}


/* ---------------------------------------------------------------------- */

int FixSsaTsdpdBuoyancy::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}



/* ---------------------------------------------------------------------- */

void FixSsaTsdpdBuoyancy::setup(int vflag)
{
    post_force(vflag);
}


/* ---------------------------------------------------------------------- */

/*
Syntax of buoyancy term (Boussinesq approximation)
f_buoyancy[i][rank_coordinate] = mass[itype] * acceleration * (C[i][rank_buoyancy] - C_ref)
*/

/*
Syntax of acceleration term (e.g., gravitational force)
f_buoyancy[i][rank_coordinate] = mass[itype] * acceleration
*/

void FixSsaTsdpdBuoyancy::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rho = atom->rho;
  double **C = atom->C;
  int **Cd = atom->Cd;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (boussinesq_sdpd_flag == 1) {
        f[i][rank_coordinate] += mass[type[i]]* acceleration * (C[i][rank_buoyancy] - C_ref);
      }
      else if (boussinesq_ssa_flag == 1) {
        f[i][rank_coordinate] += mass[type[i]]* acceleration * (Cd[i][rank_buoyancy] / atom->concentration_conversion / mass[type[i]]- C_ref) ;
      }
      else if (gravity_flag == 1) {
        f[i][rank_coordinate] += mass[type[i]]* acceleration;
      }

    }
  }
}
