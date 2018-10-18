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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fix_etd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "comm.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixETD::FixETD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;

}

/* ---------------------------------------------------------------------- */

int FixETD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixETD::init()
{

random = NULL;	
MPI_Barrier(world);
int seed=int(time(0));
random = new RanMars(lmp, (seed+ comm->me) % 900000000 );



}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixETD::initial_integrate(int vflag)
{

  // update v and x of atoms in group

  double dtf = update->dt;
  double dtv = update->dt;
  double dtfm;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **Aetd = atom->Aetd;
  double **Betd = atom->Betd;
  double **Cetd = atom->Cetd;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}



/* ---------------------------------------------------------------------- */

void FixETD::post_integrate()
{

  // update v of atoms in group
 
  double dtf = update->dt;
  double dtv = update->dt;
  double dtfm;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **Aetd = atom->Aetd;
  double **Betd = atom->Betd;
  double **Cetd = atom->Cetd;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double randnum;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
	v[i][0] = v[i][0]*exp(Cetd[i][0]*dtv) + (exp(Cetd[i][0]*dtv) - 1.0) * (Aetd[i][0] + Betd[i][0] ) / Cetd[i][0];
        v[i][1] = v[i][1]*exp(Cetd[i][1]*dtv) + (exp(Cetd[i][1]*dtv) - 1.0) * (Aetd[i][1] + Betd[i][1] ) / Cetd[i][1];
        v[i][2] = v[i][2]*exp(Cetd[i][2]*dtv) + (exp(Cetd[i][2]*dtv) - 1.0) * (Aetd[i][2] + Betd[i][2] ) / Cetd[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
 	v[i][0] = v[i][0]*exp(Cetd[i][0]*dtv) + (exp(Cetd[i][0]*dtv) - 1.0) * (Aetd[i][0] + Betd[i][0] ) / Cetd[i][0];
        v[i][1] = v[i][1]*exp(Cetd[i][1]*dtv) + (exp(Cetd[i][1]*dtv) - 1.0) * (Aetd[i][1] + Betd[i][1] ) / Cetd[i][1];
        v[i][2] = v[i][2]*exp(Cetd[i][2]*dtv) + (exp(Cetd[i][2]*dtv) - 1.0) * (Aetd[i][2] + Betd[i][2] ) / Cetd[i][2];
      }
  }
}








/* ---------------------------------------------------------------------- */

void FixETD::final_integrate()
{

  // update v of atoms in group
 
  double dtf = update->dt;
  double dtv = update->dt;
  double dtfm;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **Aetd = atom->Aetd;
  double **Betd = atom->Betd;
  double **Cetd = atom->Cetd;
  double *rmass = atom->rmass;
  double *mass = atom->mass;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];

       // x[i][0] += dtv * v[i][0];
       // x[i][1] += dtv * v[i][1];
       // x[i][2] += dtv * v[i][2];

	//v[i][0] += dtfm * f[i][0];
	//v[i][1] += dtfm * f[i][1];
	//v[i][2] += dtfm * f[i][2];

        //v[i][0] += v[i][0]*exp(Cetd[i][0]*dtfm) + (exp(Cetd[i][0]*dtfm) - 1.0) * (Aetd[i][0] + Betd[i][0]) / Cetd[i][0];   //dtfm * f[i][0];
        //v[i][1] += v[i][1]*exp(Cetd[i][1]*dtfm) + (exp(Cetd[i][1]*dtfm) - 1.0) * (Aetd[i][1] + Betd[i][1]) / Cetd[i][1];  //dtfm * f[i][1];
        //v[i][2] += v[i][2]*exp(Cetd[i][2]*dtfm) + (exp(Cetd[i][2]*dtfm) - 1.0) * (Aetd[i][2] + Betd[i][2]) / Cetd[i][2];  //dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];

       //v[i][0] = v[i][0]*exp(Cetd[i][0]*dtfm) + (exp(Cetd[i][0]*dtfm) - 1.0) * (Aetd[i][0] + Betd[i][0] ) / Cetd[i][0];
       //v[i][1] = v[i][1]*exp(Cetd[i][1]*dtfm) + (exp(Cetd[i][1]*dtfm) - 1.0) * (Aetd[i][1] + Betd[i][1] ) / Cetd[i][1];
       //v[i][2] = v[i][2]*exp(Cetd[i][2]*dtfm) + (exp(Cetd[i][2]*dtfm) - 1.0) * (Aetd[i][2] + Betd[i][2] ) / Cetd[i][2];

	//v[i][0] = v[i][0]*exp(Cetd[i][0]*dtv) + (exp(Cetd[i][0]*dtv) - 1.0) * (Aetd[i][0] + Betd[i][0] ) / Cetd[i][0];
       // v[i][1] = v[i][1]*exp(Cetd[i][1]*dtv) + (exp(Cetd[i][1]*dtv) - 1.0) * (Aetd[i][1] + Betd[i][1] ) / Cetd[i][1];
       // v[i][2] = v[i][2]*exp(Cetd[i][2]*dtv) + (exp(Cetd[i][2]*dtv) - 1.0) * (Aetd[i][2] + Betd[i][2] ) / Cetd[i][2];


// 	v[i][0] += dtv*(Cetd[i][0]*v[i][0] + Aetd[i][0] + Betd[i][0] );
//        v[i][1] += dtv*(Cetd[i][1]*v[i][1] + Aetd[i][1] + Betd[i][1] );
//        v[i][2] += dtv*(Cetd[i][2]*v[i][2] + Aetd[i][2] + Betd[i][2] );


//	v[i][0] += dtfm * f[i][0];
//	v[i][1] += dtfm * f[i][1];
//	v[i][2] += dtfm * f[i][2];

//        x[i][0] += dtv * v[i][0];
//        x[i][1] += dtv * v[i][1];
//        x[i][2] += dtv * v[i][2];



      }
  }
}



