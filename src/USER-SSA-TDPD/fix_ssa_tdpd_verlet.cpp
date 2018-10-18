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
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ssa_tdpd_verlet.h"
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

FixSsaTdpdVerlet::FixSsaTdpdVerlet(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"ssa_tdpd_verlet") != 0 && narg != 4)
    error->all(FLERR,"Illegal fix SsaTdpdVerlet command");

  time_integrate = 1;
  printf("caution:fix_ssa_tdpd_verlet, cpu is %d\n", comm->me);
  //MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

FixSsaTdpdVerlet::~FixSsaTdpdVerlet()
{
}

int FixSsaTdpdVerlet::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;


}

/* ---------------------------------------------------------------------- */

void FixSsaTdpdVerlet::init()
{
}

/* ---------------------------------------------------------------------- */

void FixSsaTdpdVerlet::initial_integrate(int vflag)
{
	double dtCm;

	double **C = atom->C;
	double **Q = atom->Q;



/*
	if (update->ntimestep == 2) {
	//printf(" i = %d, C[i][0] = %g \n", i, C[i][0] );
	//printf(" i = %d, Q[i][k] = %g \n", i, Q[i][k] );
	printf("Q[0][0] = %g \n",Q[0][0]);
	getchar();
         	                           }
*/

	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit){
			dtCm = 0.5*update->dt;

			for(int k = 0; k < atom->num_tdpd_species; k++)
			{	
		//		if(i==1) printf("i=%d  t=%d  T[i][10]=%f\n",i, update->ntimestep,T[i][10]);
				C[i][k] += Q[i][k] * dtCm ;
				C[i][k] = C[i][k] > 0 ? C[i][k] : 0.0;

			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixSsaTdpdVerlet::final_integrate()
{
  double dtCm;

  double **C = atom->C;
  double **Q = atom->Q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
		dtCm = 0.5*update->dt;                

		for(int k = 0; k < atom->num_tdpd_species; k++)
		{
			C[i][k] += dtCm * Q[i][k];
			C[i][k] = C[i][k] > 0 ? C[i][k] : 0.0;

	         }


	                       	 }
    }
/*
if (update->ntimestep == 2) {
	printf(" i = %d, C[i][0] = %g \n", i, C[i][0] );
	getchar();
         	                           }

  }
*/




}

/* ---------------------------------------------------------------------- */

