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
#include <string.h>
#include "fix_ssa_tsdpd_verlet.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSsaTsdpd::FixSsaTsdpd(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix ssa_tsdpd/verlet command requires atom_style with both energy and density");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix ssa_tsdpd/verlet command");

  time_integrate = 1;

  // seed is immune to underflow/overflow because it is unsigned
  seed = comm->nprocs + comm->me + atom->nlocal;
  //if (narg == 3) seed += force->inumeric (FLERR, arg[2]);
  random = new RanMars (lmp, seed);

}

/* ---------------------------------------------------------------------- */

int FixSsaTsdpd::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpd::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

void FixSsaTsdpd::setup_pre_force(int vflag)
{
  // set vest equal to v
  double **v = atom->v;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixSsaTsdpd::initial_integrate(int vflag) {
  // update v and x and rho and e of atoms in group
 
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  double dtCm;
  double **C = atom->C;
  double **Q = atom->Q;
  int **Cd = atom->Cd;
  int **Qd = atom->Qd;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }


        //extrapolate velocity for use with velocity-dependent potentials, e.g. SPH
      vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
      vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
      vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

	
      dtCm = 0.5*update->dt;      
      for (int k = 0; k < atom->num_tdpd_species; k++){
		C[i][k] += Q[i][k] *dtf;
		C[i][k] = C[i][k] > 0 ? C[i][k] : 0.0;
      }

      //e[i] += dtf * de[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density
 



    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpd::final_integrate() {

  // update v, rho, and e of atoms in group
  // printf("pair_ssa_tsdpd_verlet.cpp final integrate at t = %d \n", update->ntimestep);

  double **v = atom->v;
  double **f = atom->f;
  double *e = atom->e;
  double *de = atom->de;
  double *rho = atom->rho;
  double *drho = atom->drho;

  double dtCm;
  double **C = atom->C;
  double **Q = atom->Q;
  int **Cd = atom->Cd;
  int **Qd = atom->Qd;


  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;
  int r,ro,k;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];


      dtCm = 0.5*update->dt;      
      for (k = 0; k < atom->num_tdpd_species; k++){
		C[i][k] += Q[i][k] *dtf;
		C[i][k] = C[i][k] > 0 ? C[i][k] : 0.0;  // enforce positivity, but this method can lead to instabilities
      }


      // TODO: Convert Qd to Cd flux here
      for (int s=0; s<atom->num_ssa_species; s++){
          Cd[i][s] += Qd[i][s];
	  Cd[i][s] = Cd[i][s] > 0 ? Cd[i][s] : 0;
          //printf("Cd[%d][%d] = %d, Qd[%d][%d] = %d \n",i,s,Cd[i][s],i,s,Qd[i][s] );
      }


      // Calculate SSA reactions here
      
       if (atom->num_ssa_species > 0) {
        double tt=0;
        double a0 = 0.0;
        for(r=0;r<atom->num_ssa_reactions;r++) a0 += atom->ssa_rxn_propensity[i][r];
        if(a0 > 0.0){
            double r1 = random->uniform();
            double r2 = random->uniform();
            tt += -log(1.0-r1)/a0;
            double old_a;
            double delta_a0=0.0;
            double a_sum = 0;

            while(tt < update->dt){
              //find next reaction to fire
              for(r=0;r<atom->num_ssa_reactions;r++){
                  if((a_sum += atom->ssa_rxn_propensity[i][r]) > r2*a0) break;
              }
              // Change species populations for reaction r
              for(int s=0;s<atom->num_ssa_species;s++){
                Cd[i][s] += atom->ssa_stoich_matrix[i][r][s];
                // Change reaction propensities
                for(ro=0;ro<atom->num_ssa_reactions;ro++){
                  if(atom->d_ssa_rxn_prop_d_c[i][ro][s] != 0){
                    old_a = atom->ssa_rxn_propensity[i][ro];
                    atom->ssa_rxn_propensity[i][ro] += atom->d_ssa_rxn_prop_d_c[i][ro][s];
                    delta_a0 += (atom->ssa_rxn_propensity[i][ro] - old_a);
                  }
                }
              }
              // Update total propensity
              a0 += delta_a0;
              // Roll new random numbers
              r1 = random->uniform();
              r2 = random->uniform();
              // Calcuate time to next reaction
              tt += -log(1.0-r1)/a0;
          }
        }
      }
      //e[i] += dtf * de[i];
      rho[i] += dtf * drho[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpd::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
