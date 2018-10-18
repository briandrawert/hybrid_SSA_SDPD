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
#include "fix_ssa_tsdpd_stationary.h"
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

FixSsaTsdpdStationary::FixSsaTsdpdStationary(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix ssa_tsdpd/stationary command requires atom_style with both energy and density, e.g. ssa_tsdpd");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix ssa_tsdpd/stationary command");

  time_integrate = 0;
  
  seed = comm->nprocs + comm->me + atom->nlocal;
  random = new RanMars (lmp, seed);
}

/* ---------------------------------------------------------------------- */

int FixSsaTsdpdStationary::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpdStationary::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixSsaTsdpdStationary::initial_integrate(int vflag) {

  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;
  
  double **C = atom->C;
  double **Q = atom->Q;

  int **Cd = atom->Cd;
  int **Qd = atom->Qd;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
  //    e[i] += dtf * de[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density
      for (int k = 0; k < atom->num_tdpd_species; k++){ // ...and concentrations
           C[i][k] += Q[i][k] *dtf;
           C[i][k] = C[i][k] > 0 ? C[i][k] : 0.0;
      }


    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpdStationary::final_integrate() {

//  printf("fix_ssa_tsdpd_stationary.cpp at t = %d", update->ntimestep);
  double *e = atom->e;
  double *de = atom->de;
  double *rho = atom->rho;
  double *drho = atom->drho;
 
  double **C = atom->C;
  double **Q = atom->Q;

  int **Cd = atom->Cd;
  int **Qd = atom->Qd;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nmax = atom->nmax;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
//      e[i] += dtf * de[i];
      rho[i] += dtf * drho[i];
      for (int k = 0; k < atom->num_tdpd_species; k++){
                C[i][k] += Q[i][k] *dtf;
                C[i][k] = C[i][k] > 0 ? C[i][k] : 0.0;
      }
    
      for (int s=0; s<atom->num_ssa_species; s++){
        Cd[i][s] += Qd[i][s];
        Cd[i][s] = Cd[i][s] > 0 ? Cd[i][s] : 0;
        //printf("Cd[%d][%d] = %d, Qd[%d][%d] = %d \n",i,s,Cd[i][s],i,s,Qd[i][s] );
      }

      // Calculate SSA reactions here
      double tt=0;
      double a0 = 0.0;
      int r,ro,k,s;
      for(r=0;r<atom->num_ssa_reactions;r++) a0 += atom->ssa_rxn_propensity[i][r];
      if(a0 > 0.0){
        double r1 = random->uniform();
        double r2 = random->uniform();
        tt += -log(1.0-r1)/a0;
        double old_a;
        double delta_a0=0.0;
        double a_sum = 0;

        while(tt < update->dt){ 
          delta_a0=0.0;
          //find next reaction to fire
          for(r=0;r<atom->num_ssa_reactions;r++){
            if((a_sum += atom->ssa_rxn_propensity[i][r]) > r2*a0) break;
          }
          // Change species populations for reaction r
          for(s=0;s<atom->num_ssa_species;s++){
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
  }
}

/* ---------------------------------------------------------------------- */

void FixSsaTsdpdStationary::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
