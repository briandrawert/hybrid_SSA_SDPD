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

#include <string.h>
#include <stdlib.h>
#include "atom_vec_ssa_tsdpd.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSsaTsdpd::AtomVecSsaTsdpd(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  comm_x_only = 0; // we communicate not only x forward but also vest ...
  comm_f_only = 0; // we also communicate de and drho in reverse direction
  size_forward = 12; // 3 + rho + e + vest[3], that means we may only communicate 5 in hybrid (was 8 before)
  size_reverse = 5; // 3 + drho + de
  size_border = 16; // 6 + rho + e + vest[3] + cv (was 12 before)
  size_velocity = 3;
  size_data_atom = 8;
  size_data_vel = 4;
  xcol_data = 6;

  atom->e_flag = 1;
  atom->rho_flag = 1;
  atom->cv_flag = 1;
  atom->vest_flag = 1;
  atom->tsdpd_flag = 1;


}



/* ----------------------------------------------------------------------
   process additional args
   set size_forward and size_border to max sizes
------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::process_args(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Invalid atom_style body command");

  atom->num_tdpd_species = atoi(arg[0]);
  atom->num_ssa_species = atoi(arg[1]);
  if(narg < 3){
    atom->num_ssa_reactions = 0;
  }else{
    atom->num_ssa_reactions = atoi(arg[2]);
  }

  atom->concentration_conversion = 1.0;
 
  if (narg == 4){  	 
    if (strcmp(arg[3],"concentration") == 0)  atom->Cd_concentration_flag = 1;
    else if (strcmp(arg[3],"population") == 0) atom->Cd_concentration_flag = 0;
   }

 if (narg == 5){
    if (strcmp(arg[3],"concentration") == 0)  atom->Cd_concentration_flag = 1;
    else if (strcmp(arg[3],"population") == 0) atom->Cd_concentration_flag = 0;
    atom->concentration_conversion = atof(arg[4]);
 }
  
  size_forward   += atom->num_tdpd_species + atom->num_ssa_species + 2 * atom->num_ssa_reactions * atom->num_ssa_species + atom->num_ssa_reactions;
  size_reverse   += atom->num_tdpd_species + atom->num_ssa_species;
  size_border    += atom->num_tdpd_species + atom->num_ssa_species + 2 * atom->num_ssa_reactions * atom->num_ssa_species + atom->num_ssa_reactions;
  size_data_atom += atom->num_tdpd_species + atom->num_ssa_species + 2 * atom->num_ssa_reactions * atom->num_ssa_species + atom->num_ssa_reactions;





}


/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
   ------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::grow(int n)
{

//  printf("in AtomVecSsaTsdpd::grow\n");

  int num_ssa_reactions = atom->num_ssa_reactions;
  int num_tdpd_species = atom->num_tdpd_species;
  int num_ssa_species = atom->num_ssa_species;
  

  if (num_ssa_species == 0) num_ssa_species = 1;
  if (num_ssa_reactions == 0) num_ssa_reactions = 1;
  
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag, nmax, "atom:tag");
  type = memory->grow(atom->type, nmax, "atom:type");
  mask = memory->grow(atom->mask, nmax, "atom:mask");
  image = memory->grow(atom->image, nmax, "atom:image");
  x = memory->grow(atom->x, nmax, 3, "atom:x");
  v = memory->grow(atom->v, nmax, 3, "atom:v");
  f = memory->grow(atom->f, nmax*comm->nthreads, 3, "atom:f");

  rho = memory->grow(atom->rho, nmax, "atom:rho");
  drho = memory->grow(atom->drho, nmax*comm->nthreads, "atom:drho");
  e = memory->grow(atom->e, nmax, "atom:e");
  de = memory->grow(atom->de, nmax*comm->nthreads, "atom:de");
  vest = memory->grow(atom->vest, nmax, 3, "atom:vest");
  cv = memory->grow(atom->cv, nmax, "atom:cv");
 
  C = memory->grow(atom->C,nmax,atom->num_tdpd_species,"atom:C"); //added (grow C)
  Q = memory->grow(atom->Q,nmax*comm->nthreads,atom->num_tdpd_species,"atom:Q"); //added (grow Q)

  Cd = memory->grow(atom->Cd,nmax,num_ssa_species,"atom:Cd"); //added (grow Cd)
  Qd = memory->grow(atom->Qd,nmax*comm->nthreads,num_ssa_species,"atom:Qd"); //added (grow Qd)
  ssa_rxn_propensity = memory->grow(atom->ssa_rxn_propensity,nmax*comm->nthreads,num_ssa_reactions,"atom:ssa_rxn_propensity"); //added (grow ssa_rxn_propensity)
  d_ssa_rxn_prop_d_c = memory->grow(atom->d_ssa_rxn_prop_d_c,nmax*comm->nthreads,num_ssa_reactions,num_ssa_species,"atom:d_ssa_rxn_prop_d_c"); //added (grow d_ssa_rxn_prop_d_c)
  ssa_stoich_matrix = memory->grow(atom->ssa_stoich_matrix,nmax*comm->nthreads,num_ssa_reactions,num_ssa_species,"atom:ssa_stoich_matrix"); //added (grow ssa_stoich_matrix)
  dfsp_D_matrix = memory->grow(atom->dfsp_D_matrix,nmax*nmax,"atom:dfsp_D_matrix");
  dfsp_D_diag = memory->grow(atom->dfsp_D_diag,nmax,"atom:dfsp_D_diag");
  dfsp_Diffusion_coeff = memory->grow(atom->dfsp_Diffusion_coeff,nmax*nmax*num_ssa_species,"atom:dfsp_Diffusion_coeff");
  dfsp_a_i = memory->grow(atom->dfsp_a_i, nmax,"atom:dfsp_a_i");


  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
   ------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::grow_reset() {

//  printf("in AtomVecSsaTsdpd::reset\n");

  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;
  rho = atom->rho;
  drho = atom->drho;
  e = atom->e;
  de = atom->de;
  vest = atom->vest;
  cv = atom->cv;
  C = atom->C; Q = atom->Q; //added 


  Cd = atom->Cd; Qd = atom->Qd; //added
  ssa_rxn_propensity = atom->ssa_rxn_propensity; //added 

  d_ssa_rxn_prop_d_c = atom->d_ssa_rxn_prop_d_c; //added

  ssa_stoich_matrix = atom->ssa_stoich_matrix;  //added
  dfsp_D_matrix = atom->dfsp_D_matrix; //added
  dfsp_D_diag = atom->dfsp_D_diag; //added
  dfsp_Diffusion_coeff = atom->dfsp_Diffusion_coeff; //added
  dfsp_a_i = atom->dfsp_a_i; //added


}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::copy(int i, int j, int delflag) {

//  printf("in AtomVecSsaTsdpd::copy\n");


  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  rho[j] = rho[i];
  drho[j] = drho[i];
  e[j] = e[i];
  de[j] = de[i];
  cv[j] = cv[i];
  vest[j][0] = vest[i][0];
  vest[j][1] = vest[i][1];
  vest[j][2] = vest[i][2];
  for (int k = 0; k < atom->num_tdpd_species; k++)  C[j][k] = C[i][k]; //added

  for (int k = 0; k < atom->num_ssa_species; k++)  Cd[j][k] = Cd[i][k]; //added

  for (int r = 0; r < atom->num_ssa_reactions; r++)  ssa_rxn_propensity[j][r] = ssa_rxn_propensity[i][r]; //added

  for (int r = 0; r < atom->num_ssa_reactions; r++)
    for (int k = 0; k < atom->num_ssa_species; k++)
      d_ssa_rxn_prop_d_c[j][r][k] = d_ssa_rxn_prop_d_c[i][r][k]; //added

  for (int r = 0; r < atom->num_ssa_reactions; r++)
    for (int k = 0; k < atom->num_ssa_species; k++)
      ssa_stoich_matrix[j][r][k] = ssa_stoich_matrix[i][r][k]; //added

  dfsp_D_matrix[j] = dfsp_D_matrix[i]; //added
  dfsp_D_diag[j] = dfsp_D_diag[i]; //added
  dfsp_Diffusion_coeff[j] = dfsp_Diffusion_coeff[i]; //added
  dfsp_a_i[j] = dfsp_a_i[i]; //added

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i, j,delflag);

}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::force_clear(int n, size_t nbytes)
{

//  printf("in AtomVecSsaTsdpd::force_clear\n");

  memset(&de[n],0,nbytes);
  memset(&drho[n],0,nbytes);
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_comm_hybrid(int n, int *list, double *buf) {

//  printf("in AtomVecSsaTsdpd::pack_comm_hybrid\n");

  int i, j, m, r, k;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
    buf[m++] = e[j];
    buf[m++] = vest[j][0];
    buf[m++] = vest[j][1];
    buf[m++] = vest[j][2];
    for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = C[j][k];

    for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Cd[j][k];

    for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

   buf[m++] = dfsp_D_matrix[j];
   buf[m++] = dfsp_D_diag[j];
   buf[m++] = dfsp_Diffusion_coeff[j];
   buf[m++] = dfsp_a_i[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::unpack_comm_hybrid(int n, int first, double *buf) {

//  printf("in AtomVecSsaTsdpd::unpack_comm_hybrid\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    rho[i] = buf[m++];
    e[i] = buf[m++];
    vest[i][0] = buf[m++];
    vest[i][1] = buf[m++];
    vest[i][2] = buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++)  C[i][k] = buf[m++];
  
    for (int k = 0; k < atom->num_ssa_species; k++)  Cd[i][k] = (int) buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[i][r] = buf[m++];


    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            d_ssa_rxn_prop_d_c[i][r][k] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            ssa_stoich_matrix[i][r][k] = buf[m++];


   dfsp_D_matrix[i]=buf[m++];
   dfsp_D_diag[i]=buf[m++];
   dfsp_Diffusion_coeff[i]=buf[m++];
   dfsp_a_i[i]=buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_border_hybrid(int n, int *list, double *buf) {
  
//  printf("in AtomVecSsaTsdpd::pack_border_hybrid\n");

  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
    buf[m++] = e[j];
    buf[m++] = vest[j][0];
    buf[m++] = vest[j][1];
    buf[m++] = vest[j][2];
    for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = C[j][k];


    for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Cd[j][k];

    for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];


    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

     for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];
 
  buf[m++] = dfsp_D_matrix[j];
  buf[m++] = dfsp_D_diag[j];
  buf[m++] = dfsp_Diffusion_coeff[j];
  buf[m++] = dfsp_a_i[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::unpack_border_hybrid(int n, int first, double *buf) {

//  printf("in AtomVecSsaTsdpd::unpack_border_hybrid\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    rho[i] = buf[m++];
    e[i] = buf[m++];
    vest[i][0] = buf[m++];
    vest[i][1] = buf[m++];
    vest[i][2] = buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++)  C[i][k] = buf[m++];


    for (int k = 0; k < atom->num_ssa_species; k++)  Cd[i][k] = (int) buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[i][r] = buf[m++];


    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            d_ssa_rxn_prop_d_c[i][r][k] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            ssa_stoich_matrix[i][r][k] = buf[m++];

  dfsp_D_matrix[i] = buf[m++];
  dfsp_D_diag[i] = buf[m++];
  dfsp_Diffusion_coeff[i] = buf[m++];
  dfsp_a_i[i] = buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_reverse_hybrid(int n, int first, double *buf) {

//  printf("in AtomVecSsaTsdpd::pack_reverse_hybrid\n");

  int i,r, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = drho[i];
    buf[m++] = de[i];
    for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = Q[i][k];
    for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Qd[i][k];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::unpack_reverse_hybrid(int n, int *list, double *buf) {

//  printf("in AtomVecSsaTsdpd::unpack_reverse_hybrid\n");

  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    drho[j] += buf[m++];
    de[j] += buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++)  C[j][k] += buf[m++];

    for (int k = 0; k < atom->num_ssa_species; k++)  Cd[j][k] += (int) buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[j][r] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            d_ssa_rxn_prop_d_c[j][r][k] = buf[m++];


    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            ssa_stoich_matrix[j][r][k] = buf[m++];

  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_comm(int n, int *list, double *buf, int pbc_flag,
                           int *pbc) {

//  printf("in AtomVecSsaTsdpd::pack_comm\n");

  int i, j, m, r, k;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

     buf[m++] = dfsp_D_matrix[j];
     buf[m++] = dfsp_D_diag[j];
     buf[m++] = dfsp_Diffusion_coeff[j];
     buf[m++] = dfsp_a_i[j];

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

      buf[m++] = dfsp_D_matrix[j];
      buf[m++] = dfsp_D_diag[j];
      buf[m++] = dfsp_Diffusion_coeff[j];
      buf[m++] = dfsp_a_i[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_comm_vel(int n, int *list, double *buf, int pbc_flag,
                               int *pbc) {

//  printf("in AtomVecSsaTsdpd::pack_comm_vel\n");

  int i, j, m, r, k;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

      buf[m++] = dfsp_D_matrix[j];
      buf[m++] = dfsp_D_diag[j];
      buf[m++] = dfsp_Diffusion_coeff[j];
      buf[m++] = dfsp_a_i[j];
      
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

     buf[m++] = dfsp_D_matrix[j];
     buf[m++] = dfsp_D_diag[j];
     buf[m++] = dfsp_Diffusion_coeff[j];
     buf[m++] = dfsp_a_i[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::unpack_comm(int n, int first, double *buf) {

//  printf("in AtomVecSsaTsdpd::unpack_comm\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    rho[i] = buf[m++];
    e[i] = buf[m++];
    vest[i][0] = buf[m++];
    vest[i][1] = buf[m++];
    vest[i][2] = buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++) C[i][k] = buf[m++];

    for (int k = 0; k < atom->num_ssa_species; k++) Cd[i][k] = (int) buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[i][r] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            d_ssa_rxn_prop_d_c[i][r][k] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            ssa_stoich_matrix[i][r][k] = buf[m++];

   dfsp_D_matrix[i] = buf[m++];
   dfsp_D_diag[i] = buf[m++];
   dfsp_Diffusion_coeff[i] = buf[m++];
   dfsp_a_i[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::unpack_comm_vel(int n, int first, double *buf) {

//  printf("in AtomVecSsaTsdpd::unpack_comm_vel\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    rho[i] = buf[m++];
    e[i] = buf[m++];
    vest[i][0] = buf[m++];
    vest[i][1] = buf[m++];
    vest[i][2] = buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++) C[i][k] = buf[m++];

    for (int k = 0; k < atom->num_ssa_species; k++) Cd[i][k] = (int) buf[m++];

//    for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[i][r]; modified (backwards!)
    for (int r = 0; r < atom->num_ssa_reactions; r++)  ssa_rxn_propensity[i][r] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
		d_ssa_rxn_prop_d_c[i][r][k] = buf[m++];
       //     buf[m++] = d_ssa_rxn_prop_d_c[i][r][k];  modified (backwards!)

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
		ssa_stoich_matrix[i][r][k] = buf[m++];
         //   buf[m++] = ssa_stoich_matrix[i][r][k];   modified (backwards!)
  
   dfsp_D_matrix[i] = buf[m++];
   dfsp_D_diag[i] = buf[m++];
   dfsp_Diffusion_coeff[i] = buf[m++];
   dfsp_a_i[i] = buf[m++];

  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_reverse(int n, int first, double *buf) {

//  printf("in AtomVecSsaTsdpd::pack_reverse\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = drho[i];
    buf[m++] = de[i];
    for (int k = 0; k < atom->num_tdpd_species; k++)  buf[m++] = Q[i][k];  //added (this packs source term for tdpd model)
    for (int k = 0; k < atom->num_ssa_species; k++)  buf[m++] = (double) Qd[i][k];  //added (this packs source term for tdpd model)
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::unpack_reverse(int n, int *list, double *buf) {

//  printf("in AtomVecSsaTsdpd::unpack_reverse\n");

  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    drho[j] += buf[m++];
    de[j] += buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++) Q[j][k] += buf[m++]; //added (this unpacks source term for tdpd model)
    for (int k = 0; k < atom->num_ssa_species; k++) Qd[j][k] += (int) buf[m++]; //added (this unpacks source term for tdpd model)

  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_border(int n, int *list, double *buf, int pbc_flag,
                             int *pbc) {

//  printf("in AtomVecSsaTsdpd::pack_border\n");

  int i, j, m, r, k;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = cv[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

      buf[m++] = dfsp_D_matrix[j];
      buf[m++] = dfsp_D_diag[j];
      buf[m++] = dfsp_Diffusion_coeff[j];
      buf[m++] = dfsp_a_i[j];

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = cv[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

     for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];
     
     buf[m++] = dfsp_D_matrix[j]; 
     buf[m++] = dfsp_D_diag[j]; 
     buf[m++] = dfsp_Diffusion_coeff[j]; 
     buf[m++] = dfsp_a_i[j]; 
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_border_vel(int n, int *list, double *buf, int pbc_flag,
                                 int *pbc)
{
  int i,j,m, r,k;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = rho[j];
      buf[m++] = e[j];
      buf[m++] = cv[j];
      buf[m++] = vest[j][0];
      buf[m++] = vest[j][1];
      buf[m++] = vest[j][2];
      for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[j][k];

      for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[j][k];

      for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

      for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

     for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

     buf[m++] = dfsp_D_matrix[j];
     buf[m++] = dfsp_D_diag[j];
     buf[m++] = dfsp_Diffusion_coeff[j];
     buf[m++] = dfsp_a_i[j];

    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = rho[j];
        buf[m++] = e[j];
        buf[m++] = cv[j];
        buf[m++] = vest[j][0];
        buf[m++] = vest[j][1];
        buf[m++] = vest[j][2];
        for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[j][k];

        for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[j][k];

        for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

        for (int r = 0; r < atom->num_ssa_reactions; r++)
          for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

        for (int r = 0; r < atom->num_ssa_reactions; r++)
          for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];
        
        buf[m++] = dfsp_D_matrix[j];
        buf[m++] = dfsp_D_diag[j];
        buf[m++] = dfsp_Diffusion_coeff[j];
        buf[m++] = dfsp_a_i[j];
      
      }
    } else {
      dvx = pbc[0] * h_rate[0] + pbc[5] * h_rate[5] + pbc[4] * h_rate[4];
      dvy = pbc[1] * h_rate[1] + pbc[3] * h_rate[3];
      dvz = pbc[2] * h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
          buf[m++] = vest[j][0] + dvx;
          buf[m++] = vest[j][1] + dvy;
          buf[m++] = vest[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = vest[j][0];
          buf[m++] = vest[j][1];
          buf[m++] = vest[j][2];
        }
        buf[m++] = rho[j];
        buf[m++] = e[j];
        buf[m++] = cv[j];
        for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[j][k];

        for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[j][k];

        for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[j][r];

        for (int r = 0; r < atom->num_ssa_reactions; r++)
          for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = d_ssa_rxn_prop_d_c[j][r][k];

        for (int r = 0; r < atom->num_ssa_reactions; r++)
          for (int k = 0; k < atom->num_ssa_species; k++)
            buf[m++] = ssa_stoich_matrix[j][r][k];

        buf[m++] = dfsp_D_matrix[j];
        buf[m++] = dfsp_D_diag[j];
        buf[m++] = dfsp_Diffusion_coeff[j];
        buf[m++] = dfsp_a_i[j];

      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::unpack_border(int n, int first, double *buf) {
  
//  printf("in AtomVecSsaTsdpd::unpack_border\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax)
      grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    rho[i] = buf[m++];
    e[i] = buf[m++];
    cv[i] = buf[m++];
    vest[i][0] = buf[m++];
    vest[i][1] = buf[m++];
    vest[i][2] = buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++) C[i][k] = buf[m++];

    for (int k = 0; k < atom->num_ssa_species; k++) Cd[i][k] = (int) buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[i][r] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            d_ssa_rxn_prop_d_c[i][r][k] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            ssa_stoich_matrix[i][r][k] = buf[m++];

    dfsp_D_matrix[i] = buf[m++];
    dfsp_D_diag[i] = buf[m++];
    dfsp_Diffusion_coeff[i] = buf[m++];
    dfsp_a_i[i] = buf[m++];

  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecSsaTsdpd::unpack_border_vel(int n, int first, double *buf) {

  //printf("in AtomVecSsaTsdpd::unpack_border_vel\n");

  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax)
      grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    vest[i][0] = buf[m++];
    vest[i][1] = buf[m++];
    vest[i][2] = buf[m++];
    rho[i] = buf[m++];
    e[i] = buf[m++];
    cv[i] = buf[m++];
    for (int k = 0; k < atom->num_tdpd_species; k++) C[i][k] = buf[m++];

    for (int k = 0; k < atom->num_ssa_species; k++) Cd[i][k] = (int) buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[i][r] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            d_ssa_rxn_prop_d_c[i][r][k] = buf[m++];

    for (int r = 0; r < atom->num_ssa_reactions; r++)
        for (int k = 0; k < atom->num_ssa_species; k++)
            ssa_stoich_matrix[i][r][k] = buf[m++];

    dfsp_D_matrix[i] = buf[m++];
    dfsp_D_diag[i] = buf[m++];
    dfsp_Diffusion_coeff[i] = buf[m++];
    dfsp_a_i[i] = buf[m++];

  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
   ------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_exchange(int i, double *buf) {
  
//  printf("in AtomVecSsaTsdpd::pack_exchange\n");

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = rho[i];
  buf[m++] = e[i];
  buf[m++] = cv[i];
  buf[m++] = vest[i][0];
  buf[m++] = vest[i][1];
  buf[m++] = vest[i][2];
  for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[i][k];

  for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[i][k];

  for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[i][r];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
    for (int k = 0; k < atom->num_ssa_species; k++)
        buf[m++] = d_ssa_rxn_prop_d_c[i][r][k];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
    for (int k = 0; k < atom->num_ssa_species; k++)
        buf[m++] = ssa_stoich_matrix[i][r][k];

  buf[m++] = dfsp_D_matrix[i];
  buf[m++] = dfsp_D_diag[i];
  buf[m++] = dfsp_Diffusion_coeff[i];
  buf[m++] = dfsp_a_i[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i, &buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSsaTsdpd::unpack_exchange(double *buf) {

  //printf("in AtomVecSsaTsdpd::unpack_exchange\n");

  int nlocal = atom->nlocal;
  if (nlocal == nmax)
    grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  rho[nlocal] = buf[m++];
  e[nlocal] = buf[m++];
  cv[nlocal] = buf[m++];
  vest[nlocal][0] = buf[m++];
  vest[nlocal][1] = buf[m++];
  vest[nlocal][2] = buf[m++];
  for (int k = 0; k < atom->num_tdpd_species; k++) C[nlocal][k] = buf[m++];

  for (int k = 0; k < atom->num_ssa_species; k++) Cd[nlocal][k] = (int) buf[m++];

  for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[nlocal][r] = buf[m++];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
         d_ssa_rxn_prop_d_c[nlocal][r][k] = buf[m++];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
         ssa_stoich_matrix[nlocal][r][k] = buf[m++];

  dfsp_D_matrix[nlocal] = buf[m++];
  dfsp_D_diag[nlocal] = buf[m++];
  dfsp_Diffusion_coeff[nlocal] = buf[m++];
  dfsp_a_i[nlocal] = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]-> unpack_exchange(nlocal,
                                                                   &buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
   ------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::size_restart() {
  int i;

  int nlocal = atom->nlocal;
  int n = ( 17 +  atom->num_tdpd_species + atom->num_ssa_species + 2 * atom->num_ssa_reactions * atom->num_ssa_species + atom->num_ssa_reactions) * nlocal; // 11 + rho + e + cv + vest[3]
  
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
   ------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_restart(int i, double *buf) {
  int m = 1;
  int k,r;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = rho[i];
  buf[m++] = e[i];
  buf[m++] = cv[i];
  buf[m++] = vest[i][0];
  buf[m++] = vest[i][1];
  buf[m++] = vest[i][2];
  for (int k = 0; k < atom->num_tdpd_species; k++) buf[m++] = C[i][k];

  for (int k = 0; k < atom->num_ssa_species; k++) buf[m++] = (double) Cd[i][k];

  for (int r = 0; r < atom->num_ssa_reactions; r++)  buf[m++] = ssa_rxn_propensity[i][r];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
    for (int k = 0; k < atom->num_ssa_species; k++)
        buf[m++] = d_ssa_rxn_prop_d_c[i][r][k];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
     for (int k = 0; k < atom->num_ssa_species; k++)
        buf[m++] = ssa_stoich_matrix[i][r][k];


  buf[m++] = dfsp_D_matrix[i];
  buf[m++] = dfsp_D_diag[i];
  buf[m++] = dfsp_Diffusion_coeff[i];
  buf[m++] = dfsp_a_i[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i, &buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
   ------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::unpack_restart(double *buf) {
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra, nmax, atom->nextra_store, "atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  rho[nlocal] = buf[m++];
  e[nlocal] = buf[m++];
  cv[nlocal] = buf[m++];
  vest[nlocal][0] = buf[m++];
  vest[nlocal][1] = buf[m++];
  vest[nlocal][2] = buf[m++];
  for(int k = 0; k < atom->num_tdpd_species; k++) C[nlocal][k] = buf[m++]; //added

  for(int k = 0; k < atom->num_ssa_species; k++) Cd[nlocal][k] = buf[m++]; //added

  for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[nlocal][r] = buf[m++];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
         d_ssa_rxn_prop_d_c[nlocal][r][k] = buf[m++];

  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
         ssa_stoich_matrix[nlocal][r][k] = buf[m++];

  dfsp_D_matrix[nlocal] = buf[m++];
  dfsp_D_diag[nlocal] = buf[m++];
  dfsp_Diffusion_coeff[nlocal] = buf[m++];
  dfsp_a_i[nlocal] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++)
      extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
   ------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::create_atom(int itype, double *coord) {
  int nlocal = atom->nlocal;
  if (nlocal == nmax)
    grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  rho[nlocal] = 0.0;
  e[nlocal] = 0.0;
  cv[nlocal] = 1.0;
  vest[nlocal][0] = 0.0;
  vest[nlocal][1] = 0.0;
  vest[nlocal][2] = 0.0;
  de[nlocal] = 0.0;
  drho[nlocal] = 0.0;
  for (int k = 0; k < atom->num_tdpd_species; k++) C[nlocal][k] = 0.0;

  for (int k = 0; k < atom->num_ssa_species; k++) Cd[nlocal][k] = 0;

  for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[nlocal][r] = 0.0;

  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
         d_ssa_rxn_prop_d_c[nlocal][r][k] = 0.0;

  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
         ssa_stoich_matrix[nlocal][r][k] = 0.0;


  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
   ------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::data_atom(double *coord, imageint imagetmp, char **values) {
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  rho[nlocal] = atof(values[2]);
  e[nlocal] = atof(values[3]);
  cv[nlocal] = atof(values[4]);

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  

  int m = 5;

  for (int k = 0; k < atom->num_tdpd_species; k++) C[nlocal][k] = atof(values[m+k]);
  m += atom->num_tdpd_species;

  for (int k = 0; k < atom->num_ssa_species; k++) Cd[nlocal][k] = atof(values[m+k]);
  m += atom->num_ssa_species;

  for (int r = 0; r < atom->num_ssa_reactions; r++) ssa_rxn_propensity[nlocal][r] = atof(values[m+r]);
  m += atom->num_ssa_reactions;


  for (int r = 0; r < atom->num_ssa_reactions; r++)
      for (int k = 0; k < atom->num_ssa_species; k++)
        d_ssa_rxn_prop_d_c[nlocal][r][k] = atof(values[m+r*atom->num_ssa_reactions+k]);
  m += atom->num_ssa_reactions*atom->num_ssa_species;



  //printf("rho=%f, e=%f, cv=%f, x=%f\n", rho[nlocal], e[nlocal], cv[nlocal], x[nlocal][0]);

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  vest[nlocal][0] = 0.0;
  vest[nlocal][1] = 0.0;
  vest[nlocal][2] = 0.0;

  de[nlocal] = 0.0;
  drho[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
   ------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::data_atom_hybrid(int nlocal, char **values) {

  rho[nlocal] = atof(values[0]);
  e[nlocal] = atof(values[1]);
  cv[nlocal] = atof(values[2]);

  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = rho[i];
    buf[i][3] = e[i];
    buf[i][4] = cv[i];
    buf[i][5] = x[i][0];
    buf[i][6] = x[i][1];
    buf[i][7] = x[i][2];
    buf[i][8] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][10] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::pack_data_hybrid(int i, double *buf)
{
  buf[0] = rho[i];
  buf[1] = e[i];
  buf[2] = cv[i];
  return 3;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e "
            "%d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],
            buf[i][5],buf[i][6],buf[i][7],
            (int) ubuf(buf[i][8]).i,(int) ubuf(buf[i][9]).i,
            (int) ubuf(buf[i][10]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecSsaTsdpd::property_atom(char *name)
{
  if (strcmp(name,"rho") == 0) return 0;
  if (strcmp(name,"drho") == 0) return 1;
  if (strcmp(name,"e") == 0) return 2;
  if (strcmp(name,"de") == 0) return 3;
  if (strcmp(name,"cv") == 0) return 4;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecSsaTsdpd::pack_property_atom(int index, double *buf,
                                     int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = rho[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = drho[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = e[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = de[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = cv[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
   ------------------------------------------------------------------------- */

bigint AtomVecSsaTsdpd::memory_usage() {
  bigint bytes = 0;

//  printf("in AtomVecSsaTsdpd::memory_usage\n");

  if (atom->memcheck("tag"))
    bytes += memory->usage(tag, nmax);
  if (atom->memcheck("type"))
    bytes += memory->usage(type, nmax);
  if (atom->memcheck("mask"))
    bytes += memory->usage(mask, nmax);
  if (atom->memcheck("image"))
    bytes += memory->usage(image, nmax);
  if (atom->memcheck("x"))
    bytes += memory->usage(x, nmax, 3);
  if (atom->memcheck("v"))
    bytes += memory->usage(v, nmax, 3);
  if (atom->memcheck("f"))
    bytes += memory->usage(f, nmax*comm->nthreads, 3);
  if (atom->memcheck("rho"))
    bytes += memory->usage(rho, nmax);
  if (atom->memcheck("drho"))
    bytes += memory->usage(drho, nmax*comm->nthreads);
  if (atom->memcheck("e"))
    bytes += memory->usage(e, nmax);
  if (atom->memcheck("de"))
    bytes += memory->usage(de, nmax*comm->nthreads);
  if (atom->memcheck("cv"))
    bytes += memory->usage(cv, nmax);
  if (atom->memcheck("vest"))
    bytes += memory->usage(vest, nmax);
  if (atom->memcheck("C")) 
    bytes += memory->usage(C,nmax,atom->num_tdpd_species); //added
  if (atom->memcheck("Q")) 
    bytes += memory->usage(Q,nmax*comm->nthreads,atom->num_tdpd_species); //added

  if (atom->memcheck("Cd")) 
    bytes += memory->usage(Cd,nmax,atom->num_ssa_species); //added
  if (atom->memcheck("Qd")) 
    bytes += memory->usage(Qd,nmax*comm->nthreads,atom->num_ssa_species); //added

  if (atom->memcheck("ssa_rxn_propensity")) 
    bytes += memory->usage(ssa_rxn_propensity,nmax*comm->nthreads,atom->num_ssa_reactions); //added

  if (atom->memcheck("d_ssa_rxn_prop_d_c")) 
    bytes += memory->usage(d_ssa_rxn_prop_d_c,nmax*comm->nthreads,atom->num_ssa_reactions,atom->num_ssa_species); //added

  if (atom->memcheck("ssa_stoich_matrix")) 
    bytes += memory->usage(ssa_stoich_matrix,nmax*comm->nthreads,atom->num_ssa_reactions,atom->num_ssa_species); //added

  if (atom->memcheck("dfsp_D_matrix"))
    bytes += memory->usage(dfsp_D_matrix, nmax); //added

   if (atom->memcheck("dfsp_D_diag"))
    bytes += memory->usage(dfsp_D_diag, nmax); //added

  if (atom->memcheck("dfsp_Diffusion_coeff"))
    bytes += memory->usage(dfsp_Diffusion_coeff, nmax*nmax*atom->num_ssa_species); //added

  if (atom->memcheck("dfsp_a_i"))
    bytes += memory->usage(dfsp_a_i, nmax); //added

  return bytes;

}
