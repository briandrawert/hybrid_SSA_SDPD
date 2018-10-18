#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_ssa_tsdpd_ssa_rxn_mass_action.h"
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
// Example reaction of form "A+B -k-> C"
//#                                            rxn_index rxn_rate #reactants reactant1 reactant2 #products product1
//fix rxn1 all ssa_tsdpd/chem_rxn_mass_action  0          0.1      2          0         1         1         2
/* ---------------------------------------------------------------------- */

FixSsaTsdpdSsaRxnMassAction::FixSsaTsdpdSsaRxnMassAction(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"ssa_tsdpd_ssa_rxn_mass_action") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix ssa_tsdpd_ssa_rxn_mass_action command, first error.");

  int i;
  int arg_index = 3;
  
  rxn_index = atof(arg[arg_index++]);
  
  k_rate = atof(arg[arg_index++]);
  num_reactants = atoi(arg[arg_index++]);
  if(num_reactants > atom->num_tdpd_species) 
        error->all(FLERR,"Illegal fix ssa_tsdpd_ssa_rxn_mass_action command -- number of reactant species greater than number of species.\n");
  if(num_reactants > 2) 
        error->all(FLERR,"Illegal fix ssa_tsdpd_ssa_rxn_mass_action command -- mass action reactions can have at most 2 reactants.\n");
  if(num_reactants > 0){
    for(i=0;i<num_reactants;i++){
      reactants[i] = atoi(arg[arg_index++]);
    }
  }
  num_products = atoi(arg[arg_index++]);
  if(num_products > atom->num_tdpd_species) 
        error->all(FLERR,"Illegal fix ssa_tsdpd_ssa_rxn_mass_action command -- number of product species greater than number of species.\n");
  if(num_products > 4) 
        error->all(FLERR,"Illegal fix ssa_tsdpd_ssa_rxn_mass_action command -- maximum number of product limited to 4 .\n");
  if(num_products > 0){
    for(i=0;i<num_products;i++){
      products[i] = atoi(arg[arg_index++]);
    }
  }

  MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */
FixSsaTsdpdSsaRxnMassAction::~FixSsaTsdpdSsaRxnMassAction()
{
}

/* ---------------------------------------------------------------------- */
int FixSsaTsdpdSsaRxnMassAction::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */
void FixSsaTsdpdSsaRxnMassAction::init()
{
}

/* ---------------------------------------------------------------------- */
void FixSsaTsdpdSsaRxnMassAction::post_force(int vflag)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int **C = atom->Cd;
  int **Q = atom->Qd;
  
  double *rho = atom->rho;
  double *mass = atom->mass;
  int *type = atom->type;

  double **ssa_rxn_propensity = atom->ssa_rxn_propensity;
  double ***d_ssa_rxn_prop_d_c = atom->d_ssa_rxn_prop_d_c;
  int ***ssa_stoich_matrix = atom->ssa_stoich_matrix;
  
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;


  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
        //volume = 1.0; // TODO: find the volume of each voxel
        itype = type[i];
	volume = mass[itype] / rho[i];

      if(num_reactants==2){
        if(reactants[0] == reactants[1]){
            ssa_rxn_propensity[i][rxn_index] = k_rate/volume/2.0*C[i][reactants[0]]*(C[i][reactants[0]] - 1);
            d_ssa_rxn_prop_d_c[i][rxn_index][reactants[0]] = k_rate/volume;
            ssa_stoich_matrix[i][rxn_index][reactants[0]] = -2;
        }else{
            ssa_rxn_propensity[i][rxn_index] = k_rate/volume/2.0*C[i][reactants[0]]*C[i][reactants[1]];
            d_ssa_rxn_prop_d_c[i][rxn_index][reactants[0]] = k_rate/volume/2.0;
            d_ssa_rxn_prop_d_c[i][rxn_index][reactants[1]] = k_rate/volume/2.0;
            ssa_stoich_matrix[i][rxn_index][reactants[0]] = -1;
            ssa_stoich_matrix[i][rxn_index][reactants[1]] = -1;
        }
      }else if(num_reactants==1){
        ssa_rxn_propensity[i][rxn_index] = k_rate*C[i][reactants[0]];
        d_ssa_rxn_prop_d_c[i][rxn_index][reactants[0]] = k_rate;
        ssa_stoich_matrix[i][rxn_index][reactants[0]] = -1;
      }else if(num_reactants==0){
        ssa_rxn_propensity[i][rxn_index]  = k_rate*volume;
      }else{
        error->all(FLERR,"Illegal fix ssa_tsdpd_ssa_rxn_mass_action post_integrate, illegal number of reactants.");
      }
      if(num_products > 0){
        for(int j=0;j<num_products;j++){
          ssa_stoich_matrix[i][rxn_index][products[j]] = 0;
        }
        for(int j=0;j<num_products;j++){
          ssa_stoich_matrix[i][rxn_index][products[j]] += 1;
        }
      }
    }
  }
}
