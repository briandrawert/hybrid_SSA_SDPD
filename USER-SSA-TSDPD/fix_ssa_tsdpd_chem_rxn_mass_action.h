#ifdef FIX_CLASS

FixStyle(ssa_tsdpd/chem_rxn_mass_action,FixSsaTsdpdChemRxnMassAction)

#else

#ifndef FIX_SSA_TSDPD_MASS_ACTION_H
#define FIX_SSA_TSDPD_MASS_ACTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTsdpdChemRxnMassAction : public Fix {
 public:
  FixSsaTsdpdChemRxnMassAction(class LAMMPS *, int, char **);
  virtual  ~FixSsaTsdpdChemRxnMassAction();
  int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  int rxn_index;
  int num_reactants;
  int num_products;
  int reactants[2];
  int products[4];
  double k_rate;
};

}

#endif
#endif
