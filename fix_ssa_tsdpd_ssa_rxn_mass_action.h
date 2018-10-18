#ifdef FIX_CLASS

FixStyle(ssa_tsdpd/ssa_rxn_mass_action,FixSsaTsdpdSsaRxnMassAction)

#else

#ifndef FIX_SSA_TSDPD_SSA_MASS_ACTION_H
#define FIX_SSA_TSDPD_SSA_MASS_ACTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTsdpdSsaRxnMassAction : public Fix {
 public:
  FixSsaTsdpdSsaRxnMassAction(class LAMMPS *, int, char **);
  virtual  ~FixSsaTsdpdSsaRxnMassAction();
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
  int itype;
  double volume;
};

}

#endif
#endif
