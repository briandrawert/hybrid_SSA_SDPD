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

#ifdef FIX_CLASS

FixStyle(bc_force,FixBCForce)

#else

#ifndef FIX_BC_FORCE_H
#define FIX_BC_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBCForce : public Fix {
 public:
  FixBCForce(class LAMMPS *, int, char **);
  ~FixBCForce();
  virtual int setmask();
  virtual void init();
  virtual void post_force(int);


 private:
  int index, order, dim, dim_fix;
  int lo_hi, ctype;
  double kBT, U0, C0, kappa, rho, loc_lo, loc_hi;
  double lo_h[2];
  double *coe;
  double cut; 

  class RanMars *random;

};

}

#endif
#endif
