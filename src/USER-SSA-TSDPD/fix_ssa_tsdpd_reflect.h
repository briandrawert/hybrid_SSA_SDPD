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

FixStyle(ssa_tsdpd/reflect,FixSsaTsdpdReflect)

#else

#ifndef FIX_SSA_TSDPD_REFLECT_H
#define FIX_SSA_TSDPD_REFLECT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTsdpdReflect : public Fix {
 public:
  FixSsaTsdpdReflect(class LAMMPS *, int, char **);
  int setmask();
  void post_integrate();

 private:
  double lvalue, hvalue;
  int xflag,yflag,zflag;
  int dim;
};

}

#endif
#endif
