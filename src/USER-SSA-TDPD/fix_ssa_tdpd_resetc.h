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

FixStyle(ssa_tdpd/resetC,FixSsaTdpdResetC)

#else

#ifndef FIX_SSA_TDPD_RESETC_H
#define FIX_SSA_TDPD_RESETC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTdpdResetC : public Fix {
 public:
  FixSsaTdpdResetC(class LAMMPS *, int, char **);
  virtual  ~FixSsaTdpdResetC();
  int setmask();
  virtual void init();
  virtual void post_integrate();

 protected:
  int step;
  int ctype;
  int index;
  double center[2], radius, length, width;
  double value;
};

}

#endif
#endif

