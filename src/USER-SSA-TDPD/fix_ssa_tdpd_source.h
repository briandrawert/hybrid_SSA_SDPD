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

FixStyle(ssa_tdpd/source,FixSsaTdpdSource)

#else

#ifndef FIX_SSA_TDPD_SOURCE_H
#define FIX_SSA_TDPD_SOURCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTdpdSource : public Fix {
 public:
  FixSsaTdpdSource(class LAMMPS *, int, char **);
  ~FixSsaTdpdSource();
  virtual int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  int step;
  int typec;
  int wdim;
  int index;
  double center[2], radius, length, width;
  double value;
};

}

#endif
#endif
