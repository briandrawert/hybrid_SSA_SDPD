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

FixStyle(ssa_tdpd/velT,FixSsaTdpdVelT)

#else

#ifndef FIX_SSA_TDPD_VELT_H
#define FIX_SSA_TDPD_VELT_H

#include "fix_stat.h"

namespace LAMMPS_NS {

class FixSsaTdpdVelT : public FixStat {
 public:
  FixSsaTdpdVelT(class LAMMPS *, int, char **);
  ~FixSsaTdpdVelT();
  int setmask();
  void end_of_step();
  void write_stat(int);

 private:
  double ***num;
  double ***vx, ***vy, ***vz;
  double ***aC;
};

}

#endif
#endif