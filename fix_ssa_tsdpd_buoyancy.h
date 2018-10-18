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

FixStyle(ssa_tsdpd/buoyancy,FixSsaTsdpdBuoyancy)

#else

#ifndef FIX_SSA_TSDPD_BUOYANCY_H
#define FIX_SSA_TSDPD_BUOYANCY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTsdpdBuoyancy : public Fix {
 public:
  FixSsaTsdpdBuoyancy(class LAMMPS *, int, char **);
  ~FixSsaTsdpdBuoyancy();
  int setmask();
  void setup(int);
  void post_force(int);

 private:
  int boussinesq_ssa_flag, boussinesq_sdpd_flag, gravity_flag, rank_coordinate, rank_buoyancy;
  int xflag, yflag, zflag;
  double acceleration,C_ref;
};

}

#endif
#endif
