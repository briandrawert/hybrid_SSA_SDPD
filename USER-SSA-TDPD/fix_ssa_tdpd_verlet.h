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
FixStyle(ssa_tdpd_verlet,FixSsaTdpdVerlet)

#else

#ifndef FIX_SSA_TDPD_VERLET_H
#define FIX_SSA_TDPD_VERLET_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSsaTdpdVerlet : public Fix {
 public:
  FixSsaTdpdVerlet(class LAMMPS *, int, char **);
  ~FixSsaTdpdVerlet();
  virtual int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:

};

}

#endif
#endif
