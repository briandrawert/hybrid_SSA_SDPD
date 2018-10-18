/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(etd,FixETD)

#else

#ifndef LMP_FIX_ETD_H
#define LMP_FIX_ETD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixETD : public Fix {
 public:
  FixETD(class LAMMPS *, int, char **);
  virtual ~FixETD() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void post_integrate();
  virtual void final_integrate();
 
protected:
  double dtv,dtf;
  double *step_respa;
  int mass_require;
  class RanMars *random;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
