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

#ifdef PAIR_CLASS

PairStyle(ssa_tsdpd/purediffusion,PairSsaTsdpdPureDiffusion)

#else

#ifndef LMP_PAIR_SSA_TSDPD_PUREDIFFUSION_H
#define LMP_PAIR_SSA_TSDPD_PUREDIFFUSION_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSsaTsdpdPureDiffusion : public Pair {
 public:
  PairSsaTsdpdPureDiffusion(class LAMMPS *);
  virtual ~PairSsaTsdpdPureDiffusion();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double **cut, **alpha;
  void allocate();
};

}

#endif
#endif
