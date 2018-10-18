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
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_stat.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStat::FixStat(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 10+1+1 && narg != 16+1+1) error->all(FLERR,"Illegal fix stat command");

  nx = atoi(arg[3]);
  ny = atoi(arg[4]);
  nz = atoi(arg[5]);
  st_start = atoi(arg[6]);
  nevery = atoi(arg[7]);
  dump_each = atoi(arg[8]);
  if ((nx < 1)||(ny < 1)||(nz < 1)) error->all(FLERR,"Illegal division for statistics");
  if (narg == 16){
    xlo = atof(arg[9]);
    xhi = atof(arg[10]);
    ylo = atof(arg[11]);
    yhi = atof(arg[12]);
    zlo = atof(arg[13]);
    zhi = atof(arg[14]);
    sprintf(fname,arg[15]);
    index = atoi(arg[16]);
    nrepeat = atoi(arg[17]);
  }
  else{
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    sprintf(fname,arg[9]);
    index = atoi(arg[10]);
    nrepeat = atoi(arg[11]);
  }
  if ((xlo >= xhi)||(ylo >= yhi)||(zlo >= zhi)) error->all(FLERR,"Illegal coordinates for statistics");
  xs = xhi - xlo;
  ys = yhi - ylo;
  zs = zhi - zlo;
  xper = domain->xperiodic;
  yper = domain->yperiodic;
  zper = domain->zperiodic;
  dxlo = domain->boxlo[0];
  dxhi = domain->boxhi[0];
  dylo = domain->boxlo[1];
  dyhi = domain->boxhi[1];
  dzlo = domain->boxlo[2];
  dzhi = domain->boxhi[2];
  dxs = domain->xprd;
  dys = domain->yprd;
  dzs = domain->zprd;

  num_step = 0;
}

/* ---------------------------------------------------------------------- */

int FixStat::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixStat::map_index(double x, double y, double z)
{
  int ind = 0;

  if (x<dxlo || x>=dxhi || y<dylo || y>=dyhi || z<dzlo || z>=dzhi){
    if (xper) {
      while (x >= dxhi)
        x -= dxs;
      while (x < dxlo)
        x += dxs;
    }
    if (yper) {
      while (y >= dyhi)
        y -= dys;
      while (y < dylo)
        y += dys;
    }
    if (zper) {
      while (z >= dzhi)
        z -= dzs;
      while (z < dzlo)
        z += dzs;
    }
  }
  if (x>=xlo && x<xhi && y>=ylo && y<yhi && z>=zlo && z<zhi){
    is = static_cast<int> ((x - xlo)*nx/xs);
    js = static_cast<int> ((y - ylo)*ny/ys);
    ks = static_cast<int> ((z - zlo)*nz/zs); 
    ind = 1;
    if (is > nx-1) is--;
    if (js > ny-1) js--;
    if (ks > nz-1) ks--;
  }
  return ind;
}

/* ---------------------------------------------------------------------- */

int FixStat::map_index_cyl(double x, double y, double z)
{
  double theta, theta1,rr;
  int ind = 0;
  
  if (x<dxlo || x>=dxhi || y<dylo || y>=dyhi || z<dzlo || z>=dzhi){
    if (xper) {
      while (x >= dxhi) 
        x -= dxs;
      while (x < dxlo) 
        x += dxs;
    } 
    if (yper) {
      while (y >= dyhi) 
        y -= dys;
      while (y < dylo) 
        y += dys;
    }
    if (zper) {
      while (z >= dzhi) 
        z -= dzs;
      while (z < dzlo) 
        z += dzs;
    }
  }
   
  rr = sqrt((y-ylo)*(y-ylo) + (z-zlo)*(z-zlo));
  
  if (x>=xlo && x<xhi && rr<yhi){ 
    theta = acos((y-ylo)/rr);
    theta1 = asin((z-zlo)/rr);
    if (theta1 < 0.0)
      theta = 2.0*M_PI - theta; 
    is = static_cast<int> ((x - xlo)*nx/xs);
    js = static_cast<int> (rr*ny/yhi);
    ks = static_cast<int> (0.5*theta*nz/M_PI);
    ind = 1;
    if (is > nx-1) is--;
    if (js > ny-1) js--;
    if (ks > nz-1) ks--;
  }
  return ind;
} 

/* ---------------------------------------------------------------------- */


