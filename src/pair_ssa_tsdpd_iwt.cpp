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
 
#include <math.h>
#include <stdlib.h>
#include "pair_ssa_tsdpd_iwt.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "update.h"
#include "random_mars.h"
#include <set>
#include <unistd.h>
#include <time.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSsaTsdpdIwt::PairSsaTsdpdIwt(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairSsaTsdpdIwt::~PairSsaTsdpdIwt() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
    memory->destroy(kappa);
    memory->destroy(cutc);
  }
    if (random) delete random;
}


void PairSsaTsdpdIwt::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  
  //printf("PairSsaTsdpdIwt::compute() inum=%i\n",inum);
    
  
  double randnum;  

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq, velx, vely, velz;
  double rsq, tmp, wfd, wf, delVdotDelR, deltaE, mu;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;


  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *de = atom->de;
  double *e = atom->e;
  double *drho = atom->drho;
  double **C = atom->C;
  double **Q = atom->Q;
  int **Cd = atom->Cd;
  int **Qd = atom->Qd;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int dimension = domain->dimension;
  double kBoltzmann = force->boltz;
  double dtinv = 1.0 / update->dt;
  double *dfsp_D_matrix = atom->dfsp_D_matrix;
  double *dfsp_D_diag = atom->dfsp_D_diag;
  double *dfsp_Diffusion_coeff = atom->dfsp_Diffusion_coeff;
  double *dfsp_a_i = atom->dfsp_a_i;
  int nmax = atom->nmax;


  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 1.e-32) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SsaTsdpd particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
                  i, j, sqrt(cutsq[i][j]));
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


// New version of allocation (dfsp_D_matrix is now allocated at the atom_vec_ssa_tsdpd file)
  std::set<int>* dfsp_D_matrix_index = new std::set<int>[nmax];   // set for each column to store which elements are non-zero
  int inumsq = inum*inum;
  
   //allocate M
   double *M11 = new double[nmax];
   double *M12 = new double[nmax];
   double *M21 = new double[nmax];
   double *M22 = new double[nmax];

  
 //TODO: Create dfsp_D_diag, dfsp_Diffusion_coeff at atom_vec_ssa_tsdpd file. Also, we have to think in a better way of working with dfsp_D_matrix_index (std library seems slow)

 // loop over neighbors of my atoms

  //printf("\tStarting i loop\n");
  
  
  //compute kernel correction matrix
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii]; 
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

     for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
  
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];


      if (rsq < cutsq[itype][jtype] ) {
        h = cut[itype][jtype];     // for Lucy kernel
        double r = sqrt(rsq);


        if (domain->dimension == 3) { // Kernel, 3d (1/r * dwdr) 
          //Lucy kernel (3D)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
          wf = h - sqrt(rsq);
          wf  = 2.088908628081126 * wf * wf * wf * ihsq * ihsq * (h + 3.*r);

        } else if (domain->dimension == 2){ // Kernel, 2d (1/r * dwdr)
          ///*
          //Lucy kernel (2D)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
          wf = h - sqrt(rsq);
          wf  = 1.591549430918954 * wf * wf * wf * ihsq * ihsq * (h + 3.*r);
          //*/

          /*
          // Wendland C2 (2d)
          h = 0.5 * h;
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wf  = 2.*h - r;
          wfd = -0.348151438013521 * ihsq * ihsq * ihsq * ih * wf * wf * wf;
          wf  = 0.034815143801352 * ihsq * ihsq* ihsq * ih * wf * wf* wf * wf * (2.*r + h); 
          */

          /*
          // Wendland C4 (2d)
          h = 0.5 * h;
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq); 
          wf  = 2.*h - r;
          wfd = -0.007460387957433 * 7. * ihsq * ihsq * ihsq * ihsq * ihsq * wf * wf * wf * wf * wf *(2.*h + 5.*r);
          wf  = 0.011190581936149 * ihsq * ihsq* ihsq * ihsq * ihsq * wf * wf* wf * wf * wf * wf * ( (35./12.)*rsq + 3.*r*h + h*h); 
          */

          /*
          // Wendland C6 (2d)
          h = 0.5 * h;
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wf  = 2.*h - r;
          wfd = 0.886720397226274 * ihsq*ihsq*ihsq*ihsq*ihsq*ihsq*ih* ( -5.5*h*h*h*h*h*h*h*h*h + 16.5*h*h*h*h*h*h*h*rsq -43.3125*h*h*h*h*h*rsq*rsq + 57.75*h*h*h*h*rsq*rsq*r -36.0938*h*h*h*rsq*rsq*rsq +12.375*h*h*rsq*rsq*rsq*r -2.25586*h*rsq*rsq*rsq*rsq + 0.171875*rsq*rsq*rsq*rsq*r );
          wf  = 2.*h - r;
          wf  = 0.003463751551665 * ihsq * ihsq* ihsq * ihsq * ihsq * ihsq * ih * wf * wf* wf * wf * wf * wf * wf * wf * (h*h*h + 4.*r*h*h + 6.25*rsq*h + 4.*rsq*r);
          */

          /*
          // quintic Wendland (2d)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq); 
          wf  = 2.*h - r;
          wfd = - 0.0348151438013521 * ihsq * 10.* r * (2.*h-r) * (2.*h-r) * (2.*h-r) *ihsq*ihsq*ih/(r + 1e-12);
          wf  = 2.-r*ih;
          wf  = 0.0348151438013521 * ihsq * wf*wf*wf*wf * (1.+2.*r*ih);
          */

        } else if (domain->dimension == 1){ // Kernel, 1d (1/r * dwdr)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wfd = -15.0 * wfd * wfd * ihsq * ihsq * ih; //Lucy (1d)
        }

        
        //Kernel correction: Xu & Deng (2016)
        M11[i] += jmass*delx*delx*wf/rho[j];
        M12[i] += jmass*dely*delx*wf/rho[j];
        M21[i] += jmass*delx*dely*wf/rho[j];
        M22[i] += jmass*dely*dely*wf/rho[j];
       

        /*
        //Kernel correction: Oger et al. (2007)
        M11[i] += -jmass*delx*delx*wfd/rho[j];
        M12[i] += -jmass*delx*dely*wfd/rho[j];
        M21[i] += -jmass*dely*delx*wfd/rho[j];
        M22[i] += -jmass*dely*dely*wfd/rho[j];
        */
        
        if (newton_pair || j < nlocal) {
          ///*
          //Kernel correction: Xu & Deng (2016)
     	  M11[j] += imass*delx*delx*wf/rho[i];
          M12[j] += imass*dely*delx*wf/rho[i];
          M21[j] += imass*delx*dely*wf/rho[i];
          M22[j] += imass*dely*dely*wf/rho[i];
          //*/

          /*
          //Kernel correction: Oger et al. (2007)
     	  M11[j] += -imass*delx*delx*wfd/rho[i];
          M12[j] += -imass*delx*dely*wfd/rho[i];
          M21[j] += -imass*dely*delx*wfd/rho[i];
          M22[j] += -imass*dely*dely*wfd/rho[i];
          */
       }
     }
   }
 }
  
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii]; 

    //printf("\t\ti=%i\n",i);
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];


    // compute pressure of atom i with Tait EOS
    tmp = rho[i] / rho0[itype];
    fi = tmp * tmp * tmp;
    //fi = B[itype] * (fi * fi * tmp - 1.0)  / (rho[i] * rho[i]); 
    fi = B[itype] * (fi * fi * tmp - 1.0)  / (rho[i] * rho[j]); 


     for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
  
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];


      if (rsq < cutsq[itype][jtype] ) {
        h = cut[itype][jtype];     // for Lucy kernel
        double r = sqrt(rsq);


        if (domain->dimension == 3) { // Kernel, 3d (1/r * dwdr) 
          //Lucy kernel (3D)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
          wf = h - sqrt(rsq);
          wf  = 2.088908628081126 * wf * wf * wf * ihsq * ihsq * (h + 3.*r);

        } else if (domain->dimension == 2){ // Kernel, 2d (1/r * dwdr)
          ///*
          //Lucy kernel (2D)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
          wf = h - sqrt(rsq);
          wf  = 1.591549430918954 * wf * wf * wf * ihsq * ihsq * (h + 3.*r);
          //*/

          /*
          // Wendland C2 (2d)
          h = 0.5 * h;
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wf  = 2.*h - r;
          wfd = -0.348151438013521 * ihsq * ihsq * ihsq * ih * wf * wf * wf;
          wf  = 0.034815143801352 * ihsq * ihsq* ihsq * ih * wf * wf* wf * wf * (2.*r + h); 
          */

          /*
          // Wendland C4 (2d)
          h = 0.5 * h;
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq); 
          wf  = 2.*h - r;
          wfd = -0.007460387957433 * 7. * ihsq * ihsq * ihsq * ihsq * ihsq * wf * wf * wf * wf * wf *(2.*h + 5.*r);
          wf  = 0.011190581936149 * ihsq * ihsq* ihsq * ihsq * ihsq * wf * wf* wf * wf * wf * wf * ( (35./12.)*rsq + 3.*r*h + h*h); 
          */

          /*
          // Wendland C6 (2d)
          h = 0.5 * h;
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wf  = 2.*h - r;
          wfd = 0.886720397226274 * ihsq*ihsq*ihsq*ihsq*ihsq*ihsq*ih* ( -5.5*h*h*h*h*h*h*h*h*h + 16.5*h*h*h*h*h*h*h*rsq -43.3125*h*h*h*h*h*rsq*rsq + 57.75*h*h*h*h*rsq*rsq*r -36.0938*h*h*h*rsq*rsq*rsq +12.375*h*h*rsq*rsq*rsq*r -2.25586*h*rsq*rsq*rsq*rsq + 0.171875*rsq*rsq*rsq*rsq*r );
          wf  = 2.*h - r;
          wf  = 0.003463751551665 * ihsq * ihsq* ihsq * ihsq * ihsq * ihsq * ih * wf * wf* wf * wf * wf * wf * wf * wf * (h*h*h + 4.*r*h*h + 6.25*rsq*h + 4.*rsq*r);
          */

          /*
          // quintic Wendland (2d)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq); 
          wf  = 2.*h - r;
          wfd = - 0.0348151438013521 * ihsq * 10.* r * (2.*h-r) * (2.*h-r) * (2.*h-r) *ihsq*ihsq*ih/(r + 1e-12);
          wf  = 2.-r*ih;
          wf  = 0.0348151438013521 * ihsq * wf*wf*wf*wf * (1.+2.*r*ih);
          */

        } else if (domain->dimension == 1){ // Kernel, 1d (1/r * dwdr)
          ih = 1.0 / h;
          ihsq = ih * ih;
          wfd = h - sqrt(rsq);
          wfd = -15.0 * wfd * wfd * ihsq * ihsq * ih; //Lucy (1d)
        }



        //inverse of the kernel correction matrix Li
        double factor = (M11[i] * M22[i] - M12[i] * M21[i]);
        
        if (abs(factor)<1e-16) { //near singular, so use identity matrix
	  factor = 1;
          M11[i] = 1;
          M22[i] = 1;
          M12[i] = 0;
          M21[i] = 0;
        }
	
        double det_M = 1. / factor;
 
        double L11,L12,L21,L22;
        L11 = det_M * ( M22[i]);
        L12 = det_M * (-M12[i]);
        L21 = det_M * (-M21[i]);
        L22 = det_M * ( M11[i]);

        //correct delx and dely
        double delx_corr_i = (L11*delx + L12*dely);
        double dely_corr_i = (L21*delx + L22*dely);
        double delz_corr_i = delz;



        //inverse of the kernel correction matrix Lj
        factor = (M11[j] * M22[j] - M12[j] * M21[j]);
        
        if (abs(factor)<1e-16) { //near singular, so use identity matrix
	  factor = 1;
          M11[j] = 1;
          M22[j] = 1;
          M12[j] = 0;
          M21[j] = 0;
        }
	
        det_M = 1. / factor;
 
        L11 = det_M * ( M22[j]);
        L12 = det_M * (-M12[j]);
        L21 = det_M * (-M21[j]);
        L22 = det_M * ( M11[j]);

        //correct delx and dely
        double delx_corr_j = (L11*(-delx) + L12*(-dely));
        double dely_corr_j = (L21*(-delx) + L22*(-dely));
        double delz_corr_j = -delz;



        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
        //fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);
        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[i] * rho[j]);
        //if (fj < 0.0) fj = 0;

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;
	
	/*
        // Artificial viscosity (Managhan, 1992)
	if (delVdotDelR < 0.) {
	  mu = h * delVdotDelR / (rsq + 0.01 * h * h);
	  fvisc = -viscosity[itype][jtype] * (soundspeed[itype]
		  + soundspeed[jtype]) * mu / (rho[i] + rho[j]);
	} else {
	  fvisc = 0.;
	}
	fvisc *= imass * jmass * wfd / (rho[i] * rho[j]);
        */

	
        // Artificial viscosity (Managhan, 1992)
	if (delVdotDelR < 0.) {
	  mu = delVdotDelR / (rsq + 0.01 * h * h);
	  fvisc = -8.*viscosity[itype][jtype] * (soundspeed[itype]
		  + soundspeed[jtype]) * mu / (rho[i] + rho[j]) ;
	} else {
	  fvisc = 0.;
	}
	fvisc *= imass * jmass * wfd / ( 0.5*(rho[i] + rho[j]) * 0.5 *( soundspeed[itype] + soundspeed[jtype] ) );
        

        // total pair force
        fpair = imass * jmass * (-fi + fj) * wfd;

        
        // random force calculation
        // independent increments of a Wiener process matrix
        double wiener[3][3] = {0};
        for (int l=0; l<dimension; l++){
            for (int m=0; m<dimension; m++){
                wiener[l][m] = random->gaussian();
            }
        }


        // symmetric part
        wiener[0][1] = wiener[1][0] = (wiener[0][1] + wiener[1][0]) / 2.;
        wiener[0][2] = wiener[2][0] = (wiener[0][2] + wiener[2][0]) / 2.;
        wiener[1][2] = wiener[2][1] = (wiener[1][2] + wiener[2][1]) / 2.;

        // traceless part
        double trace_over_dim = (wiener[0][0] + wiener[1][1] + wiener[2][2]) / dimension;
        wiener[0][0] -= trace_over_dim;
        wiener[1][1] -= trace_over_dim;
        wiener[2][2] -= trace_over_dim;

        double prefactor = sqrt (-4. * kBoltzmann* e[i] * ( imass * jmass * wfd / (rho[i] * rho[j])  ) * dtinv) / (r+0.01*h);
        double f_random[3] = {0};

        for (int l=0; l<dimension; ++l)  f_random[l] = prefactor * (wiener[l][0]*delx + wiener[l][1]*dely + wiener[l][2]*delz);


        //Momentum evaluation
        ///*
        //kernel correction applied to the model of artificial viscosity (Monaghan, 1992) (Oger et al., 2007)
        f[i][0] += -delx_corr_i * fpair - delx_corr_i * fvisc ;//+ f_random[0];
        f[i][1] += -dely_corr_i * fpair - dely_corr_i * fvisc ;//+ f_random[1];
        f[i][2] += -delz_corr_i * fpair - delz_corr_i * fvisc ;//+ f_random[2];
        //*/
        /*
        //kernel correction applied to the model of artificial viscosity + XSPH term (Monaghan, 1992)
        double eps_xsph = 0.5;
        f[i][0] += -xcorr * fpair - xcorr * fvisc + f_random[0] - eps_xsph * imass*jmass*velx * wf/(0.5* (rho[i] + rho[j])); 
        f[i][1] += -ycorr * fpair - ycorr * fvisc + f_random[1] - eps_xsph * imass*jmass*vely * wf/(0.5* (rho[i] + rho[j]));
        f[i][2] += -delz  * fpair - delz  * fvisc + f_random[2] - eps_xsph * imass*jmass*vely * wf/(0.5* (rho[i] + rho[j]));
        */          
        
        //Density evaluation
        /*
        //artificial density diffusion: Molteni (2009) (disregards singularities)
        drho[i] += jmass * delVdotDelR * wfd - 0.1 * h * soundspeed[itype] * jmass * 2.0*(rho[j]/rho[i] - 1.0)  * wfd;
        */
        /*
        //kernel correction applied to the classical density formulation (Oger et al., 2007)
        drho[i] += rho[i] * jmass * (velx*xcorr + vely*ycorr + velz*delz) * wfd / rho[j];
        */
        ///*
        //kernel correction applied to the artificial density diffusion: Molteni (2009)
        drho[i] += rho[i] * jmass * (velx*delx_corr_i + vely*dely_corr_i + velz*delz_corr_i) * wfd / rho[j] - 0.1 * h * soundspeed[itype] * jmass * 2.0*( ((imass/rho[i]) / ( jmass/rho[j] )) - 1.0) * ( (delx*delx_corr_i + dely*dely_corr_i + delz*delz_corr_i ) /(rsq+0.01*h*h)) * wfd;
        //*/ 

        // Energy evaluation
        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));
        de[i] += deltaE;
        

        // Reactions in neighbors (j particles)
        if (newton_pair || j < nlocal) {


          //inverse of the kernel correction matrix
          double factor = (M11[j] * M22[j] - M12[j] * M21[j]);
        
          //Momentum evaluation
          ///*
          //kernel correction applied to the model of artificial viscosity (Monaghan, 1992) (Oger et al., 2007)
     	  f[j][0] += -delx_corr_j * (-fpair) - delx_corr_j * fvisc; //- f_random[0];
          f[j][1] += -dely_corr_j * (-fpair) - dely_corr_j * fvisc; //- f_random[1];
          f[j][2] += -delz_corr_j * (-fpair) - delz_corr_j * fvisc; //- f_random[2];
          //*/
          /*
          // kernel correction applied to the model of artificial viscosity + XSPH term (Monaghan, 1992)
          double eps_xsph = 0.5;
     	  f[j][0] -= (-xcorr * (-fpair) - xcorr * fvisc + f_random[0] - eps_xsph * imass*jmass*velx * wf/(0.5* (rho[i] + rho[j])) );  
          f[j][1] -= (-ycorr * (-fpair) - ycorr * fvisc + f_random[1] - eps_xsph * imass*jmass*vely * wf/(0.5* (rho[i] + rho[j])) );
          f[j][2] -= (-delz  * (-fpair) - delz  * fvisc + f_random[2] - eps_xsph * imass*jmass*vely * wf/(0.5* (rho[i] + rho[j])) );
          */         

          //Density evaluation
          /*
          // artificial density diffusion: Molteni (2009) (disregards singularities)
          drho[j] += imass * delVdotDelR * wfd - 0.1 * h * soundspeed[jtype] * imass * 2.0*(rho[i]/rho[j] - 1.0) * wfd; 
          */
          /*
          // kernel correction applied to the classical density formulation (Oger et al., 2007)
          drho[j] += rho[j] * imass * (velx*xcorr + vely*ycorr + velz*delz) * wfd / rho[i];
          */        
          ///*
          // kernel correction applied to the artificial density diffusion: Molteni (2009)
          drho[j] += rho[j] * imass * (-velx*delx_corr_j - vely*dely_corr_j - velz*delz_corr_j) * wfd / rho[i] - 0.1 * h * soundspeed[jtype] * imass * 2.0*( ((jmass/rho[j]) / ( imass/rho[i] )) - 1.0) * ( (-delx*delx_corr_j - dely*dely_corr_j - delz*delz_corr_j) /(rsq+0.01*h*h)) * wfd;
          //*/
                 

          //Energy evaluation
          de[j] += deltaE;
        }

        // transport of species
        if (r < cutc[itype][jtype]) {

          double r = sqrt(rsq);
          h = cutc[itype][jtype];

          if (domain->dimension == 3) { // Kernel, 3d (1/r * dwdr) 
            //Lucy kernel (3D)
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq);
            wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
            wf = h - sqrt(rsq);
            wf  = 2.088908628081126 * wf * wf * wf * ihsq * ihsq * (h + 3.*r);

          } else if (domain->dimension == 2){ // Kernel, 2d (1/r * dwdr)
            ///*
            //Lucy kernel (2D)
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq);
            wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
            wf = h - sqrt(rsq);
            wf  = 1.591549430918954 * wf * wf * wf * ihsq * ihsq * (h + 3.*r);
            //*/
            /*
            // Wendland C2 (2d)
            h = 0.5 * h;
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq);
            wf  = 2.*h - r;
            wfd = -0.348151438013521 * ihsq * ihsq * ihsq * ih * wf * wf * wf;
            wf  = 0.034815143801352 * ihsq * ihsq* ihsq * ih * wf * wf* wf * wf * (2.*r + h); 
            */
            /*
            // Wendland C4 (2d)
            h = 0.5 * h;
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq); 
            wf  = 2.*h - r;
            wfd = -0.007460387957433 * 7. * ihsq * ihsq * ihsq * ihsq * ihsq * wf * wf * wf * wf * wf *(2.*h + 5.*r);
            wf  = 0.011190581936149 * ihsq * ihsq* ihsq * ihsq * ihsq * wf * wf* wf * wf * wf * wf * ( (35./12.)*rsq + 3.*r*h + h*h); 
            */
            /*
            // Wendland C6 (2d)
            h = 0.5 * h;
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq);
            wf  = 2.*h - r;
            wfd = 0.886720397226274 * ihsq*ihsq*ihsq*ihsq*ihsq*ihsq*ih* ( -5.5*h*h*h*h*h*h*h*h*h + 16.5*h*h*h*h*h*h*h*rsq -43.3125*h*h*h*h*h*rsq*rsq + 57.75*h*h*h*h*rsq*rsq*r -36.0938*h*h*h*rsq*rsq*rsq +12.375*h*h*rsq*rsq*rsq*r -2.25586*h*rsq*rsq*rsq*rsq + 0.171875*rsq*rsq*rsq*rsq*r );
            wf  = 2.*h - r;
            wf  = 0.003463751551665 * ihsq * ihsq* ihsq * ihsq * ihsq * ihsq * ih * wf * wf* wf * wf * wf * wf * wf * wf * (h*h*h + 4.*r*h*h + 6.25*rsq*h + 4.*rsq*r);
            */
            /*
            // quintic Wendland (2d)
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq); 
            wf  = 2.*h - r;
            wfd = - 0.0348151438013521 * ihsq * 10.* r * (2.*h-r) * (2.*h-r) * (2.*h-r) *ihsq*ihsq*ih/(r + 1e-12);
            wf  = 2.-r*ih;
            wf  = 0.0348151438013521 * ihsq * wf*wf*wf*wf * (1.+2.*r*ih);
            */

          } else if (domain->dimension == 1){ // Kernel, 1d (1/r * dwdr)
            ih = 1.0 / h;
            ihsq = ih * ih;
            wfd = h - sqrt(rsq);
            wfd = -15.0 * wfd * wfd * ihsq * ihsq * ih; //Lucy (1d)
          }


          double dQc_base = 2.0* ((imass*jmass)/(imass+jmass)) * ((rho[i]+rho[j])/(rho[i]*rho[j])) * (delx*delx_corr_i + dely*dely_corr_i + delz*delz_corr_i) * wfd  / (rsq + 0.01*h*h); // (Tartakovsky et. al., 2007, JCP)


          //printf("\t\t\tAssigning D matrix[%i][%i]\n",i,j);
          if (atom->num_ssa_species>0){              
            dfsp_D_matrix[i*inum+j] = - dQc_base;
            dfsp_D_matrix_index[i].insert(j);
            dfsp_D_matrix[j*inum+i] = - dQc_base;
            dfsp_D_matrix_index[j].insert(i);
          }


	  // TODO: dfsp_Diffusion_coeff seems redundant here. kappa[itype][jtype][s] is already allocated.
          for(int s=0;s<atom->num_ssa_species;s++){
            dfsp_Diffusion_coeff[i*inum+j+s*inumsq] = kappa[itype][jtype][s];
            dfsp_Diffusion_coeff[j*inum+i+s*inumsq] = kappa[itype][jtype][s];
         }
            
        
          for(int k=0; k < atom->num_tdpd_species; ++k){
            double dQc = (kappa[itype][jtype][k]) * ( C[i][k] - C[j][k] ) * dQc_base;
            Q[i][k] += (dQc);
            if (newton_pair || j < nlocal) {

              double dQc_base_j = 2.0* ((imass*jmass)/(imass+jmass)) * ((rho[i]+rho[j])/(rho[i]*rho[j])) * (-delx*delx_corr_j - dely*dely_corr_j - delz*delz_corr_j) * wfd  / (rsq + 0.01*h*h); // (Tartakovsky et. al., 2007, JCP)
              double dQc_j = -(kappa[itype][jtype][k]) * ( C[i][k] - C[j][k] ) * dQc_base;
              Q[j][k] += (dQc_j);
            } 
          }
        }
  
 
        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }

   }

  }

  
  if (vflag_fdotr) virial_fdotr_compute();
  //printf("\tend of  i loop\n");


  
  //printf("Starting SSA diffusion\n");
  // Second Step: Calculate SSA Diffusion
  // Find Diagional element values
  if (atom->num_ssa_species > 0){
    std::set<int>::iterator it;
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      double total = 0;
      for (it=dfsp_D_matrix_index[i].begin(); it!=dfsp_D_matrix_index[i].end(); ++it){  // Using a set::iterator
          total += *it;
      }
      dfsp_D_diag[i] = total;
    }    

  double tt,a0,sum_d,sum_d2,r1,r2,r3,a_i_src_orig,a_i_dest_orig;
//  double* a_i = new double[inum];
  int k;
  for(int s=0;s<atom->num_ssa_species;s++){  // Calculate each species seperatly
    // sum each voxel propensity to get total propensity
    a0 = 0;
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      dfsp_a_i[i] = 0.0;
      for (it=dfsp_D_matrix_index[i].begin(); it!=dfsp_D_matrix_index[i].end(); ++it){  // iterator over outbound connections
        j = *it;

        dfsp_a_i[i] += dfsp_Diffusion_coeff[i*inum+j+s*inumsq] * dfsp_D_matrix[i*inum+j]; // Note, "dfsp_a_i" is the per-voxel base propensity (must multiply Cd[i][s]);
      }
      a0 += dfsp_a_i[i] * Cd[i][s];
    }
    // Find time to first reaction
    tt=0;
    r1 = random->uniform();
    tt += -log(1.0-r1)/a0 ;
    //printf("SSA_Diffusion[s=%i] Cd[%i][%i]=%i tt=%e a0=%e r1=%e\n",s,i,s,Cd[i][s],tt,a0,r1);
    // Loop over time
    while(tt < update->dt){
        // find which voxel the diffusion event occured
        r2 = a0 * random->uniform();
        sum_d = 0;
        for(k=0; k < inum; k++){
            sum_d += dfsp_a_i[k] * Cd[k][s];
            if(sum_d > r2) break;
            //printf("\t k=%i sum_d=%e r2=%e\n",k,sum_d,r2);
        }
        int src_vox = k;
        // find which voxel it moved to
        r3 = dfsp_a_i[k] * random->uniform();
        sum_d2=0;
        for (it=dfsp_D_matrix_index[src_vox].begin(); it!=dfsp_D_matrix_index[src_vox].end(); ++it){
            j = *it;
            sum_d2 += dfsp_Diffusion_coeff[src_vox*inum+j+s*inumsq] * dfsp_D_matrix[src_vox*inum+j];
            if(sum_d2 > r3) break;
        }
        int dest_vox = j;
        // Move molecule  //TODO: which way is better? Computing Cd directly, or passing molecule diffusion via Qd? (I guess Qd is better).
        //Cd[src_vox][s]--;
        //Cd[dest_vox][s]++;
   	Qd[src_vox][s]--;
        Qd[dest_vox][s]++;
        // Find delta in propensities
        a0 += (dfsp_a_i[dest_vox] - dfsp_a_i[src_vox]);
        
        //printf("SSA_Diffusion[s=%i] tt=%e %i -> %i.  Now a0=%e",s,tt,src_vox,dest_vox,a0);
        // Find time to next reaction
        r1 = random->uniform();
        tt += -log(1.0-r1)/(a0);
        //printf(" tt=%e\n",tt);
    }
  }

//  delete[] a_i;
  
  }
  // Dealocate the arrays.  TODO: move to the destructor
  delete[] dfsp_D_matrix_index;
  delete[] M11;
  delete[] M12;
  delete[] M21;
  delete[] M22;
  

}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSsaTsdpdIwt::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");

  memory->create(kappa,n+1,n+1, atom->num_tdpd_species + atom->num_ssa_species,"pair:kappa"); //added  
  memory->create(cutc, n + 1, n + 1, "pair:cutc");

}

/* ----------------------------------------------------------------------
   global settings
 ------------------------------------------------------------------------- */

void PairSsaTsdpdIwt::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style ssa_tsdpd/iwt");

  // seed is immune to underflow/overflow because it is unsigned
//  seed = comm->nprocs + comm->me + atom->nlocal;
//  if (narg == 3) seed += force->inumeric (FLERR, arg[2]);
//  random = new RanMars (lmp, seed);

  srand(clock());
  seed = comm->nprocs + comm->me + atom->nlocal + rand()%100;
  random = new RanMars (lmp,seed);

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSsaTsdpdIwt::coeff(int narg, char **arg) {
  if (narg != 7 + atom->num_tdpd_species + atom->num_ssa_species)
    error->all(FLERR,
        "Incorrect args for pair_style ssa_tsdpd/iwt coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  double rho0_one = force->numeric(FLERR,arg[2]);
  double soundspeed_one = force->numeric(FLERR,arg[3]);
  double viscosity_one = force->numeric(FLERR,arg[4]);
  double cut_one = force->numeric(FLERR,arg[5]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;
  double cutc_one = force->numeric(FLERR,arg[6]);

  //added
  double kappa_one[atom->num_tdpd_species + atom->num_ssa_species];
  for (int k=0; k < atom->num_tdpd_species + atom->num_ssa_species; k++){
    kappa_one[k] = atof(arg[7+k]);
  }


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      cut[i][j] = cut_one;
      cutc[i][j] = cutc_one;

      for (int k=0; k < atom->num_tdpd_species + atom->num_ssa_species; k++) {
        kappa[i][j][k] = kappa_one[k];
      }

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSsaTsdpdIwt::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair ssa_tsdpd/iwt coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];

  for(int k=0; k < atom->num_tdpd_species + atom->num_ssa_species; k++) {
    kappa[j][i][k] = kappa[i][j][k];
  }

  cutc[j][i] = cutc[i][j];


  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSsaTsdpdIwt::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
