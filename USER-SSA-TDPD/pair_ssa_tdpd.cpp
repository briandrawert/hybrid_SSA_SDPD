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

/* ----------------------------------------------------------------------
   Contributing author: Kurt Smith (U Pittsburgh)
------------------------------------------------------------------------- */

//modified version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pair_ssa_tdpd.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairSsaTdpd::PairSsaTdpd(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairSsaTdpd::~PairSsaTdpd()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cutinv); //added

    memory->destroy(cutc);

    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(sigma);
    memory->destroy(s1); //added
    memory->destroy(kC); //added
    memory->destroy(kappa); //added
    memory->destroy(s2); //added

  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairSsaTdpd::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,dot_vi,dot_vj,wd,randnum,factor_dpd;
  double wC,wD,wR;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **C = atom->C; //added
  double **Q = atom->Q; //added
 
  double **Aetd = atom->Aetd; //added (ETD)
  double **Betd = atom->Betd; //added (ETD)
  double **Cetd = atom->Cetd; //added (ETD)

  int *type = atom->type;
  double *mass = atom->mass; //added
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost; //added
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;


  double dtinvsqrt = 1.0/sqrt(update->dt);

  const double kboltz = 1.0; //added
  const double r_on = 0.01*cut_global; //added
  const double r_off = cut_global; //added

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
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


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      dot_vi = delx*vxtmp + dely*vytmp + delz*vztmp;


      double SigmaIJ = sigma[itype][jtype];
      double GammaIJ = gamma[itype][jtype];
      double sf = s1[itype][jtype];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
	
        if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx*delvx + dely*delvy + delz*delvz;
	dot_vj = delx*v[j][0] + dely*v[j][1] + delz*v[j][2];
        //wd = 1.0 - r/cut[itype][jtype];
        randnum = random->gaussian();


	//added
	double wC,wR, WD;
	if (r >= r_off) {wC = wR = wD = 0.0;}
	else if (r < r_on) {wC = wR = wD = 1.0;}
	else
	{
		double r_rc = r* cutinv[itype][jtype];
		wC = 1.0 - r_rc;
		wC = MAX(0,wC);
		wR = pow(wC,0.5*sf);
		wD = wR*wR;
	} 



        // conservative force = a0 * wd
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;
	
	/*
        fpair = a0[itype][jtype]*wd; 
        fpair -= gamma[itype][jtype]*wd*wd*dot*rinv; 
        fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt; 
        fpair *= factor_dpd*rinv; 
	*/

	//added
	//conservative force
	double fC, eC;
	//fC = a0[itype][jtype]*wC;
	eC = -0.5*a0[itype][jtype]*cut[itype][jtype]*wC*wC;

	fpair = a0[itype][jtype]*wC; 

	//dissipative force
	fpair -= gamma[itype][jtype]*wD*dot*rinv; 
	
	//random force
        fpair += sigma[itype][jtype]*wR*randnum*dtinvsqrt; 
        fpair *= factor_dpd*rinv; 

        f[i][0] += delx*fpair; 
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        //printf("i = %d,  j = %d, nlocal = %d \n" , i, j, nlocal);

      
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }


       
	//Compute C_ETD (Eq. 28, Mai-Duy (2013))
	Cetd[i][0] += -gamma[itype][jtype]*wD*rinv*rinv*delx*delx/mass[itype] - 1e-16; //- EPSILON;
        Cetd[i][1] += -gamma[itype][jtype]*wD*rinv*rinv*dely*dely/mass[itype] - 1e-16; //- EPSILON;
        Cetd[i][2] += -gamma[itype][jtype]*wD*rinv*rinv*delz*delz/mass[itype] - 1e-16; // - EPSILON;


        //Compute A_ETD  (Eq. 28, Mai-Duy (2013))
	Aetd[i][0] += -gamma[itype][jtype]*wD*dely*delx*vytmp*rinv*rinv/mass[itype] - gamma[itype][jtype]*wD*delz*delx*vztmp*rinv*rinv/mass[itype];
        Aetd[i][1] += -gamma[itype][jtype]*wD*delx*dely*vxtmp*rinv*rinv/mass[itype] - gamma[itype][jtype]*wD*delz*dely*vztmp*rinv*rinv/mass[itype];
        Aetd[i][2] += -gamma[itype][jtype]*wD*delx*delz*vxtmp*rinv*rinv/mass[itype] - gamma[itype][jtype]*wD*dely*delz*vytmp*rinv*rinv/mass[itype];

        Aetd[i][0] += gamma[itype][jtype]*wD*dot_vj*delx*rinv*rinv/mass[itype];
        Aetd[i][1] += gamma[itype][jtype]*wD*dot_vj*dely*rinv*rinv/mass[itype];
        Aetd[i][2] += gamma[itype][jtype]*wD*dot_vj*delz*rinv*rinv/mass[itype];

        Aetd[i][0] += a0[itype][jtype]*wC*delx*rinv/mass[itype];
        Aetd[i][1] += a0[itype][jtype]*wC*dely*rinv/mass[itype];
        Aetd[i][2] += a0[itype][jtype]*wC*delz*rinv/mass[itype];


        //Compute B_ETD  (Eq. 28, Mai-Duy (2013))
        Betd[i][0] += sigma[itype][jtype]*wR*delx*dtinvsqrt*randnum*rinv/mass[itype];
        Betd[i][1] += sigma[itype][jtype]*wR*dely*dtinvsqrt*randnum*rinv/mass[itype];
        Betd[i][2] += sigma[itype][jtype]*wR*delz*dtinvsqrt*randnum*rinv/mass[itype];


        if (newton_pair || j < nlocal) {

	  Cetd[j][0] += -gamma[itype][jtype]*wD*rinv*rinv*delx*delx/mass[jtype] - 1e-16;
          Cetd[j][1] += -gamma[itype][jtype]*wD*rinv*rinv*dely*dely/mass[jtype] - 1e-16;
          Cetd[j][2] += -gamma[itype][jtype]*wD*rinv*rinv*delz*delz/mass[jtype] - 1e-16;

   	  Aetd[j][0] += (-gamma[itype][jtype]*wD*dely*delx*v[j][1]*rinv*rinv/mass[jtype] - gamma[itype][jtype]*wD*delz*delx*v[j][2]*rinv*rinv/mass[jtype]);
          Aetd[j][1] += (-gamma[itype][jtype]*wD*delx*dely*v[j][0]*rinv*rinv/mass[jtype] - gamma[itype][jtype]*wD*delz*dely*v[j][2]*rinv*rinv/mass[jtype]);
          Aetd[j][2] += (-gamma[itype][jtype]*wD*delx*delz*v[j][0]*rinv*rinv/mass[jtype] - gamma[itype][jtype]*wD*dely*delz*v[j][1]*rinv*rinv/mass[jtype]);

          Aetd[j][0] += gamma[itype][jtype]*wD*dot_vi*delx*rinv*rinv/mass[jtype];
          Aetd[j][1] += gamma[itype][jtype]*wD*dot_vi*dely*rinv*rinv/mass[jtype];
          Aetd[j][2] += gamma[itype][jtype]*wD*dot_vi*delz*rinv*rinv/mass[jtype];

          Aetd[j][0] -= a0[itype][jtype]*wC*delx*rinv/mass[jtype];
          Aetd[j][1] -= a0[itype][jtype]*wC*dely*rinv/mass[jtype];
          Aetd[j][2] -= a0[itype][jtype]*wC*delz*rinv/mass[jtype];

          Betd[j][0] -= sigma[itype][jtype]*wR*delx*dtinvsqrt*randnum*rinv/mass[jtype];
          Betd[j][1] -= sigma[itype][jtype]*wR*dely*dtinvsqrt*randnum*rinv/mass[jtype];
          Betd[j][2] -= sigma[itype][jtype]*wR*delz*dtinvsqrt*randnum*rinv/mass[jtype];

        }


       
	//chemical concentration transport //Zhen Li's algorithm
	if( r < cutc[itype][jtype])
        {
        	for(int k=0; k < atom->num_tdpd_species; ++k)
                {
                	double sC   = s2[itype][jtype][k];
                        double wCR = 1.0 - r/cutc[itype][jtype];
                        wCR = MAX(0,wCR);
                        wCR = pow(wCR, 0.5*sC);

                        double dQc  = -kappa[itype][jtype][k] * wCR*wCR * ( C[i][k] - C[j][k] );  // q_cond

                        Q[i][k] += (dQc );
                        if (newton_pair || j < nlocal)
                            Q[j][k]  -= ( dQc);

        	}
        }


          if (eflag) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
          evdwl *= factor_dpd;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             eC,0.0,fpair,delx,dely,delz);

	//if (evflag) ev_tally(i,j,nlocal,newton_pair,
        //                     evdwl,0.0,fpair,delx,dely,delz);

      }
    }



/*
if (update->ntimestep == 10) {
	printf(" i = %d, Q[i][0] = %f \n", i, Q[i][0] );
	getchar();
         	                           }
*/

  //double dtfm = update->dt / mass[itype];
  //f[i][0] = v[i][0]*exp(Cetd[i][0]*dtfm) + (exp(Cetd[i][0]*dtfm) - 1.0) * (Aetd[i][0] + Betd[i][0] ) / Cetd[i][0];
  //f[i][1] = v[i][1]*exp(Cetd[i][1]*dtfm) + (exp(Cetd[i][1]*dtfm) - 1.0) * (Aetd[i][1] + Betd[i][1] ) / Cetd[i][1];
  //f[i][2] = v[i][2]*exp(Cetd[i][2]*dtfm) + (exp(Cetd[i][2]*dtfm) - 1.0) * (Aetd[i][2] + Betd[i][2] ) / Cetd[i][2];

// randnum = random->gaussian();
// Betd[i][0] *= randnum;
// Betd[i][1] *= randnum;
// Betd[i][2] *= randnum;

  }

	


  if (vflag_fdotr) virial_fdotr_compute();

     
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSsaTdpd::allocate()
{

  int i,j;
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  
  memory->create(cutinv,n+1,n+1,"pair:cutinv"); //added
  memory->create(cutc,n+1,n+1,"pair:cutc"); //added

  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(s1,n+1,n+1,"pair:s1");  //added
  memory->create(kC,n+1,n+1,atom->num_tdpd_species,"pair:kC"); //added
  memory->create(kappa,n+1,n+1,atom->num_tdpd_species,"pair:kappa"); //added
  memory->create(s2,n+1,n+1,atom->num_tdpd_species,"pair:s2"); //added

  for (i = 0; i <= atom->ntypes; i++)
    for (j = 0; j <= atom->ntypes; j++)
      sigma[i][j] = gamma[i][j] = 0.0;

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSsaTdpd::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  //temperature = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[0]);
  seed = force->inumeric(FLERR,arg[1]);


  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) 
	{
	cut[i][j] = cut_global;
	cutinv[i][j] = 1.0 / cut[i][j];
	}
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSsaTdpd::coeff(int narg, char **arg)
{
  if (narg != 7+1+3*atom->num_tdpd_species)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double a0_one = force->numeric(FLERR,arg[2]);
  double gamma_one = force->numeric(FLERR,arg[3]);
  double sigma_one = force->numeric(FLERR,arg[4]);
  double s1_one = force->numeric(FLERR,arg[5]);
  double cut_one = force->numeric(FLERR,arg[6]);
  double cut_two = force->numeric(FLERR,arg[7]);
  double kC_one[atom->num_tdpd_species], kappa_one[atom->num_tdpd_species], s2_one[atom->num_tdpd_species];

 // double cut_one = cut_global;
//  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  //added
  for (int k=0; k<atom->num_tdpd_species; k++)
  {
  	kC_one[k]      = atof(arg[8+3*k]);
  	kappa_one[k]   = atof(arg[9+3*k]);
  	s2_one[k]      = atof(arg[10+3*k]);
  }


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      sigma[i][j] = sigma_one;
      s1[i][j] = s1_one;
      cut[i][j] = cut_one; //added
      cutc[i][j] = cut_two; //added
      for(int k=0; k<atom->num_tdpd_species; k++)  //added
      {
	  kC [i][j][k] = kC_one[k];
    	  kappa [i][j][k] = kappa_one[k];
          s2 [i][j][k] = s2_one[k];
      }

      cutinv [i][j]    = 1.0 / cut_one; //added
      setflag[i][j]    = 1;
      count++;   
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSsaTdpd::init_style()
{
printf("Here processor %d!!\n",comm->me);
printf("comm->ghost_velocity = %d ",comm->ghost_velocity);

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair ssa_tdpd requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers


  if (force->newton_pair == 0 && comm->me == 0) error->warning(FLERR,
      "Pair ssa_tdpd needs newton pair on for momentum and heat conservation");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSsaTdpd::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  //sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut[j][i] = cut[i][j];
  cutinv[j][i] = cutinv[i][j];
  cutc[j][i] = cutc[i][j];
  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j]; //added
  s1[j][i] = s1[i][j]; //added

  for(int k=0; k<atom->num_tdpd_species; k++) //added
  {
  	kC[j][i][k]      = kC[i][j][k];
 	kappa[j][i][k]   = kappa[i][j][k];
  	s2[j][i][k]      = s2[i][j][k];
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSsaTdpd::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) 
	{
	    fwrite(&a0     [i][j],sizeof(double),1,fp);
            fwrite(&gamma  [i][j],sizeof(double),1,fp);
            fwrite(&sigma  [i][j],sizeof(double),1,fp);
            fwrite(&s1     [i][j],sizeof(double),1,fp);
            fwrite(&cut    [i][j],sizeof(double),1,fp);
            fwrite(&cutc   [i][j],sizeof(double),1,fp);
            for(int k=0; k<atom->num_tdpd_species; k++)
            {
                fwrite(&kC   [i][j][k],sizeof(double),1,fp);
                fwrite(&kappa[i][j][k],sizeof(double),1,fp);
                fwrite(&s2   [i][j][k],sizeof(double),1,fp);
            }
     }
   }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSsaTdpd::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;


  for (i = 1; i <= atom->ntypes; i++)
  for (j = i; j <= atom->ntypes; j++) 
  {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) 
      {
        if (me == 0)
        {  
        	fread(&a0     [i][j],sizeof(double),1,fp);
                fread(&gamma  [i][j],sizeof(double),1,fp);
                fread(&sigma  [i][j],sizeof(double),1,fp);
                fread(&s1     [i][j],sizeof(double),1,fp);
                fread(&cut    [i][j],sizeof(double),1,fp);
                fread(&cutc   [i][j],sizeof(double),1,fp);
                for(int k=0; k<atom->num_tdpd_species; k++)
                {
                    fread(&kC   [i][j][k],sizeof(double),1,fp);
                    fread(&kappa[i][j][k],sizeof(double),1,fp);
                    fread(&s2   [i][j][k],sizeof(double),1,fp);
                }

 	}
	
        MPI_Bcast(&a0     [i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma  [i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma  [i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&s1     [i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut    [i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cutc   [i][j],1,MPI_DOUBLE,0,world);
	
	for(int k=0; k<atom->num_tdpd_species; k++)
            {
                MPI_Bcast(&kC   [i][j][k],1,MPI_DOUBLE,0,world);
                MPI_Bcast(&kappa[i][j][k],1,MPI_DOUBLE,0,world);
                MPI_Bcast(&s2   [i][j][k],1,MPI_DOUBLE,0,world);
            }
            cutinv[i][j] = 1.0 / cut[i][j];
         }
    }

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSsaTdpd::write_restart_settings(FILE *fp)
{
  //fwrite(&temperature,sizeof(double),1,fp); modified
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSsaTdpd::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    //fread(&temperature,sizeof(double),1,fp); modified
    fread(&cut_global,sizeof(double),1,fp);
    fread(&seed,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  //MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world); modified
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);


  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
  // seed = int(time(0));
  // random = new RanMars(lmp, (seed+ comm->me) % 900000000);
  //seed_h = (comm->me + 3)*seed;
  //seed48(&seed_h);

}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairSsaTdpd::write_data(FILE *fp)
{

  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,a0[i][i],gamma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairSsaTdpd::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,a0[i][j],gamma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairSsaTdpd::single(int i, int j, int itype, int jtype, double rsq,
                       double factor_coul, double factor_dpd, double &fforce)
{
  double r,rinv,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];
  fforce = a0[itype][jtype]*wd * factor_dpd*rinv;

  phi = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
  return factor_dpd*phi;

}
