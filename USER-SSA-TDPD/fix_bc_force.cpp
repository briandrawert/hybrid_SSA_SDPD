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
   Contributing author: Zhen Li (zhen_li@brown.edu)
   the CRUNCH group, Division of Applied Mathematics, Brown University

   This fix was designed to apply boundary forces and fluxes in tDPD
------------------------------------------------------------------------- */

#include "mpi.h"
#include "time.h"
#include "math.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_bc_force.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "update.h"
#include "neighbor.h"
#include "random_mars.h"
#include "comm.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBCForce::FixBCForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{ 
	if (narg < 4) error->all(FLERR,"Illegal fix bcforce command");

	if (strcmp(arg[3],"FD") == 0)		index = 0;
	else if (strcmp(arg[3],"FC") == 0) 	index = 1;
	else if (strcmp(arg[3],"ccD") == 0) 	index = 2;
	else if (strcmp(arg[3],"ccN") == 0) 	index = 3;
	else error->all(FLERR,"Illegal fix bcforce/Name command");
	
	cut = atof( arg[4] );

        if (strcmp(arg[5],"lo") == 0)           lo_hi = 0;
	else if (strcmp(arg[5],"hi") == 0)      lo_hi = 1;
	else error->all(FLERR,"Illegal fix bcforce/lo_hi command");

	order = atoi( arg[6] );
	coe   = (double *) memory->smalloc(order*sizeof(double),"fix_bc_force:coe");

	switch (index){
		case 0 :
                        if (narg != 11+order ) error->all(FLERR,"Illegal fix bcforce/FD command");

		        if (strcmp(arg[7],"x") == 0)		dim = 0;
		        else if (strcmp(arg[7],"y") == 0)	dim = 1;
		        else if (strcmp(arg[7],"z") == 0)	dim = 2;
		        else error->all(FLERR,"Illegal fix bcforce/FD command");

		        if (strcmp(arg[8],"x") == 0)		dim_fix = 0;
		        else if (strcmp(arg[8],"y") == 0)	dim_fix = 1;
			else if (strcmp(arg[8],"z") == 0)	dim_fix = 2;
			else error->all(FLERR,"Illegal fix bcforce/FD command");

			for(int i=0; i<order; ++i)
				coe[i] = atof(arg[9+i]);

			kBT	= atof( arg[9 +order] );
			U0	= atof( arg[10 +order] );

		break;
		
		case 1 :
			if (narg != 8+order ) error->all(FLERR,"Illegal fix bcforce/FC command");
			
			if (strcmp(arg[7],"x") == 0)		dim = 0;
			else if (strcmp(arg[7],"y") == 0)	dim = 1;
			else if (strcmp(arg[7],"z") == 0)	dim = 2;
			else error->all(FLERR,"Illegal fix bcforce/FC command");

			for(int i=0; i<order; ++i)
				coe[i] = atof(arg[8+i]);

		break;

		case 2 :
			if (narg != 11+order ) error->all(FLERR,"Illegal fix bcforce/ccD command");

			if (strcmp(arg[7],"x") == 0)		dim = 0;
			else if (strcmp(arg[7],"y") == 0)	dim = 1;
			else if (strcmp(arg[7],"z") == 0)	dim = 2;
			else error->all(FLERR,"Illegal fix bcforce/ccD command");

			for(int i=0; i<order; ++i)
				coe[i] = atof(arg[8+i]);

			ctype = atoi( arg[8+order] );
			kappa = atof( arg[9+order] ); 
			C0    = atof( arg[10+order]);   

		break;

		case 3 :
			if (narg != 14+order ) error->all(FLERR,"Illegal fix bcforce/QV command");

			if (strcmp(arg[7],"x") == 0)		dim = 0;
			else if (strcmp(arg[7],"y") == 0)	dim = 1;
			else if (strcmp(arg[7],"z") == 0)	dim = 2;
			else error->all(FLERR,"Illegal fix bcforce/ccN command");

			if (strcmp(arg[8],"x") == 0)            dim_fix = 0;
                        else if (strcmp(arg[8],"y") == 0)       dim_fix = 1;
                        else if (strcmp(arg[8],"z") == 0)       dim_fix = 2;
                        else error->all(FLERR,"Illegal fix bcforce/FD command");


   			for(int i=0; i<order; ++i)
				coe[i] = atof(arg[9+i]);
			
			ctype	= atoi( arg[9+order] );
			
			U0	= atof( arg[10+order] );
			rho     = atof( arg[11+order] );
			loc_lo  = atof( arg[12+order] );
			loc_hi  = atof( arg[13+order] );

		break;
	}

	lo_h[0]= domain->boxlo[dim];
	lo_h[1]= domain->boxhi[dim]; 

	random = NULL;

	MPI_Barrier(world);

	if(index == 0){
		if(lo_hi==0) 		printf("caution:fix_bcforce/FD_lo, cpu is %d\n", comm->me);
		else if(lo_hi==1)	printf("caution:fix_bcforce/FD_hi, cpu is %d\n", comm->me);
		else error->all(FLERR,"Illegal fix bcforce/FD:lo_hi command");
	}
	else if(index == 1){
		if(lo_hi==0) 		printf("caution:fix_bcforce/FC_lo, cpu is %d\n", comm->me);
		else if(lo_hi==1)	printf("caution:fix_bcforce/FC_hi, cpu is %d\n", comm->me);
		else error->all(FLERR,"Illegal fix bcforce/FC:lo_hi command");
	}
	else if(index == 2){
		if(lo_hi==0) 		printf("caution:fix_bcforce/ccD_lo, cpu is %d\n", comm->me);
		else if(lo_hi==1)	printf("caution:fix_bcforce/ccD_hi, cpu is %d\n", comm->me);
		else error->all(FLERR,"Illegal fix bcforce/ccD:lo_hi command");
	}
	else if(index == 3){
		if(lo_hi==0) 		printf("caution:fix_bcforce/ccN_lo, cpu is %d\n", comm->me);
		else if(lo_hi==1)	printf("caution:fix_bcforce/ccN_hi, cpu is %d\n", comm->me);
		else error->all(FLERR,"Illegal fix bcforce/ccN:lo_hi command");
	}
	else error->all(FLERR,"Illegal fix bcforce/Name command");
		
	int seed=int(time(0));
	random = new RanMars(lmp, (seed+ comm->me) % 900000000 );
}

FixBCForce::~FixBCForce()
{
	memory->sfree(coe);	
	delete random;
}

int FixBCForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

void FixBCForce::init()
{

}


void FixBCForce::post_force(int vflag)
{
	int i, j;
	double **x = atom->x;
	double **f = atom->f;
	double **v = atom->v;
	double **C = atom->C;
	double **Q = atom->Q;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	double h, h1, h2, gammah, sigmah, fc, fd, fr, Qc, Qv, Qr, TT, prefactor, randnum;


	double delta_Qw = 0; 
	for (i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	{
	    h1 = x[i][dim]-lo_h[0];
	    h2 = lo_h[1]-x[i][dim];
		
            if( (lo_hi==0 && fabs(h1)<cut) || (lo_hi==1 && fabs(h2)<cut) ){ 
		randnum = random->gaussian();
		randnum = MAX(-3.0, MIN(randnum,3.0));

		if(lo_hi==0)
		    h = 1 - h1/cut;
		else if(lo_hi==1)
		    h = 1 - h2/cut;
		else error->all(FLERR,"Illegal fix bcforce/lo_hi command");
		
		if(index==0)
		{
			gammah = 0.0;
			for(j=0; j<order; ++j)
				gammah += coe[j]*pow(h,order-j); 

			if(h < 0.95) prefactor = 1 + 0.187*h/(1-h)-0.093*h*h*h;
			else if(h > 0.95) prefactor = 4.4733;
			else error->all(FLERR,"Illegal fix bcforce/lo_hi command, h > rc");

			gammah *= prefactor;

			sigmah = sqrt(2.0*kBT*fabs(gammah));
			fd   = -gammah * ( v[i][dim_fix]-U0 );
			fr   = sigmah * randnum;

			f[i][dim_fix]+=fd+fr;

		}
		else if(index==1)
		{
			fc = 0.0;
			for(j=0; j<order; ++j)
			    fc += coe[j]*pow(h,order-j);
			if(lo_hi==1) fc *= -1; 

			f[i][dim] += fc;
		}
		else if(index==2)
		{
			Qc = 0.0;
			for(j=0; j<order; ++j)
				Qc += coe[j]*pow(h,order-j);

                        if(h < 0.95) prefactor = 1 + 0.2673*h/(1-h)-0.2891*h*h*h;
                        else if(h > 0.95) prefactor = 5.8308;
                        else error->all(FLERR,"Illegal fix bcforce/lo_hi command, h > rc");

			Qc *= prefactor*kappa*(C0 - C[i][ctype]);

			Q[i][ctype] += Qc;
		}
		else if(index==3)
		{
			if(x[i][dim_fix] > loc_lo && x[i][dim_fix] < loc_hi)
			{
				Qv = 0.0;
				for(j=0; j<order; ++j)
			    		Qv += coe[j]*pow(h,order-j);
				Q[i][ctype] += U0*Qv/rho;
			}
		}
		else error->all(FLERR,"Illegal fix bcforce/name command");
	    }
	}
}

