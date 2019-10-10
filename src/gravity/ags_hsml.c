#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/*! \file ags_hsml.c
 *  \brief kernel length determination for non-gas particles
 *
 *  This file contains a loop modeled on the gas density computation which 
 *    determines softening lengths (and appropriate correction terms) 
 *    for all particle types, to make softenings fully adaptive
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#define AGS_DSOFT_TOL (0.75)    // amount by which softening lengths are allowed to vary in single timesteps //

/*! this routine is called by the adaptive gravitational softening neighbor search and forcetree (for application 
    of the appropriate correction terms), to determine which particle types "talk to" which other particle types 
    (i.e. which particle types you search for to determine the softening radii for gravity). For effectively volume-filling
    fluids like gas or dark matter, it makes sense for this to be 'matched' to particles of the same type. For other 
    particle types like stars or black holes, it's more ambiguous, and requires some judgement on the part of the user. 
    The routine specifically returns a bitflag which defines all valid particles to which a particle of type 'primary' 
    can 'see': i.e. SUM(2^n), where n are all the particle types desired for neighbor finding,
    so e.g. if you want particle types 0 and 4, set the bitmask = 17 = 1 + 16 = 2^0 + 2^4
 */
int ags_gravity_kernel_shared_BITFLAG(short int particle_type_primary)
{
    /* gas particles see gas particles */
    if(particle_type_primary == 0) {return 1;}

#ifdef ADAPTIVE_GRAVSOFT_FORALL
#ifdef BLACK_HOLES
    /* black hole particles see gas */
    if(particle_type_primary == 5) {return 1;}
#endif
#ifdef GALSF
    /* stars see baryons (any type) */
    if(All.ComovingIntegrationOn)
    {
        if(particle_type_primary == 4) {return 17;} // 2^0+2^4
    } else {
        if((particle_type_primary == 4)||(particle_type_primary == 2)||(particle_type_primary == 3)) {return 29;} // 2^0+2^2+2^3+2^4
    }
#endif
#ifdef DM_SIDM
    /* SIDM particles see other SIDM particles */
    if((1 << particle_type_primary) & (DM_SIDM)) {return DM_SIDM;}
#endif
    /* if we haven't been caught by one of the above checks, we simply return whether or not we see 'ourselves' */
    return (1 << particle_type_primary);
#endif
    
    return 0;
}



#ifdef ADAPTIVE_GRAVSOFT_FORALL

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct ags_densdata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat AGS_Hsml;
  int NodeList[NODELISTLENGTH];
  int Type;
}
 *AGS_DensDataIn, *AGS_DensDataGet;

static struct ags_densdata_out
{
    MyLongDouble Ngb;
    MyLongDouble DhsmlNgb;
    MyLongDouble AGS_zeta;
    MyLongDouble AGS_vsig;
    MyLongDouble Particle_DivVel;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    MyLongDouble NV_T[3][3];
#endif
}
 *AGS_DensDataResult, *AGS_DensDataOut;

void ags_particle2in_density(struct ags_densdata_in *in, int i);
void ags_out2particle_density(struct ags_densdata_out *out, int i, int mode);

void ags_particle2in_density(struct ags_densdata_in *in, int i)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    in->AGS_Hsml = PPP[i].AGS_Hsml;
    in->Type = P[i].Type;
}

void ags_out2particle_density(struct ags_densdata_out *out, int i, int mode)
{
    ASSIGN_ADD(PPP[i].NumNgb, out->Ngb, mode);
    ASSIGN_ADD(PPPZ[i].AGS_zeta, out->AGS_zeta,   mode);
    if(out->AGS_vsig > PPP[i].AGS_vsig) {PPP[i].AGS_vsig = out->AGS_vsig;}
    ASSIGN_ADD(P[i].Particle_DivVel, out->Particle_DivVel,   mode);
    ASSIGN_ADD(PPP[i].DhsmlNgbFactor, out->DhsmlNgb, mode);
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    {int j,k; for(k = 0; k < 3; k++) {for(j = 0; j < 3; j++) {ASSIGN_ADD(P[i].NV_T[k][j], out->NV_T[k][j], mode);}}}
#endif
}

struct kernel_density
{
    double dp[3],dv[3],r;
    double wk, dwk;
    double hinv, hinv3, hinv4;
};


void ags_density(void)
{
  MyFloat *Left, *Right, *AGS_Prev;
  int i, j, k, ndone, ndone_flag, npleft, iter = 0;
  int ngrp, recvTask, place;
  long long ntot;
  double fac, fac_lim;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported = 0;
  int redo_particle;
  int particle_set_to_minhsml_flag = 0;
  int particle_set_to_maxhsml_flag = 0;

  CPU_Step[CPU_AGSDENSMISC] += measure_time();
  AGS_Prev = (MyFloat *) mymalloc("AGS_Prev", NumPart * sizeof(MyFloat));
    
  long long NTaskTimesNumPart;
  NTaskTimesNumPart = maxThreads * NumPart;
  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(ags_density_isactive(i))
      {
          Left[i] = Right[i] = 0;
          AGS_Prev[i] = PPP[i].AGS_Hsml;
          PPP[i].AGS_vsig = 0;
#ifdef WAKEUP
          P[i].wakeup = 0;
#endif
      }
    }

  /* allocate buffers to arrange communication */
  size_t MyBufferSize = All.BufferSize;
  All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct ags_densdata_in) + sizeof(struct ags_densdata_out) +
					     sizemax(sizeof(struct ags_densdata_in),sizeof(struct ags_densdata_out))));
  DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = my_second();

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {

      NextParticle = FirstActiveParticle;	/* begin with this index */

      do
	{
	  BufferFullFlag = 0;
	  Nexport = 0;
	  save_NextParticle = NextParticle;

	  tstart = my_second();

#ifdef PTHREADS_NUM_THREADS
	  pthread_t mythreads[PTHREADS_NUM_THREADS - 1];

	  int threadid[PTHREADS_NUM_THREADS - 1];

	  pthread_attr_t attr;

	  pthread_attr_init(&attr);
	  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	  pthread_mutex_init(&mutex_nexport, NULL);
	  pthread_mutex_init(&mutex_partnodedrift, NULL);

	  TimerFlag = 0;

	  for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
	    {
	      threadid[j] = j + 1;
	      pthread_create(&mythreads[j], &attr, ags_density_evaluate_primary, &threadid[j]);
	    }
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
#ifdef _OPENMP
	    int mainthreadid = omp_get_thread_num();
#else
	    int mainthreadid = 0;
#endif
	    ags_density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
	  }

#ifdef PTHREADS_NUM_THREADS
	  for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);
#endif

	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);

	  if(BufferFullFlag)
	    {
	      int last_nextparticle = NextParticle;

	      NextParticle = save_NextParticle;

	      while(NextParticle >= 0)
		{
		  if(NextParticle == last_nextparticle)
		    break;

		  if(ProcessedFlag[NextParticle] != 1)
		    break;

		  ProcessedFlag[NextParticle] = 2;

		  NextParticle = NextActiveParticle[NextParticle];
		}

	      if(NextParticle == save_NextParticle)
		{
		  /* in this case, the buffer is too small to process even a single particle */
		  printf("ags-Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
			 P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);

		  endrun(111008);
		}


	      int new_export = 0;

	      for(j = 0, k = 0; j < Nexport; j++)
		if(ProcessedFlag[DataIndexTable[j].Index] != 2)
		  {
		    if(k < j + 1)
		      k = j + 1;

		    for(; k < Nexport; k++)
		      if(ProcessedFlag[DataIndexTable[k].Index] == 2)
			{
			  int old_index = DataIndexTable[j].Index;

			  DataIndexTable[j] = DataIndexTable[k];
			  DataNodeList[j] = DataNodeList[k];
			  DataIndexTable[j].IndexGet = j;
			  new_export++;

			  DataIndexTable[k].Index = old_index;
			  k++;
			  break;
			}
		  }
		else
		  new_export++;

	      Nexport = new_export;

	    }


	  n_exported += Nexport;

	  for(j = 0; j < NTask; j++)
	    Send_count[j] = 0;
	  for(j = 0; j < Nexport; j++)
	    Send_count[DataIndexTable[j].Task]++;

	  MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

	  tstart = my_second();

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

	  tend = my_second();
	  timewait1 += timediff(tstart, tend);

	  for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      Nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  AGS_DensDataGet = (struct ags_densdata_in *) mymalloc("AGS_DensDataGet", Nimport * sizeof(struct ags_densdata_in));
	  AGS_DensDataIn = (struct ags_densdata_in *) mymalloc("AGS_DensDataIn", Nexport * sizeof(struct ags_densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      ags_particle2in_density(&AGS_DensDataIn[j], place);

	      memcpy(AGS_DensDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	    }
	  /* exchange particle data */
	  tstart = my_second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&AGS_DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct ags_densdata_in), MPI_BYTE,
				   recvTask, TAG_AGS_DENS_A,
				   &AGS_DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct ags_densdata_in), MPI_BYTE,
				   recvTask, TAG_AGS_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = my_second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(AGS_DensDataIn);
	  AGS_DensDataResult = (struct ags_densdata_out *) mymalloc("AGS_DensDataResult", Nimport * sizeof(struct ags_densdata_out));
	  AGS_DensDataOut = (struct ags_densdata_out *) mymalloc("AGS_DensDataOut", Nexport * sizeof(struct ags_densdata_out));

	  /* now do the particles that were sent to us */

	  tstart = my_second();

	  NextJ = 0;

#ifdef PTHREADS_NUM_THREADS
	  for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
	    pthread_create(&mythreads[j], &attr, ags_density_evaluate_secondary, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
#ifdef _OPENMP
	    int mainthreadid = omp_get_thread_num();
#else
	    int mainthreadid = 0;
#endif
	    ags_density_evaluate_secondary(&mainthreadid);
	  }

#ifdef PTHREADS_NUM_THREADS
	  for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
	    pthread_join(mythreads[j], NULL);

	  pthread_mutex_destroy(&mutex_partnodedrift);
	  pthread_mutex_destroy(&mutex_nexport);
	  pthread_attr_destroy(&attr);
#endif

	  tend = my_second();
	  timecomp2 += timediff(tstart, tend);

	  if(NextParticle < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  tstart = my_second();
	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  tend = my_second();
	  timewait2 += timediff(tstart, tend);


	  /* get the result */
	  tstart = my_second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&AGS_DensDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct ags_densdata_out),
				   MPI_BYTE, recvTask, TAG_AGS_DENS_B,
				   &AGS_DensDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct ags_densdata_out),
				   MPI_BYTE, recvTask, TAG_AGS_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = my_second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = my_second();
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
	      ags_out2particle_density(&AGS_DensDataOut[j], place, 1);
	    }
	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(AGS_DensDataOut);
	  myfree(AGS_DensDataResult);
	  myfree(AGS_DensDataGet);
	}
      while(ndone < NTask);


      /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(ags_density_isactive(i))
            {
                if(PPP[i].NumNgb > 0)
                {
                    PPP[i].DhsmlNgbFactor *= PPP[i].AGS_Hsml / (NUMDIMS * PPP[i].NumNgb);
                    P[i].Particle_DivVel /= PPP[i].NumNgb;
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    PPP[i].NumNgb *= NORM_COEFF * pow(PPP[i].AGS_Hsml,NUMDIMS);
                } else {
                    PPP[i].NumNgb = PPP[i].DhsmlNgbFactor = P[i].Particle_DivVel = 0;
                }
                
                // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(PPP[i].DhsmlNgbFactor > -0.5)	/* note: this would be -1 if only a single particle at zero lag is found */
                    PPP[i].DhsmlNgbFactor = 1 / (1 + PPP[i].DhsmlNgbFactor);
                else
                    PPP[i].DhsmlNgbFactor = 1;
                P[i].Particle_DivVel *= PPP[i].DhsmlNgbFactor;
                
                /* now check whether we have enough neighbours */
                redo_particle = 0;
                
                double minsoft = ags_return_minsoft(i);
                double maxsoft = ags_return_maxsoft(i);
                minsoft = DMAX(minsoft , AGS_Prev[i]*AGS_DSOFT_TOL);
                maxsoft = DMIN(maxsoft , AGS_Prev[i]/AGS_DSOFT_TOL);
                if(All.Time==All.TimeBegin)
                {
                    minsoft = All.ForceSoftening[P[i].Type];
                    maxsoft = ADAPTIVE_GRAVSOFT_FORALL * All.ForceSoftening[P[i].Type];
                }
                desnumngb = All.AGS_DesNumNgb;
                desnumngbdev = All.AGS_MaxNumNgbDeviation;
                if(All.Time==All.TimeBegin) {if(All.AGS_MaxNumNgbDeviation > 0.05) desnumngbdev=0.05;}
                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
                double ConditionNumber = do_cbe_nvt_inversion_for_faces(i); // right now we don't do anything with this, but could use to force expansion of search, as in hydro
                if(ConditionNumber > MAX_REAL_NUMBER) {printf("CNUM warning for CBE: ThisTask=%d i=%d ConditionNumber=%g desnumngb=%g NumNgb=%g iter=%d NVT=%g/%g/%g/%g/%g/%g AGS_Hsml=%g \n",ThisTask,i,ConditionNumber,desnumngb,PPP[i].NumNgb,iter,P[i].NV_T[0][0],P[i].NV_T[1][1],P[i].NV_T[2][2],P[i].NV_T[0][1],P[i].NV_T[0][2],P[i].NV_T[1][2],PPP[i].AGS_Hsml);}
                if(iter > 10) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*((double)iter - 9.)) );}
#else
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}
#endif
                
                
                /* check if we are in the 'normal' range between the max/min allowed values */
                if((PPP[i].NumNgb < (desnumngb - desnumngbdev) && PPP[i].AGS_Hsml < 0.99*maxsoft) ||
                   (PPP[i].NumNgb > (desnumngb + desnumngbdev) && PPP[i].AGS_Hsml > 1.01*minsoft))
                    redo_particle = 1;
                
                /* check maximum kernel size allowed */
                particle_set_to_maxhsml_flag = 0;
                if((PPP[i].AGS_Hsml >= 0.99*maxsoft) && (PPP[i].NumNgb < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].AGS_Hsml == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        if(P[i].Type==0) redo_particle = 1;
                        PPP[i].AGS_Hsml = maxsoft;
                        particle_set_to_maxhsml_flag = 1;
                    }
                }
                
                /* check minimum kernel size allowed */
                particle_set_to_minhsml_flag = 0;
                if((PPP[i].AGS_Hsml <= 1.01*minsoft) && (PPP[i].NumNgb > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(PPP[i].AGS_Hsml == minsoft)
                    {
                        /* this means we've already done an iteration with the MinHsml value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minhsml_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        if(P[i].Type==0) redo_particle = 1;
                        PPP[i].AGS_Hsml = minsoft;
                        particle_set_to_minhsml_flag = 1;
                    }
                }
                
                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
#ifndef IO_REDUCED_MODE
                        printf("AGS: i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)\n",
                               i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, PPP[i].AGS_Hsml, PPP[i].DhsmlNgbFactor, Left[i], Right[i],
                               (float) PPP[i].NumNgb, Right[i] - Left[i], particle_set_to_maxhsml_flag, particle_set_to_minhsml_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                        fflush(stdout);
#endif
                    }
                    
                    /* need to redo this particle */
                    npleft++;
                    
                    if(Left[i] > 0 && Right[i] > 0)
                        if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                            continue;
                        }
                    
                    if((particle_set_to_maxhsml_flag==0)&&(particle_set_to_minhsml_flag==0))
                    {
                        if(PPP[i].NumNgb < (desnumngb - desnumngbdev))
                            Left[i] = DMAX(PPP[i].AGS_Hsml, Left[i]);
                        else
                        {
                            if(Right[i] != 0)
                            {
                                if(PPP[i].AGS_Hsml < Right[i])
                                    Right[i] = PPP[i].AGS_Hsml;
                            }
                            else
                                Right[i] = PPP[i].AGS_Hsml;
                        }
                        
                        // right/left define upper/lower bounds from previous iterations
                        if(Right[i] > 0 && Left[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(Right[i]/Left[i]);}
                            if(PPP[i].NumNgb > 1)
                            {
                                double jumpvar = PPP[i].DhsmlNgbFactor * log( desnumngb / PPP[i].NumNgb ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                PPP[i].AGS_Hsml *= exp(jumpvar);
                            } else {
                                PPP[i].AGS_Hsml *= 2.0;
                            }
                            if((PPP[i].AGS_Hsml<Right[i])&&(PPP[i].AGS_Hsml>Left[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(PPP[i].AGS_Hsml > Right[i] / hfac) {PPP[i].AGS_Hsml = Right[i] / hfac;}
                                    if(PPP[i].AGS_Hsml < Left[i] * hfac) {PPP[i].AGS_Hsml = Left[i] * hfac;}
                                }
                            } else {
                                if(PPP[i].AGS_Hsml>Right[i]) PPP[i].AGS_Hsml=Right[i];
                                if(PPP[i].AGS_Hsml<Left[i]) PPP[i].AGS_Hsml=Left[i];
                                PPP[i].AGS_Hsml = pow(PPP[i].AGS_Hsml * Left[i] * Right[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(Right[i] == 0 && Left[i] == 0)
                            {
                                char buf[1000];
                                sprintf(buf, "AGS: Right[i] == 0 && Left[i] == 0 && PPP[i].AGS_Hsml=%g\n", PPP[i].AGS_Hsml);
                                terminate(buf);
                            }
                            
                            if(Right[i] == 0 && Left[i] > 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                                
                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=10)
                                        if(PPP[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                    
                                    if(fac < fac_lim+0.231)
                                    {
                                        PPP[i].AGS_Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        PPP[i].AGS_Hsml *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    PPP[i].AGS_Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                            
                            if(Right[i] > 0 && Left[i] == 0)
                            {
                                if (PPP[i].NumNgb > 1)
                                    fac_lim = log( desnumngb / PPP[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
                                
                                if (fac_lim < -1.535) fac_lim = -1.535; // decreasing N_ngb by factor ~100
                                
                                if((PPP[i].NumNgb < 2*desnumngb)&&(PPP[i].NumNgb > 0.1*desnumngb))
                                {
                                    double slope = PPP[i].DhsmlNgbFactor;
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=10)
                                        if(PPP[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
                                    
                                    if(fac > fac_lim-0.231)
                                    {
                                        PPP[i].AGS_Hsml *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        PPP[i].AGS_Hsml *= exp(fac_lim-0.231); // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    PPP[i].AGS_Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if[particle_set_to_max/minhsml_flag]
                    } // closes redo_particle
                    /* resets for max/min values */
                    if(PPP[i].AGS_Hsml < minsoft) PPP[i].AGS_Hsml = minsoft;
                    if(particle_set_to_minhsml_flag==1) PPP[i].AGS_Hsml = minsoft;
                    if(PPP[i].AGS_Hsml > maxsoft) PPP[i].AGS_Hsml = maxsoft;
                    if(particle_set_to_maxhsml_flag==1) PPP[i].AGS_Hsml = maxsoft;
                }
                else
                    P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
            } //  if(ags_density_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0)
            {
                printf("ags-ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
            }
            if(iter > MAXITER)
            {
                printf("ags-failed to converge in neighbour iteration in density()\n");
                fflush(stdout);
                endrun(1155);
            }
        }
    }
    while(ntot > 0);

    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Right);
    myfree(Left);
    myfree(Ngblist);
    
    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0)
            P[i].TimeBin = -P[i].TimeBin - 1;
    }

    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(ags_density_isactive(i))
        {
            if((P[i].Mass>0)&&(PPP[i].AGS_Hsml>0)&&(PPP[i].NumNgb>0))
            {
                double minsoft = ags_return_minsoft(i);
                double maxsoft = ags_return_maxsoft(i);
                minsoft = DMAX(minsoft , AGS_Prev[i]*AGS_DSOFT_TOL);
                maxsoft = DMIN(maxsoft , AGS_Prev[i]/AGS_DSOFT_TOL);
                /* check that we're within the 'valid' range for adaptive softening terms, otherwise zeta=0 */
                if((fabs(PPP[i].NumNgb-All.AGS_DesNumNgb)/All.AGS_DesNumNgb < 0.05)
                   &&(PPP[i].AGS_Hsml <= 0.99*maxsoft)&&(PPP[i].AGS_Hsml >= 1.01*minsoft)
                   &&(PPP[i].NumNgb >= (All.AGS_DesNumNgb - All.AGS_MaxNumNgbDeviation))
                   &&(PPP[i].NumNgb <= (All.AGS_DesNumNgb + All.AGS_MaxNumNgbDeviation)))
                {
                    double ndenNGB = PPP[i].NumNgb / ( NORM_COEFF * pow(PPP[i].AGS_Hsml,NUMDIMS) );
                    PPPZ[i].AGS_zeta *= 0.5 * P[i].Mass * PPP[i].AGS_Hsml / (NUMDIMS * ndenNGB) * PPP[i].DhsmlNgbFactor;
                } else {
                    PPPZ[i].AGS_zeta = 0;
                }
                /* convert NGB to the more useful format, NumNgb^(1/NDIMS), which we can use to obtain the corrected particle sizes */
                PPP[i].NumNgb = pow(PPP[i].NumNgb , 1./NUMDIMS);
            } else {
                PPPZ[i].AGS_zeta = 0;
                PPP[i].NumNgb = 0;
            }
#ifdef PM_HIRES_REGION_CLIPPING
            if(PPP[i].NumNgb <= 0) {P[i].Mass = 0;}
            if((PPP[i].AGS_Hsml <= 0) || (PPP[i].AGS_Hsml >= PM_HIRES_REGION_CLIPPING)) {P[i].Mass = 0;}
            double vmag=0; for(k=0;k<3;k++) {vmag+=P[i].Vel[k]*P[i].Vel[k];} vmag = sqrt(vmag);
            if(vmag>5.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) {P[i].Mass=0;}
            if(vmag>1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) {for(k=0;k<3;k++) {P[i].Vel[k]*=(1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s)/vmag;}}
#endif
        }
    }
    myfree(AGS_Prev);
    
    /* collect some timing information */
    
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm;
    CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}






/*! This function represents the core of the density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int ags_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, n;
    int startnode, numngb_inbox, listindex = 0;
    double r2, h2, u;
    struct kernel_density kernel;
    struct ags_densdata_in local;
    struct ags_densdata_out out;
    memset(&out, 0, sizeof(struct ags_densdata_out));
    
    if(mode == 0)
        ags_particle2in_density(&local, target);
    else
        local = AGS_DensDataGet[target];
    
    h2 = local.AGS_Hsml * local.AGS_Hsml;
    kernel_hinv(local.AGS_Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms
    
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = AGS_DensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
    
    
    double fac_mu = -3 / (All.cf_afac3 * All.cf_atime);
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.AGS_Hsml, target, &startnode, mode, exportflag,
                                          exportnodecount, exportindex, ngblist, AGS_kernel_shared_BITFLAG);
            
            if(numngb_inbox < 0)
                return -1;
            
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Mass <= 0) continue;
                
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
                if(r2 < h2)
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);

                    out.Ngb += kernel.wk;
                    out.DhsmlNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);

                    if(kernel.r > 0)
                    {
                        out.AGS_zeta += P[j].Mass * kernel_gravity(u, kernel.hinv, kernel.hinv3, 0);
                        if(P[j].Type==0)
                        {
                            kernel.dv[0] = local.Vel[0] - SphP[j].VelPred[0];
                            kernel.dv[1] = local.Vel[1] - SphP[j].VelPred[1];
                            kernel.dv[2] = local.Vel[2] - SphP[j].VelPred[2];
                        } else {
                            kernel.dv[0] = local.Vel[0] - P[j].Vel[0];
                            kernel.dv[1] = local.Vel[1] - P[j].Vel[1];
                            kernel.dv[2] = local.Vel[2] - P[j].Vel[2];
                        }
#ifdef BOX_SHEARING
                        if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
                        if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
#endif
                        double v_dot_r = kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2];
                        double vsig = 0.5 * fabs( fac_mu * v_dot_r / kernel.r );
                        if(TimeBinActive[P[j].TimeBin]) {if(vsig > PPP[j].AGS_vsig) PPP[j].AGS_vsig = vsig;}
                        if(vsig > out.AGS_vsig) {out.AGS_vsig = vsig;}
#ifdef WAKEUP
                        if(!(TimeBinActive[P[j].TimeBin]) && (All.Time > All.TimeBegin)) {if(vsig > WAKEUP*P[j].AGS_vsig) {P[j].wakeup = 1;}}
#if defined(GALSF)
                        if((P[j].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[j].Type == 2)||(P[j].Type==3)))) {P[j].wakeup = 0;} // don't wakeup star particles, or risk 2x-counting feedback events! //
#endif
#endif
                        out.Particle_DivVel -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve */
                        
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
                        out.NV_T[0][0] +=  kernel.wk * kernel.dp[0] * kernel.dp[0];
                        out.NV_T[0][1] +=  kernel.wk * kernel.dp[0] * kernel.dp[1];
                        out.NV_T[0][2] +=  kernel.wk * kernel.dp[0] * kernel.dp[2];
                        out.NV_T[1][1] +=  kernel.wk * kernel.dp[1] * kernel.dp[1];
                        out.NV_T[1][2] +=  kernel.wk * kernel.dp[1] * kernel.dp[2];
                        out.NV_T[2][2] +=  kernel.wk * kernel.dp[2] * kernel.dp[2];
#endif
                    }
                }
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = AGS_DensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
        ags_out2particle_density(&out, target, 0);
    else
        AGS_DensDataResult[target] = out;
    
    return 0;
}



void *ags_density_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(ags_density_isactive(i))
#define EVALUATION_CALL ags_density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *ags_density_evaluate_secondary(void *p)
{
#define EVALUATION_CALL ags_density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}


/* routine to determine if we need to use ags_density to calculate Hsml */
int ags_density_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0; /* check our 'marker' for particles which have finished
                                        iterating to an Hsml solution (if they have, dont do them again) */
    if(P[i].Type==0)
    {
        PPP[i].AGS_Hsml = PPP[i].Hsml; // gas sees gas, these are identical
        return 0; // don't actually need to do the loop //
    }
    return 1;
}

/* routine to return the maximum allowed softening */
double ags_return_maxsoft(int i)
{
    double maxsoft = All.MaxHsml; // overall maximum - nothing is allowed to exceed this
#if !(EXPAND_PREPROCESSOR_(ADAPTIVE_GRAVSOFT_FORALL) == 1)
    maxsoft = DMIN(maxsoft, ADAPTIVE_GRAVSOFT_FORALL * All.ForceSoftening[P[i].Type]); // user-specified maximum
#ifdef PMGRID
    /*!< this gives the maximum allowed gravitational softening when using the TreePM method.
     *  The quantity is given in units of the scale used for the force split (ASMTH) */
    maxsoft = DMIN(maxsoft, ADAPTIVE_GRAVSOFT_FORALL * 0.5 * All.Asmth[0]); /* no more than 1/2 the size of the largest PM cell */
#endif
#else
    maxsoft = DMIN(maxsoft, 50.0 * All.ForceSoftening[P[i].Type]);
#ifdef PMGRID
    maxsoft = DMIN(maxsoft, 0.5 * All.Asmth[0]); /* no more than 1/2 the size of the largest PM cell */
#endif
#endif
#ifdef BLACK_HOLES
    if(P[i].Type == 5) {maxsoft = All.BlackHoleMaxAccretionRadius  / All.cf_atime;}   // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#endif
    return maxsoft;
}

/* routine to return the minimum allowed softening */
double ags_return_minsoft(int i)
{
    return All.ForceSoftening[P[i].Type]; // this is the user-specified minimum
}


/* routine to return effective particle sizes (inter-particle separation) based on AGS_Hsml saved values */
double INLINE_FUNC Get_Particle_Size_AGS(int i)
{
    /* in previous versions of the code, we took NumNgb^(1/NDIMS) here; however, now we
     take that when NumNgb is computed (at the end of the density routine), so we
     don't have to re-compute it each time. That makes this function fast enough to
     call -inside- of loops (e.g. hydro computations) */
#if (NUMDIMS == 1)
    return 2.00000 * PPP[i].AGS_Hsml / PPP[i].NumNgb;
#endif
#if (NUMDIMS == 2)
    return 1.25331 * PPP[i].AGS_Hsml / PPP[i].NumNgb; // sqrt(Pi/2)
#endif
#if (NUMDIMS == 3)
    return 1.61199 * PPP[i].AGS_Hsml / PPP[i].NumNgb; // (4pi/3)^(1/3)
#endif
}


/* --------------------------------------------------------------------------
 very quick sub-routine to get the particle densities from their volumes
 -------------------------------------------------------------------------- */
double get_particle_volume_ags(int j)
{
    double L_j = Get_Particle_Size_AGS(j);
#if (NUMDIMS==1)
    return L_j;
#elif (NUMDIMS==2)
    return L_j*L_j;
#else
    return L_j*L_j*L_j;
#endif
}


#ifdef AGS_FACE_CALCULATION_IS_ACTIVE

/* --------------------------------------------------------------------------
 Subroutine here exists to calculate the MFM-like effective faces for purposes of face-interaction evaluation
 -------------------------------------------------------------------------- */

/* routine to invert the NV_T matrix after neighbor pass */
double do_cbe_nvt_inversion_for_faces(int i)
{
    MyFloat NV_T[3][3]; int j,k;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k]=P[i].NV_T[j][k];}} // initialize matrix to be inverted //
    double Tinv[3][3], FrobNorm=0, FrobNorm_inv=0, detT=0;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {Tinv[j][k]=0;}}
    /* fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
    NV_T[1][0]=NV_T[0][1]; NV_T[2][0]=NV_T[0][2]; NV_T[2][1]=NV_T[1][2];
    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
     this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {FrobNorm += NV_T[j][k]*NV_T[j][k];}}
#if (NUMDIMS==1) // 1-D case //
    detT = NV_T[0][0];
    if(detT!=0 && !isnan(detT)) {Tinv[0][0] = 1/detT}; /* only one non-trivial element in 1D! */
#endif
#if (NUMDIMS==2) // 2-D case //
    detT = NV_T[0][0]*NV_T[1][1] - NV_T[0][1]*NV_T[1][0];
    if((detT != 0)&&(!isnan(detT)))
    {
        Tinv[0][0] = NV_T[1][1] / detT; Tinv[0][1] = -NV_T[0][1] / detT;
        Tinv[1][0] = -NV_T[1][0] / detT; Tinv[1][1] = NV_T[0][0] / detT;
    }
#endif
#if (NUMDIMS==3) // 3-D case //
    detT = NV_T[0][0] * NV_T[1][1] * NV_T[2][2] + NV_T[0][1] * NV_T[1][2] * NV_T[2][0] +
    NV_T[0][2] * NV_T[1][0] * NV_T[2][1] - NV_T[0][2] * NV_T[1][1] * NV_T[2][0] -
    NV_T[0][1] * NV_T[1][0] * NV_T[2][2] - NV_T[0][0] * NV_T[1][2] * NV_T[2][1];
    /* check for zero determinant */
    if((detT != 0) && !isnan(detT))
    {
        Tinv[0][0] = (NV_T[1][1] * NV_T[2][2] - NV_T[1][2] * NV_T[2][1]) / detT;
        Tinv[0][1] = (NV_T[0][2] * NV_T[2][1] - NV_T[0][1] * NV_T[2][2]) / detT;
        Tinv[0][2] = (NV_T[0][1] * NV_T[1][2] - NV_T[0][2] * NV_T[1][1]) / detT;
        Tinv[1][0] = (NV_T[1][2] * NV_T[2][0] - NV_T[1][0] * NV_T[2][2]) / detT;
        Tinv[1][1] = (NV_T[0][0] * NV_T[2][2] - NV_T[0][2] * NV_T[2][0]) / detT;
        Tinv[1][2] = (NV_T[0][2] * NV_T[1][0] - NV_T[0][0] * NV_T[1][2]) / detT;
        Tinv[2][0] = (NV_T[1][0] * NV_T[2][1] - NV_T[1][1] * NV_T[2][0]) / detT;
        Tinv[2][1] = (NV_T[0][1] * NV_T[2][0] - NV_T[0][0] * NV_T[2][1]) / detT;
        Tinv[2][2] = (NV_T[0][0] * NV_T[1][1] - NV_T[0][1] * NV_T[1][0]) / detT;
    }
#endif
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {FrobNorm_inv += Tinv[j][k]*Tinv[j][k];}}
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {P[i].NV_T[j][k]=Tinv[j][k];}} // now P[i].NV_T holds the inverted matrix elements //
    double ConditionNumber = DMAX(sqrt(FrobNorm * FrobNorm_inv) / NUMDIMS, 1); // = sqrt( ||NV_T^-1||*||NV_T|| ) :: should be ~1 for a well-conditioned matrix //
#ifdef CBE_DEBUG
    if((ThisTask==0)&&(ConditionNumber>1.0e10)) {printf("Condition number == %g (Task=%d i=%d)\n",ConditionNumber,ThisTask,i);}
#endif
    return ConditionNumber;
}

#endif





/* ------------------------------------------------------------------------------------------------------
 Everything below here is a giant block to define the sub-routines needed to calculate additional force
  terms for particle types that do not fall into the 'hydro' category.
 -------------------------------------------------------------------------------------------------------- */

/* structure for variables needed in evaluation sub-routines which must be passed from particles (sent to other processors) */
struct AGSForce_data_in
{
    double Mass;
    double AGS_Hsml;
    double Pos[3];
    double Vel[3];
    int NodeList[NODELISTLENGTH];
    int Type;
    int dt_step;
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    double NV_T[3][3];
    double V_i;
#endif
#ifdef DM_SIDM
    int dt_step_sidm;
    MyIDType ID;
#endif
}
*AGSForce_DataIn, *AGSForce_DataGet;

/* structure for variables which must be returned -from- the evaluation sub-routines */
struct AGSForce_data_out
{
#ifdef DM_SIDM
    double sidm_kick[3];
    int dt_step_sidm;
    int si_count;
#endif
}
*AGSForce_DataResult, *AGSForce_DataOut;

/* this is a temporary structure for quantities used ONLY in the loop below, for example for computing the slope-limiters (for the Reimann problem) */
/* (currently not used, but can un-comment this and add as needed, if needed for future sub-routines)
static struct temporary_data_topass
{}
*AGSForce_DataPasser;
*/

struct kernel_AGSForce
{
    double dp[3], dv[3], r, wk_i, wk_j, dwk_i, dwk_j, h_i, hinv_i, hinv3_i, hinv4_i, h_j, hinv_j, hinv3_j, hinv4_j;
};

/* routine to pass particle information to the actual evaluation sub-routines */
static inline void particle2in_AGSForce(struct AGSForce_data_in *in, int i);
static inline void particle2in_AGSForce(struct AGSForce_data_in *in, int i)
{
    in->Mass = PPP[i].Mass;
    in->AGS_Hsml = PPP[i].AGS_Hsml;
    in->Type = P[i].Type;
    in->dt_step = P[i].dt_step;
    int k,k2; k=0; k2=0;
    for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    for(k=0;k<3;k++) {in->Vel[k] = P[i].Vel[k];}
#if defined(AGS_FACE_CALCULATION_IS_ACTIVE)
    in->V_i = get_particle_volume_ags(i);
    for(k=0;k<3;k++) {for(k2=0;k2<3;k2++) {in->NV_T[k][k2] = P[i].NV_T[k][k2];}}
#endif
#ifdef DM_SIDM
    in->dt_step_sidm = P[i].dt_step_sidm;
    in->ID = P[i].ID;
#endif
}

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))

static inline void out2particle_AGSForce(struct AGSForce_data_out *out, int i, int mode);
static inline void out2particle_AGSForce(struct AGSForce_data_out *out, int i, int mode)
{
    int k,k2; k=0; k2=0;
#ifdef DM_SIDM
    for(k=0;k<3;k++) {P[i].Vel[k] += out->sidm_kick[k];}
    MIN_ADD(P[i].dt_step_sidm, out->dt_step_sidm, mode);
    P[i].NInteractions += out->si_count;
#endif

}


void AGSForce_calc(void)
{
    int i, j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
    double timecomp, timecomm, timewait, tstart, tend, t0, t1;
    long long n_exported = 0, NTaskTimesNumPart;
    /* allocate buffers to arrange communication */
    //AGSForce_DataPasser = (struct temporary_data_topass *) mymalloc("AGSForce_DataPasser",NumPart * sizeof(struct temporary_data_topass));
    NTaskTimesNumPart = maxThreads * NumPart;
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                           sizeof(struct AGSForce_data_in) + sizeof(struct AGSForce_data_out) + sizemax(sizeof(struct AGSForce_data_in),sizeof(struct AGSForce_data_out))));
    CPU_Step[CPU_AGSDENSMISC] += measure_time();
    t0 = my_second();
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

    /* before doing any operations, need to zero the appropriate memory so we can correctly do pair-wise operations */
    //for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {if(P[i].Type==0) {memset(&AGSForce_DataPasser[i], 0, sizeof(struct temporary_data_topass));}}
#ifdef DM_SIDM
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {P[i].dt_step_sidm = 1.e10*P[i].dt_step;}
#endif

    /* begin the main gradient loop */
    NextParticle = FirstActiveParticle;    /* begin with this index */
    do
    {
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        /* do local particles and prepare export list */
        tstart = my_second();
        
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            AGSForce_evaluate_primary(&mainthreadid);    /* do local particles and prepare export list */
        }
        
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            NextParticle = save_NextParticle;
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle) break;
                if(ProcessedFlag[NextParticle] != 1) break;
                ProcessedFlag[NextParticle] = 2;
                NextParticle = NextActiveParticle[NextParticle];
            }
            if(NextParticle == save_NextParticle)
            {
                endrun(123708); /* in this case, the buffer is too small to process even a single particle */
            }
            int new_export = 0;
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) {k = j + 1;}
                    for(; k < Nexport; k++)
                        if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                        {
                            int old_index = DataIndexTable[j].Index;
                            DataIndexTable[j] = DataIndexTable[k];
                            DataNodeList[j] = DataNodeList[k];
                            DataIndexTable[j].IndexGet = j;
                            new_export++;
                            DataIndexTable[k].Index = old_index;
                            k++;
                            break;
                        }
                }
                else {new_export++;}
            Nexport = new_export;
        }
        n_exported += Nexport;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
        for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
        tstart = my_second();
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        tend = my_second();
        timewait1 += timediff(tstart, tend);
        for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            Nimport += Recv_count[j];
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        
        /* prepare particle data for export */
        AGSForce_DataGet = (struct AGSForce_data_in *) mymalloc("AGSForce_DataGet", Nimport * sizeof(struct AGSForce_data_in));
        AGSForce_DataIn = (struct AGSForce_data_in *) mymalloc("AGSForce_DataIn", Nexport * sizeof(struct AGSForce_data_in));
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_AGSForce(&AGSForce_DataIn[j], place);
            memcpy(AGSForce_DataIn[j].NodeList,DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        
        /* exchange particle data */
        tstart = my_second();
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&AGSForce_DataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct AGSForce_data_in), MPI_BYTE,
                                 recvTask, TAG_GRADLOOP_A,
                                 &AGSForce_DataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct AGSForce_data_in), MPI_BYTE,
                                 recvTask, TAG_GRADLOOP_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm1 += timediff(tstart, tend);
        myfree(AGSForce_DataIn);
        AGSForce_DataResult = (struct AGSForce_data_out *) mymalloc("AGSForce_DataResult", Nimport * sizeof(struct AGSForce_data_out));
        AGSForce_DataOut = (struct AGSForce_data_out *) mymalloc("AGSForce_DataOut", Nexport * sizeof(struct AGSForce_data_out));
        
        /* now do the particles that were sent to us */
        tstart = my_second();
        NextJ = 0;
        
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            AGSForce_evaluate_secondary(&mainthreadid);
        }
        
        tend = my_second();
        timecomp2 += timediff(tstart, tend);
        
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        tstart = my_second();
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        tend = my_second();
        timewait2 += timediff(tstart, tend);
        
        /* get the result */
        tstart = my_second();
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&AGSForce_DataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct AGSForce_data_out),
                                 MPI_BYTE, recvTask, TAG_GRADLOOP_B,
                                 &AGSForce_DataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct AGSForce_data_out),
                                 MPI_BYTE, recvTask, TAG_GRADLOOP_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        tend = my_second();
        timecommsumm2 += timediff(tstart, tend);
        
        /* add the result to the local particles */
        tstart = my_second();
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_AGSForce(&AGSForce_DataOut[j], place, 1);
        }
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        myfree(AGSForce_DataOut);
        myfree(AGSForce_DataResult);
        myfree(AGSForce_DataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
    
    /* do final operations on results: these are operations that can be done after the complete set of iterations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
    }
    //myfree(AGSForce_DataPasser); /* free the temporary structure we created for the MinMax and additional data passing */
    MPI_Barrier(MPI_COMM_WORLD); // force barrier so we know the first derivatives are fully-computed //

    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm;
    CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}



int AGSForce_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    /* define variables */
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double r2, u_i, u_j;
    struct kernel_AGSForce kernel;
    struct AGSForce_data_in local;
    struct AGSForce_data_out out;
    /* zero memory and import data for local target */
    memset(&out, 0, sizeof(struct AGSForce_data_out)); memset(&kernel, 0, sizeof(struct kernel_AGSForce));
    if(mode == 0) {particle2in_AGSForce(&local, target);} else {local = AGSForce_DataGet[target];}
    /* check if we should bother doing a neighbor loop */
    if(local.Mass <= 0) return 0;
    if(local.AGS_Hsml <= 0) return 0;
    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.AGS_Hsml;
    kernel_hinv(kernel.h_i, &kernel.hinv_i, &kernel.hinv3_i, &kernel.hinv4_i);
    int AGS_kernel_shared_BITFLAG = ags_gravity_kernel_shared_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms
#ifdef DM_SIDM
    out.dt_step_sidm = local.dt_step_sidm;
#endif

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = AGSForce_DataGet[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            double search_len = local.AGS_Hsml;
#ifdef DM_SIDM
            search_len *= 3.0; // need a 'buffer' because we will consider interactions with any kernel -overlap, not just inside one or the other kernel radius
#endif
            numngb_inbox = ngb_treefind_pairs_threads_targeted(local.Pos, search_len, target, &startnode, mode, exportflag,
                                                                  exportnodecount, exportindex, ngblist, AGS_kernel_shared_BITFLAG);
            if(numngb_inbox < 0) {return -1;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n];
                if((P[j].Mass <= 0)||(PPP[j].AGS_Hsml <= 0)) continue; /* make sure neighbor is valid */
                /* calculate position relative to target */
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0]; kernel.dp[1] = local.Pos[1] - P[j].Pos[1]; kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); /*  now find the closest image in the given box size  */
#endif
                r2 = kernel.dp[0]*kernel.dp[0] + kernel.dp[1]*kernel.dp[1] + kernel.dp[2]*kernel.dp[2];
                if(r2 <= 0) continue;
                kernel.r = sqrt(r2);
                kernel.h_j = PPP[j].AGS_Hsml;
#ifdef DM_SIDM
                if(kernel.r > kernel.h_i+kernel.h_j) continue;
#else
                if(kernel.r > kernel.h_i && kernel.r > kernel.h_j) continue;
#endif
                /* calculate kernel quantities needed below */
                kernel_hinv(kernel.h_j, &kernel.hinv_j, &kernel.hinv3_j, &kernel.hinv4_j);
                u_i = kernel.r * kernel.hinv_i; u_j = kernel.r * kernel.hinv_j;
                if(u_i < 1) {kernel_main(u_i, kernel.hinv3_i, kernel.hinv4_i, &kernel.wk_i, &kernel.dwk_i, 0);} else {kernel.wk_i=kernel.dwk_i=0;}
                if(u_j < 1) {kernel_main(u_j, kernel.hinv3_j, kernel.hinv4_j, &kernel.wk_j, &kernel.dwk_j, 0);} else {kernel.wk_j=kernel.dwk_j=0;}
                for(k=0;k<3;k++)
                {
                    kernel.dv[k] = local.Vel[k] - P[j].Vel[k];
                    if(All.ComovingIntegrationOn) {kernel.dv[k] += All.cf_hubble_a * kernel.dp[k]/All.cf_a2inv;}
                }
                
#ifdef DM_SIDM
#include "../sidm/sidm_core_flux_computation.h"
#endif

            } // numngb_inbox loop
        } // while(startnode)
        /* continue to open leaves if needed */
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = AGSForce_DataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
            }
        }
    }
    
    /* Collect the result at the right place */
    if(mode == 0) {out2particle_AGSForce(&out, target, 0);} else {AGSForce_DataResult[target] = out;}
    return 0;
}

void *AGSForce_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(AGSForce_isactive(i))
#define EVALUATION_CALL AGSForce_evaluate(i,0,exportflag,exportnodecount,exportindex,ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}

void *AGSForce_evaluate_secondary(void *p)
{
#define EVALUATION_CALL AGSForce_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}


/* routine to determine if we need to apply the additional AGS-Force calculation[s] */
int AGSForce_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0; /* check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
#ifdef DM_SIDM
    if((1 << P[i].Type) & (DM_SIDM)) return 1;
#endif
    return 0; // default to no-action, need to affirm calculation above //
}





#endif // ADAPTIVE_GRAVSOFT_FORALL
