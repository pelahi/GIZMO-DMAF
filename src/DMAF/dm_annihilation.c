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

/*! \file dm_annihilation
 *  \brief DM density calculation around particles
 *
 *  This file contains a loop modeled on the standard gas density computation which
 *    determines the DARK MATTER density around a given set of particles and adjusts the smoothing length for this
*     calculation.
 *    The dark matter density is used to set the energy injection due to dark matter annihilation. 
 *
 * This file was written by Florian List, for GIZMO, based on Phil Hopkins's adaptive gravitational softening
 *    routine.
 *
 * Method 1 (activate using DMANNIHILATION)
  - DM density is calculated at each gas particle using smoothing length HsmlDM
  - from DM density, DM annihilation rate at gas particle is calculated
  - energy injection at each gas particle
 * Method 2 also uses the functionality in this file in order to determine the DM density around DM particles.
 */

#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)

struct kernel_density
{
  double dp[3],dv[3],r;
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct dm_densdata_in
{
    MyDouble Pos[3];
	MyFloat Vel[3];
    MyFloat HsmlDM;
    int NodeList[NODELISTLENGTH];
}
*DM_DensDataIn, *DM_DensDataGet;

static struct dm_densdata_out
{
    MyLongDouble NgbDM, RhoDM, DhsmlNgbDM, Particle_DivVelDM;
}
*DM_DensDataResult, *DM_DensDataOut;

void dm_particle2in_density(struct dm_densdata_in *in, int i);
void dm_out2particle_density(struct dm_densdata_out *out, int i, int mode);

void dm_particle2in_density(struct dm_densdata_in *in, int i)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        if(P[i].Type==0) {in->Vel[k]=SphP[i].VelPred[k];} else {in->Vel[k]=P[i].Vel[k];}
    }
    in->HsmlDM = P[i].HsmlDM;
}

void dm_out2particle_density(struct dm_densdata_out *out, int i, int mode)
{
	ASSIGN_ADD(P[i].NumNgbDM, out->NgbDM, mode);
	ASSIGN_ADD(P[i].DensityDM, out->RhoDM, mode);
	ASSIGN_ADD(P[i].DhsmlNgbFactorDM, out->DhsmlNgbDM, mode);
	ASSIGN_ADD(P[i].Particle_DivVelDM, out->Particle_DivVelDM, mode); 
}

void dm_density(void)
{
    MyFloat *LeftDM, *RightDM;
    int i, j, k, ndone, ndone_flag, npleft, iter = 0;
    int ngrp, recvTask, place;
    long long ntot;
    double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 = 0, timewait2 = 0;
    double timecomp, timecomm, timewait;
    double tstart, tend, t0, t1;
    double desnumngb, desnumngbdev;
    int save_NextParticle;
    long long n_exported = 0;
    int redo_particle;
    int particle_set_to_minhsml_flag_DM = 0;
    int particle_set_to_maxhsml_flag_DM = 0;
	double fac, fac_lim;
    
    CPU_Step[CPU_DMDENSCOMPUTE] += measure_time();
    
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    
    LeftDM = (MyFloat *) mymalloc("LeftDM", NumPart * sizeof(MyFloat));
    RightDM = (MyFloat *) mymalloc("RightDM", NumPart * sizeof(MyFloat));
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(dm_density_isactive(i))
        {
            P[i].NumNgbDM = 0;
            LeftDM[i] = RightDM[i] = 0;
        }
    }
    
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                           sizeof(struct dm_densdata_in) + sizeof(struct dm_densdata_out) +
                                                           sizemax(sizeof(struct dm_densdata_in),sizeof(struct dm_densdata_out))));
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
                pthread_create(&mythreads[j], &attr, dm_density_evaluate_primary, &threadid[j]);
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
                dm_density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
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
                    printf("DMAF: Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
                           P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);
                    
                    endrun(290494);
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
            
            DM_DensDataGet = (struct dm_densdata_in *) mymalloc("DM_DensDataGet", Nimport * sizeof(struct dm_densdata_in));
            DM_DensDataIn = (struct dm_densdata_in *) mymalloc("DM_DensDataIn", Nexport * sizeof(struct dm_densdata_in));
            
            /* prepare particle data for export */
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                
                dm_particle2in_density(&DM_DensDataIn[j], place);
                
                memcpy(DM_DensDataIn[j].NodeList,
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
                        MPI_Sendrecv(&DM_DensDataIn[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct dm_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DMDENS_A,
                                     &DM_DensDataGet[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct dm_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DMDENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm1 += timediff(tstart, tend);
            
            myfree(DM_DensDataIn);
            DM_DensDataResult = (struct dm_densdata_out *) mymalloc("DM_DensDataResult", Nimport * sizeof(struct dm_densdata_out));
            DM_DensDataOut = (struct dm_densdata_out *) mymalloc("DM_DensDataOut", Nexport * sizeof(struct dm_densdata_out));
            
            /* now do the particles that were sent to us */
            
            tstart = my_second();
            
            NextJ = 0;
            
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                pthread_create(&mythreads[j], &attr, dm_density_evaluate_secondary, &threadid[j]);
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
                dm_density_evaluate_secondary(&mainthreadid);
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
                        MPI_Sendrecv(&DM_DensDataResult[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct dm_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DMDENS_B,
                                     &DM_DensDataOut[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct dm_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DMDENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
                dm_out2particle_density(&DM_DensDataOut[j], place, 1);
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            
            myfree(DM_DensDataOut);
            myfree(DM_DensDataResult);
            myfree(DM_DensDataGet);
        }
        while(ndone < NTask);
        
		/* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
	        if(dm_density_isactive(i)) // This makes sure that for method 1 (2), only gas particles (DM particles) are treated
	        {
	            if(P[i].NumNgbDM > 0)
	            {
	                P[i].DhsmlNgbFactorDM *= P[i].HsmlDM / (NUMDIMS * P[i].NumNgbDM);
	                P[i].Particle_DivVelDM /= P[i].NumNgbDM;
	                /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
	                P[i].NumNgbDM *= NORM_COEFF * pow(P[i].HsmlDM,NUMDIMS);
	            } else {
	                P[i].NumNgbDM = P[i].DhsmlNgbFactorDM = P[i].Particle_DivVelDM = 0;
	            }
	            
	            // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
	            if(P[i].DhsmlNgbFactorDM > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
	                P[i].DhsmlNgbFactorDM = 1 / (1 + P[i].DhsmlNgbFactorDM);
	            else
	                P[i].DhsmlNgbFactorDM = 1;
	            P[i].Particle_DivVelDM *= P[i].DhsmlNgbFactorDM;		        		       		            
	            
				double minsoft = All.MinHsml;
            	double maxsoft = All.MaxHsml;
			    desnumngb = All.DesNumNgb; // Leave it constant
				desnumngbdev = All.MaxNumNgbDeviation;

				/* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}

				redo_particle = 0; // set redo_particle = 0, check if it needs to be set to 1 in the following

	            /* check if we are in the 'normal' range between the max/min allowed values */
	            if((P[i].NumNgbDM < (desnumngb - desnumngbdev) && P[i].HsmlDM < 0.99*maxsoft) ||
	               (P[i].NumNgbDM > (desnumngb + desnumngbdev) && P[i].HsmlDM > 1.01*minsoft))
	                redo_particle = 1;
	            
	            /* check maximum kernel size allowed */
	            particle_set_to_maxhsml_flag_DM = 0;
	            if((P[i].HsmlDM >= 0.99*maxsoft) && (P[i].NumNgbDM < (desnumngb - desnumngbdev)))
	            {
	                redo_particle = 0;
	                if(P[i].HsmlDM == maxsoft)
	                {
	                    /* iteration at the maximum value is already complete */
	                    particle_set_to_maxhsml_flag_DM = 0;
	                } else {
	                    /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
	                    if(P[i].Type==0) redo_particle = 1;
	                    P[i].HsmlDM = maxsoft;
	                    particle_set_to_maxhsml_flag_DM = 1;
	                }
	            }
	            
	            /* check minimum kernel size allowed */
	            particle_set_to_minhsml_flag_DM = 0;
	            if((P[i].HsmlDM <= 1.01*minsoft) && (P[i].NumNgbDM > (desnumngb + desnumngbdev)))
	            {
	                redo_particle = 0;
	                if(P[i].HsmlDM == minsoft)
	                {
	                    /* this means we've already done an iteration with the MinHsml value, so the
	                     neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
	                    particle_set_to_minhsml_flag_DM = 0;
	                } else {
	                    /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
	                    if(P[i].Type==0) redo_particle = 1;
	                    P[i].HsmlDM = minsoft;
	                    particle_set_to_minhsml_flag_DM = 1;
	                }
	            }		               		            
	            
	            if(redo_particle)
	            {
	                if(iter >= MAXITER - 10)
	                {
	                    printf("DMAF loop parameters:\n i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)\n",
	                           i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, P[i].HsmlDM, P[i].DhsmlNgbFactorDM, LeftDM[i], RightDM[i],
	                           (float) P[i].NumNgbDM, RightDM[i] - LeftDM[i], particle_set_to_maxhsml_flag_DM, particle_set_to_minhsml_flag_DM, minsoft,
	                           maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);							
						printf("SLOW CONVERGENCE IN DMAF LOOP!\n");
	                }
	                
	                /* need to redo this particle */
	                npleft++;
	                
	                if(LeftDM[i] > 0 && RightDM[i] > 0)
	                    if((RightDM[i] - LeftDM[i]) < 1.0e-3 * LeftDM[i])
	                    {
	                        /* this one should be ok */
	                        npleft--;
	                        P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
							redo_particle = 0;
	                        continue;
	                    }
	                
	                if((particle_set_to_maxhsml_flag_DM==0)&&(particle_set_to_minhsml_flag_DM==0))
	                {
	                    if(P[i].NumNgbDM < (desnumngb - desnumngbdev))
	                        LeftDM[i] = DMAX(P[i].HsmlDM, LeftDM[i]);
	                    else
	                    {
	                        if(RightDM[i] != 0)
	                        {
	                            if(P[i].HsmlDM < RightDM[i])
	                                RightDM[i] = P[i].HsmlDM;
	                        }
	                        else
	                            RightDM[i] = P[i].HsmlDM;
	                    }
	                    
	                    // right/left define upper/lower bounds from previous iterations
	                    if(RightDM[i] > 0 && LeftDM[i] > 0)
	                    {
	                        // geometric interpolation between right/left //
	                        double maxjump=0;
	                        if(iter>1) {maxjump = 0.2*log(RightDM[i]/LeftDM[i]);}
	                        if(P[i].NumNgbDM > 1)
	                        {
	                            double jumpvar = P[i].DhsmlNgbFactorDM * log( desnumngb / P[i].NumNgbDM ) / NUMDIMS;
	                            if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
	                            P[i].HsmlDM *= exp(jumpvar);
	                        } else {
	                            P[i].HsmlDM *= 2.0;
	                        }
	                        if((P[i].HsmlDM<RightDM[i])&&(P[i].HsmlDM>LeftDM[i]))
	                        {
	                            if(iter > 1)
	                            {
	                                double hfac = exp(maxjump);
	                                if(P[i].HsmlDM > RightDM[i] / hfac) {P[i].HsmlDM = RightDM[i] / hfac;}
	                                if(P[i].HsmlDM < LeftDM[i] * hfac) {P[i].HsmlDM = LeftDM[i] * hfac;}
	                            }
	                        } else {
	                            if(P[i].HsmlDM>RightDM[i]) P[i].HsmlDM=RightDM[i];
	                            if(P[i].HsmlDM<LeftDM[i]) P[i].HsmlDM=LeftDM[i];
	                            P[i].HsmlDM = pow(P[i].HsmlDM * LeftDM[i] * RightDM[i] , 1.0/3.0);
	                        }
	                    }
	                    else
	                    {
	                        if(RightDM[i] == 0 && LeftDM[i] == 0)
	                        {
	                            char buf[1000];
	                            sprintf(buf, "RightDM[i] == 0 && LeftDM[i] == 0 && P[i].HsmlDM=%g\n", P[i].HsmlDM);
	                            terminate(buf);
	                        }
	                        
	                        if(RightDM[i] == 0 && LeftDM[i] > 0)
	                        {
	                            if (P[i].NumNgbDM > 1)
	                                fac_lim = log( desnumngb / P[i].NumNgbDM ) / NUMDIMS; // this would give desnumgb if constant density (+0.231=2x desnumngb)
	                            else
	                                fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
	                            
	                            if((P[i].NumNgbDM < 2*desnumngb)&&(P[i].NumNgbDM > 0.1*desnumngb))
	                            {
	                                double slope = P[i].DhsmlNgbFactorDM;
	                                if(iter>2 && slope<1) slope = 0.5*(slope+1);
	                                fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
	                                if(iter>=4)
	                                    if(P[i].DhsmlNgbFactorDM==1) fac *= 10; // tries to help with being trapped in small steps
	                                
	                                if(fac < fac_lim+0.231)
	                                {
	                                    P[i].HsmlDM *= exp(fac); // more expensive function, but faster convergence
	                                }
	                                else
	                                {
	                                    P[i].HsmlDM *= exp(fac_lim+0.231);
	                                    // fac~0.26 leads to expected doubling of number if density is constant,
	                                    //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
	                                }
	                            }
	                            else
	                                P[i].HsmlDM *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
	                        }
	                        
	                        if(RightDM[i] > 0 && LeftDM[i] == 0)
	                        {
	                            if (P[i].NumNgbDM > 1)
	                                fac_lim = log( desnumngb / P[i].NumNgbDM ) / NUMDIMS; // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
	                            else
	                                fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
	                            
	                            if (fac_lim < -1.535) fac_lim = -1.535; // decreasing N_ngb by factor ~100
	                            
	                            if((P[i].NumNgbDM < 2*desnumngb)&&(P[i].NumNgbDM > 0.1*desnumngb))
	                            {
	                                double slope = P[i].DhsmlNgbFactorDM;
	                                if(iter>2 && slope<1) slope = 0.5*(slope+1);
	                                fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
	                                if(iter>=4)
	                                    if(P[i].DhsmlNgbFactorDM==1) fac *= 10; // tries to help with being trapped in small steps
	                                
	                                if(fac > fac_lim-0.231)
	                                {
	                                    P[i].HsmlDM *= exp(fac); // more expensive function, but faster convergence
	                                }
	                                else
	                                    P[i].HsmlDM *= exp(fac_lim-0.231); // limiter to prevent --too-- far a jump in a single iteration
	                            }
	                            else
	                                P[i].HsmlDM *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
	                        }
	                    } // RightDM && LeftDM > 0 
	                } // neither maxHsml or minHsml
	                /* resets for max/min values */
	                if(P[i].HsmlDM < minsoft) P[i].HsmlDM = minsoft;
	                if(particle_set_to_minhsml_flag_DM==1) P[i].HsmlDM = minsoft;
	                if(P[i].HsmlDM > maxsoft) P[i].HsmlDM = maxsoft;
	                if(particle_set_to_maxhsml_flag_DM==1) P[i].HsmlDM = maxsoft;
	            } // redo particle
	            else
				{
	                P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
					redo_particle = 0;
				}
	        } // active particle
        } // end DM loop
        
        tend = my_second();
        timecomp1 += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0)
            {
                printf("DMAF: ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
            }
            if(iter > MAXITER)
            {
                printf("DMAF: failed to converge in neighbour iteration in dm_density()\n");
                fflush(stdout);
                endrun(1156);
            }
        }
    }
    while(ntot > 0);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(RightDM);
    myfree(LeftDM);
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
        if(dm_density_isactive(i))
        {
			if(P[i].NumNgbDM > 0) {P[i].NumNgbDM=pow(P[i].NumNgbDM,1./NUMDIMS);} else {P[i].NumNgbDM=0;}
        }
    }
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_DMDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_DMDENSWAIT] += timewait;
    CPU_Step[CPU_DMDENSCOMM] += timecomm;
    CPU_Step[CPU_DMDENSMISC] += timeall - (timecomp + timewait + timecomm);
}




/*! This function represents the core of the density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int dm_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, n;
    int startnode, numngb_inbox, listindex = 0;
	double r2, h2, mass_j, u;
	struct kernel_density kernel;
    struct dm_densdata_in local;
    struct dm_densdata_out out;
    memset(&out, 0, sizeof(struct dm_densdata_out));
    
    if(mode == 0)
        dm_particle2in_density(&local, target);
    else
        local = DM_DensDataGet[target];

    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = DM_DensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
	h2 = local.HsmlDM * local.HsmlDM;

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.HsmlDM, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, 2); // search for DM particles only
            if(numngb_inbox < 0) return -1;

			kernel_hinv(local.HsmlDM, &kernel.hinv, &kernel.hinv3, &kernel.hinv4); 

            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Mass <= 0) continue;
                kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
                kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
                kernel.dp[2] = local.Pos[2] - P[j].Pos[2];
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1);
#endif
                r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
			
			if(r2 < h2)
                {
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                    mass_j = P[j].Mass;
                    kernel.mj_wk = FLT(mass_j * kernel.wk);
                    
                    out.NgbDM += kernel.wk;
                    out.RhoDM += kernel.mj_wk;
                    out.DhsmlNgbDM += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);

					/* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {                        
                        kernel.dv[0] = local.Vel[0] - P[j].Vel[0]; // Use Vel here to estimate divergence since VelPred is not calculated for DM particles
                        kernel.dv[1] = local.Vel[1] - P[j].Vel[1];
                        kernel.dv[2] = local.Vel[2] - P[j].Vel[2];
#ifdef BOX_SHEARING
                        if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
                        if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
#endif
                        out.Particle_DivVelDM -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how Hsml will evolve (particle drift) */
                                                
                    } // kernel.r > 0 //
                }

            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = DM_DensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
        dm_out2particle_density(&out, target, 0);
    else
        DM_DensDataResult[target] = out;
    
    return 0;
}



void *dm_density_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(dm_density_isactive(i))
#define EVALUATION_CALL dm_density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *dm_density_evaluate_secondary(void *p)
{
#define EVALUATION_CALL dm_density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

#endif


/* routine to determine if we need to use dm_density to calculate HsmlDM and rhoDM */
int dm_density_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0;
	if(P[i].Mass <= 0) return 0;
#if defined(DMANNIHILATION)
    if(P[i].Type != 0) return 0;// only gas particles //
#elif defined(DMANNIHILATION_DM)
	if(P[i].Type != 1) return 0;// only DM particles //
#endif
    return 1;
}


