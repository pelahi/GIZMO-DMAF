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

/*! \file dm_rho_gas
 *  \brief density calculation for dark matter particles around dark matter particles
 *
 *  This file contains a loop modeled on the standard gas density computation which
 *    determines the gas density around a given set of DARK MATTER particles and adjusts the smoothing length for this
*     calculation.
 *    The dark matter density is used to set the change of energy injection to dark matter annihilation. 
 *
 * This file was written by Florian List, for GIZMO, based on Phil Hopkins's density routine.
 * Method 2 uses the functionality in this file in order to determine the gas density around DM particles.
 */

#if defined(DMANNIHILATION_DM)

struct kernel_density
{
  double dp[3],dv[3],r;
  double wk, dwk;
  double hinv, hinv3, hinv4;
  double mj_wk, mj_dwk_r;
};

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static struct gas_densdata_in
{
    MyDouble Pos[3];
	MyFloat Vel[3];
    MyFloat Hsml;
    int NodeList[NODELISTLENGTH];
}
*GAS_DensDataIn, *GAS_DensDataGet;

static struct gas_densdata_out
{
    MyLongDouble Ngb, Rho, DhsmlNgb, Particle_DivVel;
}
*GAS_DensDataResult, *GAS_DensDataOut;

void gas_particle2in_density(struct gas_densdata_in *in, int i);
void gas_out2particle_density(struct gas_densdata_out *out, int i, int mode);

void gas_particle2in_density(struct gas_densdata_in *in, int i)
{
    int k;
    for(k = 0; k < 3; k++)
    {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    in->Hsml = P[i].Hsml;
}

void gas_out2particle_density(struct gas_densdata_out *out, int i, int mode)
{
	ASSIGN_ADD(P[i].NumNgb, out->Ngb, mode);
	ASSIGN_ADD(P[i].DensAroundDM, out->Rho, mode);
	ASSIGN_ADD(P[i].DhsmlNgbFactor, out->DhsmlNgb, mode);
	ASSIGN_ADD(P[i].Particle_DivVel, out->Particle_DivVel, mode); 
}

void gas_density_for_DM(void)
{
    MyFloat *Left, *Right;
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
    
    CPU_Step[CPU_GASDENSCOMPUTE] += measure_time();
    
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(gas_density_for_DM_isactive(i))
        {
            P[i].NumNgb = 0;
            Left[i] = Right[i] = 0;
        }
    }
    
    /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                                           sizeof(struct gas_densdata_in) + sizeof(struct gas_densdata_out) +
                                                           sizemax(sizeof(struct gas_densdata_in),sizeof(struct gas_densdata_out))));
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
                pthread_create(&mythreads[j], &attr, gas_density_evaluate_primary, &threadid[j]);
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
                gas_density_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
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
                    printf("Gas density for DM: Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
                           P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);
                    
                    endrun(290495);
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
            
            GAS_DensDataGet = (struct gas_densdata_in *) mymalloc("GAS_DensDataGet", Nimport * sizeof(struct gas_densdata_in));
            GAS_DensDataIn = (struct gas_densdata_in *) mymalloc("GAS_DensDataIn", Nexport * sizeof(struct gas_densdata_in));
            
            /* prepare particle data for export */
            for(j = 0; j < Nexport; j++)
            {
                place = DataIndexTable[j].Index;
                
                gas_particle2in_density(&GAS_DensDataIn[j], place);
                
                memcpy(GAS_DensDataIn[j].NodeList,
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
                        MPI_Sendrecv(&GAS_DensDataIn[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct gas_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DMDENS_A,
                                     &GAS_DensDataGet[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct gas_densdata_in), MPI_BYTE,
                                     recvTask, TAG_DMDENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            tend = my_second();
            timecommsumm1 += timediff(tstart, tend);
            
            myfree(GAS_DensDataIn);
            GAS_DensDataResult = (struct gas_densdata_out *) mymalloc("GAS_DensDataResult", Nimport * sizeof(struct gas_densdata_out));
            GAS_DensDataOut = (struct gas_densdata_out *) mymalloc("GAS_DensDataOut", Nexport * sizeof(struct gas_densdata_out));
            
            /* now do the particles that were sent to us */
            
            tstart = my_second();
            
            NextJ = 0;
            
#ifdef PTHREADS_NUM_THREADS
            for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
                pthread_create(&mythreads[j], &attr, gas_density_evaluate_secondary, &threadid[j]);
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
                gas_density_evaluate_secondary(&mainthreadid);
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
                        MPI_Sendrecv(&GAS_DensDataResult[Recv_offset[recvTask]],
                                     Recv_count[recvTask] * sizeof(struct gas_densdata_out),
                                     MPI_BYTE, recvTask, TAG_DMDENS_B,
                                     &GAS_DensDataOut[Send_offset[recvTask]],
                                     Send_count[recvTask] * sizeof(struct gas_densdata_out),
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
                gas_out2particle_density(&GAS_DensDataOut[j], place, 1);
            }
            tend = my_second();
            timecomp1 += timediff(tstart, tend);
            
            
            myfree(GAS_DensDataOut);
            myfree(GAS_DensDataResult);
            myfree(GAS_DensDataGet);
        }
        while(ndone < NTask);
        
		/* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        tstart = my_second();
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
	        if(gas_density_for_DM_isactive(i)) // This makes sure that for method 1 (2), only gas particles (DM particles) are treated
	        {
	            if(P[i].NumNgb > 0)
	            {
	                P[i].DhsmlNgbFactor *= P[i].Hsml / (NUMDIMS * P[i].NumNgb);
	                P[i].Particle_DivVel /= P[i].NumNgb;
	                /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
	                P[i].NumNgb *= NORM_COEFF * pow(P[i].Hsml,NUMDIMS);
	            } else {
	                P[i].NumNgb = P[i].DhsmlNgbFactor = P[i].Particle_DivVel = 0;
	            }
	            
	            // inverse of SPH volume element (to satisfy constraint implicit in Lagrange multipliers)
	            if(P[i].DhsmlNgbFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
	                P[i].DhsmlNgbFactor = 1 / (1 + P[i].DhsmlNgbFactor);
	            else
	                P[i].DhsmlNgbFactor = 1;
	            P[i].Particle_DivVel *= P[i].DhsmlNgbFactor;		        		       		            
	            
				double minsoft = All.MinHsml;
            	double maxsoft = All.MaxHsml;
			    desnumngb = All.DesNumNgb; // Leave it constant
				desnumngbdev = All.MaxNumNgbDeviation;

				/* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}

				redo_particle = 0; // set redo_particle = 0, check if it needs to be set to 1 in the following

	            /* check if we are in the 'normal' range between the max/min allowed values */
	            if((P[i].NumNgb < (desnumngb - desnumngbdev) && P[i].Hsml < 0.99*maxsoft) ||
	               (P[i].NumNgb > (desnumngb + desnumngbdev) && P[i].Hsml > 1.01*minsoft))
	                redo_particle = 1;
	            
	            /* check maximum kernel size allowed */
	            particle_set_to_maxhsml_flag_DM = 0;
	            if((P[i].Hsml >= 0.99*maxsoft) && (P[i].NumNgb < (desnumngb - desnumngbdev)))
	            {
	                redo_particle = 0;
	                if(P[i].Hsml == maxsoft)
	                {
	                    /* iteration at the maximum value is already complete */
	                    particle_set_to_maxhsml_flag_DM = 0;
	                } else {
	                    /* ok, the particle needs to be set to the maximum, and iterated one more time */
	                    redo_particle = 1;
	                    P[i].Hsml = maxsoft;
	                    particle_set_to_maxhsml_flag_DM = 1;
	                }
	            }
	            
	            /* check minimum kernel size allowed */
	            particle_set_to_minhsml_flag_DM = 0;
	            if((P[i].Hsml <= 1.01*minsoft) && (P[i].NumNgb > (desnumngb + desnumngbdev)))
	            {
	                redo_particle = 0;
	                if(P[i].Hsml == minsoft)
	                {
	                    /* this means we've already done an iteration with the MinHsml value, so the
	                     neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
	                    particle_set_to_minhsml_flag_DM = 0;
	                } else {
	                    /* ok, the particle needs to be set to the minimum, and iterated one more time */
	                    redo_particle = 1;
	                    P[i].Hsml = minsoft;
	                    particle_set_to_minhsml_flag_DM = 1;
	                }
	            }		               		            
	            
	            if(redo_particle)
	            {
	                if(iter >= MAXITER - 10)
	                {
	                    printf("Gas density for DM particle loop parameters:\n i=%d task=%d ID=%llu Type=%d Hsml=%g dhsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)\n",
	                           i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, P[i].Hsml, P[i].DhsmlNgbFactor, Left[i], Right[i],
	                           (float) P[i].NumNgb, Right[i] - Left[i], particle_set_to_maxhsml_flag_DM, particle_set_to_minhsml_flag_DM, minsoft,
	                           maxsoft, desnumngb, desnumngbdev, redo_particle, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);							
						printf("SLOW CONVERGENCE IN GAS DENSITY FOR DM PARTICLE LOOP!\n");
	                }
	                
	                /* need to redo this particle */
	                npleft++;
	                
	                if(Left[i] > 0 && Right[i] > 0)
	                    if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
	                    {
	                        /* this one should be ok */
	                        npleft--;
	                        P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
							redo_particle = 0;
	                        continue;
	                    }
	                
	                if((particle_set_to_maxhsml_flag_DM==0)&&(particle_set_to_minhsml_flag_DM==0))
	                {
	                    if(P[i].NumNgb < (desnumngb - desnumngbdev))
	                        Left[i] = DMAX(P[i].Hsml, Left[i]);
	                    else
	                    {
	                        if(Right[i] != 0)
	                        {
	                            if(P[i].Hsml < Right[i])
	                                Right[i] = P[i].Hsml;
	                        }
	                        else
	                            Right[i] = P[i].Hsml;
	                    }
	                    
	                    // right/left define upper/lower bounds from previous iterations
	                    if(Right[i] > 0 && Left[i] > 0)
	                    {
	                        // geometric interpolation between right/left //
	                        double maxjump=0;
	                        if(iter>1) {maxjump = 0.2*log(Right[i]/Left[i]);}
	                        if(P[i].NumNgb > 1)
	                        {
	                            double jumpvar = P[i].DhsmlNgbFactor * log( desnumngb / P[i].NumNgb ) / NUMDIMS;
	                            if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
	                            P[i].Hsml *= exp(jumpvar);
	                        } else {
	                            P[i].Hsml *= 2.0;
	                        }
	                        if((P[i].Hsml<Right[i])&&(P[i].Hsml>Left[i]))
	                        {
	                            if(iter > 1)
	                            {
	                                double hfac = exp(maxjump);
	                                if(P[i].Hsml > Right[i] / hfac) {P[i].Hsml = Right[i] / hfac;}
	                                if(P[i].Hsml < Left[i] * hfac) {P[i].Hsml = Left[i] * hfac;}
	                            }
	                        } else {
	                            if(P[i].Hsml>Right[i]) P[i].Hsml=Right[i];
	                            if(P[i].Hsml<Left[i]) P[i].Hsml=Left[i];
	                            P[i].Hsml = pow(P[i].Hsml * Left[i] * Right[i] , 1.0/3.0);
	                        }
	                    }
	                    else
	                    {
	                        if(Right[i] == 0 && Left[i] == 0)
	                        {
	                            char buf[1000];
	                            sprintf(buf, "Right[i] == 0 && Left[i] == 0 && P[i].Hsml=%g\n", P[i].Hsml);
	                            terminate(buf);
	                        }
	                        
	                        if(Right[i] == 0 && Left[i] > 0)
	                        {
	                            if (P[i].NumNgb > 1)
	                                fac_lim = log( desnumngb / P[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (+0.231=2x desnumngb)
	                            else
	                                fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
	                            
	                            if((P[i].NumNgb < 2*desnumngb)&&(P[i].NumNgb > 0.1*desnumngb))
	                            {
	                                double slope = P[i].DhsmlNgbFactor;
	                                if(iter>2 && slope<1) slope = 0.5*(slope+1);
	                                fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
	                                if(iter>=4)
	                                    if(P[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
	                                
	                                if(fac < fac_lim+0.231)
	                                {
	                                    P[i].Hsml *= exp(fac); // more expensive function, but faster convergence
	                                }
	                                else
	                                {
	                                    P[i].Hsml *= exp(fac_lim+0.231);
	                                    // fac~0.26 leads to expected doubling of number if density is constant,
	                                    //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
	                                }
	                            }
	                            else
	                                P[i].Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
	                        }
	                        
	                        if(Right[i] > 0 && Left[i] == 0)
	                        {
	                            if (P[i].NumNgb > 1)
	                                fac_lim = log( desnumngb / P[i].NumNgb ) / NUMDIMS; // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
	                            else
	                                fac_lim = 1.4; // factor ~66 increase in N_NGB in constant-density medium
	                            
	                            if (fac_lim < -1.535) fac_lim = -1.535; // decreasing N_ngb by factor ~100
	                            
	                            if((P[i].NumNgb < 2*desnumngb)&&(P[i].NumNgb > 0.1*desnumngb))
	                            {
	                                double slope = P[i].DhsmlNgbFactor;
	                                if(iter>2 && slope<1) slope = 0.5*(slope+1);
	                                fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
	                                if(iter>=4)
	                                    if(P[i].DhsmlNgbFactor==1) fac *= 10; // tries to help with being trapped in small steps
	                                
	                                if(fac > fac_lim-0.231)
	                                {
	                                    P[i].Hsml *= exp(fac); // more expensive function, but faster convergence
	                                }
	                                else
	                                    P[i].Hsml *= exp(fac_lim-0.231); // limiter to prevent --too-- far a jump in a single iteration
	                            }
	                            else
	                                P[i].Hsml *= exp(fac_lim); // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
	                        }
	                    } // Right && Left > 0 
	                } // neither maxHsml or minHsml
	                /* resets for max/min values */
	                if(P[i].Hsml < minsoft) P[i].Hsml = minsoft;
	                if(particle_set_to_minhsml_flag_DM==1) P[i].Hsml = minsoft;
	                if(P[i].Hsml > maxsoft) P[i].Hsml = maxsoft;
	                if(particle_set_to_maxhsml_flag_DM==1) P[i].Hsml = maxsoft;
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
                printf("Gas density for DM particles: ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
                       (int) (ntot / 1000000000), (int) (ntot % 1000000000));
            }
            if(iter > MAXITER)
            {
                printf("Gas density for DM particles: failed to converge in neighbour iteration in gas_density_for_DM()\n");
                fflush(stdout);
                endrun(1156);
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
        if(gas_density_for_DM_isactive(i))
        {
			if(P[i].NumNgb > 0) {P[i].NumNgb=pow(P[i].NumNgb,1./NUMDIMS);} else {P[i].NumNgb=0;}
        }
    }
    
    /* collect some timing information */
    t1 = WallclockTime = my_second();
    timeall += timediff(t0, t1);
    timecomp = timecomp1 + timecomp2;
    timewait = timewait1 + timewait2;
    timecomm = timecommsumm1 + timecommsumm2;
    CPU_Step[CPU_GASDENSCOMPUTE] += timecomp;
    CPU_Step[CPU_GASDENSWAIT] += timewait;
    CPU_Step[CPU_GASDENSCOMM] += timecomm;
    CPU_Step[CPU_GASDENSMISC] += timeall - (timecomp + timewait + timecomm);
}




/*! This function represents the core of the density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int gas_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int j, n, k;
    int startnode, numngb_inbox, listindex = 0;
	double r2, h2, mass_j, u;
	struct kernel_density kernel;
    struct gas_densdata_in local;
    struct gas_densdata_out out;
    memset(&out, 0, sizeof(struct gas_densdata_out));
    
    if(mode == 0)
        gas_particle2in_density(&local, target);
    else
        local = GAS_DensDataGet[target];

    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = GAS_DensDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    
	h2 = local.Hsml * local.Hsml;

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
			// Search for neighbours in a symmetric way, so, also include gas particles that can see the DM particle and not the other way round!
			// Motivation for this: https://authors.library.caltech.edu/87408/1/sty674.pdf
			// This must be the same as in dm_fb.c since DensAroundDM is used as a norming factor for the injection and the neighbours must therefore be the same. 
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist); // search for gas only
            if(numngb_inbox < 0) return -1;

			kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4); 

            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}                
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) continue; // outside kernel (in both 'directions') //
				
				// NOTE: do NOT require r2 > 0 here because it can be that DM and gas reside at same position!
	
				// Calculate kernel properties of DM particle
                kernel.r = sqrt(r2);
                u = kernel.r * kernel.hinv;
                kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
				if((kernel.wk <= 0)||(isnan(kernel.wk))) continue;

				// Calculate kernel properties of gas particle
                mass_j = P[j].Mass;
                kernel.mj_wk = FLT(mass_j * kernel.wk);				
                
				// Calculate hydro quantities
                out.Ngb += kernel.wk;
                out.Rho += kernel.mj_wk;
                out.DhsmlNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);

				/* no need to take care of avoiding particle self contribution since DM searches for gas particles */                      
                kernel.dv[0] = local.Vel[0] - P[j].Vel[0]; // Use Vel here to estimate divergence since VelPred is not calculated for DM particles
                kernel.dv[1] = local.Vel[1] - P[j].Vel[1];
                kernel.dv[2] = local.Vel[2] - P[j].Vel[2];
#ifdef BOX_SHEARING
                if(local.Pos[0] - P[j].Pos[0] > +boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
                if(local.Pos[0] - P[j].Pos[0] < -boxHalf_X) {kernel.dv[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
#endif
                out.Particle_DivVel -= kernel.dwk * (kernel.dp[0] * kernel.dv[0] + kernel.dp[1] * kernel.dv[1] + kernel.dp[2] * kernel.dv[2]) / kernel.r;
                /* this is the -particle- divv estimator, which determines how Hsml will evolve (particle drift) */				                                                  				
            }
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = GAS_DensDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                    startnode = Nodes[startnode].u.d.nextnode;	/* open it */
            }
        }
    }
    
    if(mode == 0)
        gas_out2particle_density(&out, target, 0);
    else
        GAS_DensDataResult[target] = out;
    
    return 0;
}



void *gas_density_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(gas_density_for_DM_isactive(i))
#define EVALUATION_CALL gas_density_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *gas_density_evaluate_secondary(void *p)
{
#define EVALUATION_CALL gas_density_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

#endif


/* routine to determine if we need to use gas_density to calculate Hsml and rhoDM */
int gas_density_for_DM_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0;
	if(P[i].Mass <= 0) return 0;
	if(P[i].Type != 1) return 0;// only DM particles //
    return 1;
}


