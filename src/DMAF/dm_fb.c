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

/*! \file dm_fb
 *  \brief DM annihilation feedback calculation at dark matter particles
 * This file was written by Florian List, for GIZMO, based on Phil Hopkins's feedback routines.
 * Method 2 (activate using DMANNIHILATION_DM)
  - gas density is calculated at each DM particle (in dm_rho_gas.c)
  - DM annihilation rate at DM particle is calculated (in dm_annihilation.c)
  - energy is injected into surrounding gas particles (in dm_fb.c)
 */

/* 1. Loop: AreaSumDM is calculated if solid angle weighted injection, 
/           maximum gas time bins for the next time step are computed,  
			subsequent DM time steps may be set if DMAF_LIMIT_DM_TIMESTEP.
/           (loop_fb == 0)
/  2. Loop: DMAF energy is stored (loop_fb == 1)
/  If energy distribution is calculated based on density:
/     Normalising factor is DensAroundDM, has already been computed in dm_rho_gas.c. */

#ifdef DMANNIHILATION_DM

struct kernel_addDMAF {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};

/* define structures to use below */
struct addDMAFdata_in
{
    MyDouble Pos[3], Hsml, E_DMAF, wt_sum;
    int NodeList[NODELISTLENGTH];
	int timebin_DM;
	int DMAF_MaxTimebin;
#ifndef DMAF_DENSITY_WEIGHTED
	MyDouble V_i;
#endif
}
*addDMAFDataIn, *addDMAFDataGet;

           
/* Get gas density around DM */
void particle2in_addDMAF(struct addDMAFdata_in *in, int i, short int loop_fb);
void particle2in_addDMAF(struct addDMAFdata_in *in, int i, short int loop_fb)
{
	MyDouble heat_source, k_rho_dm_times_mass; 
	MyDouble new_cf_a3inv, new_cf_atime;
    int k;
	in->Hsml=PPP[i].Hsml; 
#ifndef DMAF_DENSITY_WEIGHTED
	if(loop_fb)
		in->wt_sum = P[i].AreaSumDM;	
	double heff = PPP[i].Hsml / PPP[i].NumNgb; 
	in->V_i = heff*heff*heff;
#else
	if(loop_fb)	
		in->wt_sum = P[i].DensAroundDM;	
#endif
	in->timebin_DM = P[i].TimeBin;
#ifdef DMAF_LIMIT_DM_TIMESTEP
	in->DMAF_MaxTimebin = P[i].DMAF_MaxTimebin;
#endif
	for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];}
	
	// Calculate scale factor at the end of the coming time step
	if(All.ComovingIntegrationOn)
	{				
		new_cf_atime = All.TimeBegin * exp((All.Ti_Current + (P[i].TimeBin ? (1 << P[i].TimeBin) : 0)) * All.Timebase_interval);
		new_cf_a3inv = 1.0 / (new_cf_atime * new_cf_atime * new_cf_atime);
	}
	else
	{
		new_cf_atime = new_cf_a3inv = 1.0;
	}

	// NOTE: do NOT take into account gas mass / density here, this will be done at each gas particle!
	k_rho_dm_times_mass = P[i].DensityDM * P[i].Mass * new_cf_a3inv * All.HubbleParam;
	heat_source = All.UnitTime_in_s * k_rho_dm_times_mass * C_km_s * C_km_s * CONVERT_FAC * All.CrossSection / All.ParticleMass; // this is in units of energy / internal time
	// Set total energy rate created by DM annihilation
	in->E_DMAF = 1.0 / new_cf_atime * heat_source;

// For debugging and testing: INJECT_ENERGY_DM
#ifdef INJECT_ENERGY_DM
	if( ((P[i].ID == All.EnergyID) || (All.EnergyID == -1)) && (P[i].Type == 1) )
	{
		// constant:
		in->E_DMAF = All.EnergySource;		
		// linear in time (uncomment):
		// in->E_DMAF = All.EnergySource * (All.Time + 0.390625); // needs to be the END of the current DM time step
	}
	else
		in->E_DMAF = 0.0;
#endif

}

struct addDMAFdata_out
{
#ifdef DMAF_LIMIT_DM_TIMESTEP
    int NewBin;
#endif
#ifndef DMAF_DENSITY_WEIGHTED
	double AreaSum;
#endif	
}
*addDMAFDataResult, *addDMAFDataOut;

void out2particle_addDMAF(struct addDMAFdata_out *out, int i, int mode, short int loop_fb);
void out2particle_addDMAF(struct addDMAFdata_out *out, int i, int mode, short int loop_fb)
{
#ifndef DMAF_DENSITY_WEIGHTED
	if(!loop_fb)
		ASSIGN_ADD(P[i].AreaSumDM, out->AreaSum, mode);
#endif
#ifdef DMAF_LIMIT_DM_TIMESTEP
	// if NewBin is set: update DMAF_MaxTimebin such that time step is updated next time
	if(!loop_fb)
		P[i].DMAF_MaxTimebin = (out->NewBin) ? out->NewBin : P[i].DMAF_MaxTimebin; 
#endif
}

int addDMAF_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, short int loop_fb)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n, max_gas_bin;
    double u,r2,h2,wk,wk_loc,u_DMAF;
    struct kernel_addDMAF kernel;
    struct addDMAFdata_in local;
    struct addDMAFdata_out out;
    memset(&out, 0, sizeof(struct addDMAFdata_out));
#ifndef DMAF_DENSITY_WEIGHTED
	if(!loop_fb)
		out.AreaSum = 0;
#endif
    
    /* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_addDMAF(&local, target, loop_fb);} else {local = addDMAFDataGet[target];}
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml;    
    
#ifdef DMAF_LIMIT_DM_TIMESTEP
	out.NewBin = TIMEBINS;
#endif		

    /* Now start the actual DMAF computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = addDMAFDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {									
			// Search for neighbours in a symmetric way, so, also include gas particles that can see the DM particle and not the other way round!
			// Motivation for this: https://authors.library.caltech.edu/87408/1/sty674.pdf 
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist); // search for gas only
						
            if(numngb_inbox < 0) return -1;			

			kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

			// This is the loop where the DMAF energy is stored / AreaSum is calculated
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

				// in loop 0: save maximum time bin for gas, given by minimum time bin of injecting DM					
				if(!loop_fb)
				{
					// save maximum time bin for gas
					max_gas_bin = local.timebin_DM;
					if(max_gas_bin < P[j].DMAF_MaxTimebin && max_gas_bin)
						P[j].DMAF_MaxTimebin = max_gas_bin;
				}

#ifdef DMAF_LIMIT_DM_TIMESTEP
				// in loop 0: if gas particle has much smaller timestep than injecting DM particle: reduce next time step of DM.				
				if(P[j].TimeBin)
				{
					new_DM_bin = P[j].TimeBin + DMAF_LIMIT_DM_TIMESTEP;
					out.NewBin = DMIN(out.NewBin, new_DM_bin);
				}
#endif	

#ifndef DMAF_DENSITY_WEIGHTED
				// Calculate kernel properties of gas particle
				double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j;
				double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j;
				double V_j = P[j].Mass / SphP[j].Density;
				kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);
  				if(local.V_i<0 || isnan(local.V_i)) {local.V_i = 0;}
                if(V_j<0 || isnan(V_j)) {V_j = 0;}
				// Calculate effective face area //
				double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j);
				// Calculate corresponding geometric weight //
				if(kernel.r > 0)
				{
					wk_loc = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r)));
				}
				else if(kernel.r == 0)
				{
					wk_loc = 1.0; // total solid angle / 4 * PI
				} 
				else
				{
					continue;
				}
				// in loop 0: add contribution to the total area sum and exit loop
				if(!loop_fb)
				{									
					out.AreaSum += wk_loc;									
				}	
				// in loop 1: calculate normalized weight
				else
				{					
					wk = wk_loc / local.wt_sum; // (solid angle weighted)  
					if((wk <= 0)||(isnan(wk))) continue; // no point in going further, there's no physical weight here
				}
							 			    
#else				
				wk = P[j].Mass * kernel.wk / local.wt_sum; // normalized weight function (simply density weighted)    
#endif
				// in loop 0: continue
				if(!loop_fb)
				{
					continue; 
				}	

                /* calculate energy */
				u_DMAF = wk * local.E_DMAF; // Note that local.E_DMAF is an energy, NOT specific energy. 
														// but: do NOT divide by particle mass here, it could change over time! Conserve TOTAL energy!											
								
				// Store energy in correct bin
				if(local.timebin_DM >= P[j].TimeBin) // regular case: injecting into gas that has the same or a smaller time step currently
				{	
					SphP[j].DMAF_Dtu[local.timebin_DM] += u_DMAF;					
					SphP[j].DMAF_Dtu_tot += u_DMAF;
				}
				else // if DM encounters a gas particle with a larger time step than itself				 
				{											
					if(TimeBinActive[P[j].TimeBin]) // if gas particle has just become active again: gas should have reduced its time step -> something is wrong!
					{
						printf("Gas is active and has a larger time set than an injecting DM particle! This should not happen.\n");
						endrun(134730);
					}
					else // gas particle is currently inactive
					{
						double ratio = 1.0 * (local.timebin_DM ? (1 << local.timebin_DM) : 1) / (P[j].TimeBin ? (1 << P[j].TimeBin) : 1); 
						SphP[j].DMAF_Dtu[local.timebin_DM] += ratio * u_DMAF;
						SphP[j].DMAF_Dtu_tot += ratio * u_DMAF;
					}																									
				} // local.timebin_DM >= P[j].TimeBin																		    			
            } // for(n = 0; n < numngb; n++)									
        } // while(startnode >= 0)				
	
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
			{
                startnode = addDMAFDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}    /* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addDMAF(&out, target, 0, loop_fb);} else {addDMAFDataResult[target] = out;}
	
    return 0;
} // int addDMAF_evaluate



void DMAF_calc(short int loop_fb)
{
    int j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    long long n_exported = 0;
    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct addDMAFdata_in) + sizeof(struct addDMAFdata_out) + sizemax(sizeof(struct addDMAFdata_in),sizeof(struct addDMAFdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    NextParticle = FirstActiveParticle;	/* begin with this index */
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
            pthread_create(&mythreads[j], &attr, addDMAF_evaluate_primary, &threadid[j]);
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
            addDMAF_evaluate_primary(&mainthreadid, loop_fb);	/* do local particles and prepare export list */
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
#endif
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
            if(NextParticle == save_NextParticle) {endrun(116609);} /* in this case, the buffer is too small to process even a single particle */
            int new_export = 0;
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1) k = j + 1;
                    
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
        for(j = 0; j < NTask; j++) Send_count[j] = 0;
        for(j = 0; j < Nexport; j++) Send_count[DataIndexTable[j].Task]++;
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
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
        addDMAFDataGet = (struct addDMAFdata_in *) mymalloc("addDMAFDataGet", Nimport * sizeof(struct addDMAFdata_in));
        addDMAFDataIn = (struct addDMAFdata_in *) mymalloc("addDMAFDataIn", Nexport * sizeof(struct addDMAFdata_in));
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_addDMAF(&addDMAFDataIn[j], place, loop_fb);
            memcpy(addDMAFDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        /* exchange particle data */
        int TAG_TO_USE = TAG_DMAFLOOP_1A;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&addDMAFDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct addDMAFdata_in), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &addDMAFDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct addDMAFdata_in), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        myfree(addDMAFDataIn);
        addDMAFDataResult = (struct addDMAFdata_out *) mymalloc("addDMAFDataResult", Nimport * sizeof(struct addDMAFdata_out));
        addDMAFDataOut = (struct addDMAFdata_out *) mymalloc("addDMAFDataOut", Nexport * sizeof(struct addDMAFdata_out));
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, addDMAF_evaluate_secondary, &threadid[j]);
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
            addDMAF_evaluate_secondary(&mainthreadid, loop_fb);
        }
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) pthread_join(mythreads[j], NULL);
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        TAG_TO_USE = TAG_DMAFLOOP_1B;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&addDMAFDataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct addDMAFdata_out), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &addDMAFDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct addDMAFdata_out), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_addDMAF(&addDMAFDataOut[j], place, 1, loop_fb);
        }
        myfree(addDMAFDataOut);
        myfree(addDMAFDataResult);
        myfree(addDMAFDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);

// now: if desired, reduce the DM particle mass
// Note that at the first call, all particles are in time bin 0 -> dt = 0, no mass loss takes place.
#ifndef DMAF_KEEP_DM_MASS
	MyDouble new_cf_atime, new_cf_a3inv;
	if(loop_fb) // Compute mass loss in second loop
	{
		double dt;
		int i;
		for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
		{
			if(P[i].Type != 1) continue; // only need to consider DM particles
#ifndef WAKEUP
    	    dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; // dloga to dt_physical
#else
    	    dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a; // get particle timestep //
#endif
			// Calculate scale factor at the end of the coming time step
			if(All.ComovingIntegrationOn)
			{				
				new_cf_atime = All.TimeBegin * exp((All.Ti_Current + (P[i].TimeBin ? (1 << P[i].TimeBin) : 0)) * All.Timebase_interval);
				new_cf_a3inv = 1.0 / (new_cf_atime * new_cf_atime * new_cf_atime);
			}
			else
			{
				new_cf_atime = new_cf_a3inv = 1.0;
			}
			// Store mass loss
#ifndef INJECT_ENERGY_DM
			P[i].DMAF_dM = - dt * CONVERT_FAC * All.CrossSection * All.UnitTime_in_s * All.HubbleParam * P[i].DensityDM * P[i].Mass * new_cf_a3inv / All.ParticleMass / new_cf_atime; 
#else
			// Test case: inject given amount of energy dE -> reduce mass by dM = dE / c^2
			if( ((P[i].ID == All.EnergyID) || (All.EnergyID == -1)) && (P[i].Type == 1) )
			{
				P[i].DMAF_dM = - dt * All.EnergySource / (C_km_s * C_km_s);
			}
			else
			{
				P[i].DMAF_dM = 0;
			}
#endif
		}
	}
#endif
}

// This function injects the DMAF energy calculated at the beginning of the time step
void inject_DMAF(void)
{
	int i, k;
	MyDouble dt;

	// double check_mass_energy = 0.0;

	// inject energy
	for(i = FirstActiveParticle; i != -1; i = NextActiveParticle[i])
	{	// GAS particles
		if(P[i].Type == 0)
		{
#ifndef WAKEUP
			dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; // dloga to dt_physical
#else
			dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a; // get particle timestep
#endif
			do_DMAF_energy_injection(i, dt);
		}

		// DM particles
#ifndef DMAF_KEEP_DM_MASS
		if(P[i].Type == 1)
		{						 			
			// check_mass_energy += P[i].DMAF_dM * (C_km_s * C_km_s);

			P[i].DensityDM *= (1.0 + P[i].DMAF_dM/P[i].Mass); // inject mass at constant particle volume -> scale DM density	
			P[i].Mass += P[i].DMAF_dM;
			if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}	
			P[i].DMAF_dM = 0; // reset mass change
		}
		// printf("Energy lost as mass: %g.\n", check_mass_energy);
#endif
	}	
	
	// zero out all DMAF energy rates that will be updated now, from All.HighestActiveTimeBin downwards to P[i].TimeBin	
	// moreover, set DMAF_Dtu_tot
	for(i = 0; i < N_gas; i++)
	{		
		// this must only be done for gas particles that are active
		if(TimeBinActive[P[i].TimeBin])
		{
			SphP[i].DMAF_Dtu_tot = 0.0; // zero out DMAF_Dtu_tot, will be set now and then updated in update_DMAF
			for(k = 0; k < TIMEBINS; k++)
			{
				if(k <= All.HighestActiveTimeBin)
                {
#ifdef DMAF_INJECT_AFTER_DM_STEP
					if(TimeBinActive[k]) // if injection is only after time step of injecting DM particle: zero out energy rate only if DM particle is active
					{
#endif
						SphP[i].DMAF_Dtu[k] = 0.0;
#ifdef DMAF_INJECT_AFTER_DM_STEP
					}
#endif				
				}			
			SphP[i].DMAF_Dtu_tot += SphP[i].DMAF_Dtu[k]; 							
			}
		}					
    }
}

// This function contains the core of the DMAF energy injection
void do_DMAF_energy_injection(int i, MyDouble dt)
{
	MyDouble energy_inject;
	int k;
#ifdef DMAF_CORRECT_BIAS
	int next_time_DM;
	MyDouble correction_fac, next_a_DM;
#endif

	// inject energy
	for(k = 0; k < TIMEBINS; k++)
	{
#ifdef DMAF_INJECT_AFTER_DM_STEP			
		if(TimeBinActive[k]) // inject energy only after each DM time step
		{
			dt = (1 << k) * All.Timebase_interval / All.cf_hubble_a; // dt of injecting DM particle
#endif
			energy_inject = dt * SphP[i].DMAF_Dtu[k] / P[i].Mass; // need to divide by current particle mass since total energy is stored in DMAF_Dtu				
#ifdef DMAF_CORRECT_BIAS							
			correction_fac = 1.0;
			if(All.ComovingIntegrationOn)
			{							
				if(k >= P[i].TimeBin && All.Ti_Current)
				{
					next_time_DM = (P[i].Ti_current / (1 << k)) * (1 << k);
					while(next_time_DM < P[i].Ti_current)
					{
						next_time_DM += (1 << k);
					} 			 
											
					next_a_DM = All.TimeBegin * exp(next_time_DM * All.Timebase_interval); // scale factor at the end of DM time step (DMAF_Dtu[k] was set using this scale factor)
					correction_fac = pow(next_a_DM / All.Time, 4.0); // 3.0 because of density scaling w.r.t. expansion, 1.0 because d(log a) = 1/a * da						
				}
			}
			energy_inject *= correction_fac;			
#endif
			SphP[i].InternalEnergy     += energy_inject; 
			SphP[i].InternalEnergyPred += energy_inject;				 		
				
#ifdef DMAF_INJECT_AFTER_DM_STEP
		}			
#endif
	}		
}

// This function calculates the DMAF energy
void update_DMAF(void)
{	
	int i;

	// DMANNIHILATION at DM: calculate densities
	dm_density(); // computes dark matter density:
				  // in case of DMANNIHILATION_DM: around DM particles 
	gas_density_for_DM();  // computes gas density around DM particles, determines smoothing length depending on NumNgb 
						   // this is only needed in order to find the right search radius for energy injection

	// compute store DMAF energy		
	DMAF_calc(0); // Calculate area sum of neighbouring gas particles
	find_timesteps(1); // Reset gas time bins
	DMAF_calc(1); // Calculate energy injection from DMAF at DM particles into gas particles
	find_timesteps(1); // Update gas time bins depending on energy injected

	// in first call, delete DMAF_Dtu since time bins are not known yet - it will be called again when time bins are known.
	if(!All.Ti_Current)
	{
		for(i = 0; i < N_gas; i++)
		{		
			SphP[i].DMAF_Dtu[0] = 0.0;
			SphP[i].DMAF_Dtu_tot = 0.0;		
		}
	}	
}

// This function calculates an energy prediction such that the internal energy output corresponds to the time tend.
// This is needed because the energy is injected at the end of each time step, which means that energy is "missing"
// while large time steps have not finished.
MyDouble drift_DMAF(int i, integertime tstart, integertime tend)
{
	double dt_entr = 1.0 * (tend - tstart) * All.Timebase_interval / All.cf_hubble_a;
	int k, next_time_DM;
	MyDouble energy_inject, energy_tmp;
	energy_inject = 0;
#ifdef DMAF_CORRECT_BIAS
	MyDouble correction_fac, next_a_DM, next_a_gas;
#endif

	// predict energy at tend, if particle is not active
	for(k = 0; k < TIMEBINS; k++)
	{
		energy_tmp = dt_entr * SphP[i].DMAF_Dtu[k] / P[i].Mass; // need to divide by current particle mass since total energy is stored in DMAF_Dtu	
#if defined(DMAF_CORRECT_BIAS) && !defined(DMAF_INJECT_AFTER_DM_STEP)
		correction_fac = 1.0;
		if(All.ComovingIntegrationOn)
		{							
			if(k >= P[i].TimeBin && All.Ti_Current)
			{
				next_time_DM = (P[i].Ti_current / (1 << k)) * (1 << k);
				while(next_time_DM < P[i].Ti_current)
				{
					next_time_DM += (1 << k);
				} 
							 										
				next_a_DM = All.TimeBegin * exp(next_time_DM * All.Timebase_interval); // scale factor at the end of DM time step (DMAF_Dtu[k] was set using this scale factor)
				next_a_gas = All.TimeBegin * exp((P[i].Ti_begstep + (1 << P[i].TimeBin)) * All.Timebase_interval); // scale factor at the end of gas time step (energy will be injected using this scale factor)
				correction_fac = pow(next_a_DM / next_a_gas, 4.0); // 3.0 because of density scaling w.r.t. expansion, 1.0 because d(log a) = 1/a * da
			}			
			energy_tmp *= correction_fac;	
		}	
#endif	
		energy_inject += energy_tmp;
	}
	return energy_inject;	
}


int addDMAF_evaluate_active_check(int i);
int addDMAF_evaluate_active_check(int i)
{
    if(P[i].Type != 1) return 0; // only DM particles trigger annihilation feedback
    if(P[i].Mass <= 0) return 0;
    return 1;
}


void *addDMAF_evaluate_primary(void *p, short int loop_fb)
{
#define CONDITION_FOR_EVALUATION if(addDMAF_evaluate_active_check(i)==1)
#define EVALUATION_CALL addDMAF_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist, loop_fb)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *addDMAF_evaluate_secondary(void *p, short int loop_fb)
{
#define EVALUATION_CALL addDMAF_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist, loop_fb);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}



#endif /* DMANNIHILATION_DM */

