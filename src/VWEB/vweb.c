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

/*! \file vweb
 *  \brief Estimation of the velocity shear tensor of gas particles.
 */

#ifdef VWEB

struct kernel_vweb {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};

/* define structures to use below */
struct vwebdata_in
{
    MyDouble Pos[3], Vel[3], Hsml;
	int NodeList[NODELISTLENGTH];
}
*vwebDataIn, *vwebDataGet;

           
/* Get velocity shear tensor */
void particle2in_vweb(struct vwebdata_in *in, int i);
void particle2in_vweb(struct vwebdata_in *in, int i)
{
    int k;
	in->Hsml = PPP[i].Hsml; 
	for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k];}
	for(k=0;k<3;k++) {in->Vel[k]=SphP[i].VelPred[k];}
}

struct vwebdata_out
{
	double VelShear[3][3];
}
*vwebDataResult, *vwebDataOut;

void out2particle_vweb(struct vwebdata_out *out, int i, int mode);
void out2particle_vweb(struct vwebdata_out *out, int i, int mode)
{
	int k, l;
	for(k=0;k<3;k++)
	{
		for(l=0;l<3;l++)
		{
			ASSIGN_ADD(P[i].VelShear[k][l], out->VelShear[k][l] / SphP[i].Density, mode);	// must divide by density here
		}
	}
}

int vweb_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int startnode, numngb_inbox, listindex = 0, j, k, l, n, max_gas_bin;
    double u,r2,h2,wk;
	double v_kl, v_lk;
    struct kernel_vweb kernel;
    struct vwebdata_in local;
    struct vwebdata_out out;
    memset(&out, 0, sizeof(struct vwebdata_out));

    /* Load the data for the particle */
    if(mode == 0) {particle2in_vweb(&local, target);} else {local = vwebDataGet[target];}
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml;    
    
    /* Now start the velocity shear computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = vwebDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }

    while(startnode >= 0)
    {
        while(startnode >= 0)
        {									
			// Search for neighbours
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, 1);  // only search for gas
						
            if(numngb_inbox < 0) return -1;			

			kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

			// This is the loop where the velocity shear is calculated
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}                

                if(r2>h2) continue; // outside kernel//
				if(r2==0) continue; // skip itself //

                // Calculate kernel properties
                kernel.r = sqrt(r2);
                u = kernel.r * kernel.hinv;
                kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 0);
                if((kernel.wk <= 0)||(isnan(kernel.wk))) continue;

				// Velocity shear tensor here!		
				for(k = 0; k < 3; k++)
				{
					for(l = 0; l <= k; l++)
					{
						v_kl = P[j].Mass * (P[j].Vel[k] - local.Vel[k]) * kernel.dwk * kernel.dp[l] / kernel.r;
						v_lk = P[j].Mass * (P[j].Vel[l] - local.Vel[l]) * kernel.dwk * kernel.dp[k] / kernel.r;
						out.VelShear[k][l] += (- 1.0 / (2.0 * 100 * 1.0e-3 * All.HubbleParam)) * (v_kl + v_lk);
						if(l < k)
							out.VelShear[l][k] += (- 1.0 / (2.0 * 100 * 1.0e-3 * All.HubbleParam)) * (v_kl + v_lk);  // 100 * All.HubbleParam: H0 in km/s/Mpc -> x 1e-3 in km/s/kpc
					}
				}		
														    			
            } // for(n = 0; n < numngb; n++)									
        } // while(startnode >= 0)				
	
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
			{
                startnode = vwebDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}    /* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_vweb(&out, target, 0);} else {vwebDataResult[target] = out;}
	
    return 0;
} // int vweb_evaluate


void vweb_calc()
{
    int j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    long long n_exported = 0;
    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct vwebdata_in) + sizeof(struct vwebdata_out) + sizemax(sizeof(struct vwebdata_in),sizeof(struct vwebdata_out))));
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
            pthread_create(&mythreads[j], &attr, vweb_evaluate_primary, &threadid[j]);
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
            vweb_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
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
            if(NextParticle == save_NextParticle) {endrun(116629);} /* in this case, the buffer is too small to process even a single particle */
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
        vwebDataGet = (struct vwebdata_in *) mymalloc("vwebDataGet", Nimport * sizeof(struct vwebdata_in));
        vwebDataIn = (struct vwebdata_in *) mymalloc("vwebDataIn", Nexport * sizeof(struct vwebdata_in));
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_vweb(&vwebDataIn[j], place);
            memcpy(vwebDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        /* exchange particle data */
        int TAG_TO_USE = TAG_VWEBLOOP_1A;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&vwebDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct vwebdata_in), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &vwebDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct vwebdata_in), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        myfree(vwebDataIn);
        vwebDataResult = (struct vwebdata_out *) mymalloc("vwebDataResult", Nimport * sizeof(struct vwebdata_out));
        vwebDataOut = (struct vwebdata_out *) mymalloc("vwebDataOut", Nexport * sizeof(struct vwebdata_out));
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, vweb_evaluate_secondary, &threadid[j]);
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
            vweb_evaluate_secondary(&mainthreadid);
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
        TAG_TO_USE = TAG_VWEBLOOP_1B;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&vwebDataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct vwebdata_out), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &vwebDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct vwebdata_out), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_vweb(&vwebDataOut[j], place, 1);
        }
        myfree(vwebDataOut);
        myfree(vwebDataResult);
        myfree(vwebDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}


int vweb_evaluate_active_check(int i);
int vweb_evaluate_active_check(int i)
{
	// do the calculation for all particles
    if(P[i].Mass <= 0) return 0;
    return 1;
}


void *vweb_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(vweb_evaluate_active_check(i)==1)
#define EVALUATION_CALL vweb_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *vweb_evaluate_secondary(void *p)
{
#define EVALUATION_CALL vweb_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}

#endif /* VWEB */

