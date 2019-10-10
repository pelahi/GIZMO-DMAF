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
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

/* Routines for mechanical feedback/enrichment models: stellar winds, supernovae, etc
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef GALSF_FB_MECHANICAL

// define kernel structure (purely for convenience, will hold variables below) //
struct kernel_addFB {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};

struct addFBdata_out
{
    MyFloat M_coupled, Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS];
}
*AddFBDataResult, *AddFBDataOut;

void particle2in_addFB(struct addFBdata_in *in, int i, int fb_loop_iteration);
void particle2in_addFB_wt(struct addFBdata_in *in, int i);
void out2particle_addFB(struct addFBdata_out *out, int i, int mode, int fb_loop_iteration);

void particle2in_addFB(struct addFBdata_in *in, int i, int fb_loop_iteration)
{
    // pre-assign various values that will be used regardless of feedback physics //
    int k; for(k=0;k<3;k++) {in->Pos[k]=P[i].Pos[k]; in->Vel[k]=P[i].Vel[k];}
    double heff=PPP[i].Hsml / PPP[i].NumNgb; in->V_i=heff*heff*heff; in->Hsml = PPP[i].Hsml;
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.0;}
#endif
    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {in->Area_weighted_sum[k] = P[i].Area_weighted_sum[k];}
    in->Msne = 0; in->unit_mom_SNe = 0; in->SNe_v_ejecta = 0;
    if((P[i].DensAroundStar <= 0)||(P[i].Mass == 0)) {return;} // events not possible
    if(fb_loop_iteration < 0) {in->Msne=P[i].Mass; in->unit_mom_SNe=1.e-4; in->SNe_v_ejecta=1.0e-4; return;} // weighting loop
    particle2in_addFB_fromstars(in,i,fb_loop_iteration); // subroutine that actually deals with the assignment of feedback properties
    in->unit_mom_SNe = in->Msne * in->SNe_v_ejecta;
}

void out2particle_addFB(struct addFBdata_out *out, int i, int mode, int fb_loop_iteration)
{
    if(fb_loop_iteration < 0)
    {
        int k=0, kmin=0, kmax=7; if(fb_loop_iteration == -1) {kmin=kmax; kmax=AREA_WEIGHTED_SUM_ELEMENTS;}
#ifdef GALSF_USE_SNE_ONELOOP_SCHEME
        kmin=0; kmax=AREA_WEIGHTED_SUM_ELEMENTS;
#endif
        for(k=kmin;k<kmax;k++) {ASSIGN_ADD(P[i].Area_weighted_sum[k], out->Area_weighted_sum[k], mode);}
    } else {
        P[i].Mass -= out->M_coupled; if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}
    }
}



void mechanical_fb_calc(int fb_loop_iteration)
{
    /* allocate buffers to arrange communication */
    int j, k, ngrp, ndone, ndone_flag, recvTask, place, save_NextParticle;
    long long n_exported = 0, NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct addFBdata_in) + sizeof(struct addFBdata_out) + sizemax(sizeof(struct addFBdata_in),sizeof(struct addFBdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    NextParticle = FirstActiveParticle;    /* begin with this index */
    do
    {
        BufferFullFlag = 0; Nexport = 0; save_NextParticle = NextParticle;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0; Exportflag[j] = -1;}
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
            pthread_create(&mythreads[j], &attr, addFB_evaluate_primary, &threadid[j]);
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
            addFB_evaluate_primary(&mainthreadid, fb_loop_iteration);    /* do local particles and prepare export list */
        }
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
#endif
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            NextParticle = save_NextParticle;
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle) {break;}
                if(ProcessedFlag[NextParticle] != 1) {break;}
                ProcessedFlag[NextParticle] = 2;
                NextParticle = NextActiveParticle[NextParticle];
            }
            if(NextParticle == save_NextParticle) {endrun(116608);} /* in this case, the buffer is too small to process even a single particle */
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
            } else {new_export++;}
            Nexport = new_export;
        }
        n_exported += Nexport;
        for(j = 0; j < NTask; j++) {Send_count[j] = 0;}
        for(j = 0; j < Nexport; j++) {Send_count[DataIndexTable[j].Task]++;}
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
        AddFBDataGet = (struct addFBdata_in *) mymalloc("AddFBDataGet", Nimport * sizeof(struct addFBdata_in));
        AddFBDataIn = (struct addFBdata_in *) mymalloc("AddFBDataIn", Nexport * sizeof(struct addFBdata_in));
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_addFB(&AddFBDataIn[j], place, fb_loop_iteration);
            memcpy(AddFBDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        /* exchange particle data */
        int TAG_TO_USE = TAG_FBLOOP_1A;
        if(fb_loop_iteration==-2) {TAG_TO_USE = TAG_FBLOOP_5A;}
        if(fb_loop_iteration==-1) {TAG_TO_USE = TAG_FBLOOP_1A;}
        if(fb_loop_iteration== 0) {TAG_TO_USE = TAG_FBLOOP_2A;}
        if(fb_loop_iteration== 1) {TAG_TO_USE = TAG_FBLOOP_3A;}
        if(fb_loop_iteration== 2) {TAG_TO_USE = TAG_FBLOOP_4A;}
        if(fb_loop_iteration== 3) {TAG_TO_USE = TAG_FBLOOP_5A;}
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) /* get the particles */
                {
                    MPI_Sendrecv(&AddFBDataIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct addFBdata_in), MPI_BYTE,
                                 recvTask, TAG_TO_USE, &AddFBDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct addFBdata_in), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        myfree(AddFBDataIn);
        AddFBDataResult = (struct addFBdata_out *) mymalloc("AddFBDataResult", Nimport * sizeof(struct addFBdata_out));
        AddFBDataOut = (struct addFBdata_out *) mymalloc("AddFBDataOut", Nexport * sizeof(struct addFBdata_out));
        
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_create(&mythreads[j], &attr, addFB_evaluate_secondary, &threadid[j]);}
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
            addFB_evaluate_secondary(&mainthreadid, fb_loop_iteration);
        }
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++) {pthread_join(mythreads[j], NULL);}
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        if(NextParticle < 0) {ndone_flag = 1;} else {ndone_flag = 0;}
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /* get the result */
        TAG_TO_USE = TAG_FBLOOP_1B;
        if(fb_loop_iteration==-2) {TAG_TO_USE = TAG_FBLOOP_5B;}
        if(fb_loop_iteration==-1) {TAG_TO_USE = TAG_FBLOOP_1B;}
        if(fb_loop_iteration== 0) {TAG_TO_USE = TAG_FBLOOP_2B;}
        if(fb_loop_iteration== 1) {TAG_TO_USE = TAG_FBLOOP_3B;}
        if(fb_loop_iteration== 2) {TAG_TO_USE = TAG_FBLOOP_4B;}
        if(fb_loop_iteration== 3) {TAG_TO_USE = TAG_FBLOOP_5B;}
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0) /* send the results */
                {
                    MPI_Sendrecv(&AddFBDataResult[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct addFBdata_out), MPI_BYTE, recvTask, TAG_TO_USE,
                                 &AddFBDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct addFBdata_out), MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_addFB(&AddFBDataOut[j], place, 1, fb_loop_iteration);
        }
        myfree(AddFBDataOut);
        myfree(AddFBDataResult);
        myfree(AddFBDataGet);
    }
    while(ndone < NTask);
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}




#ifdef GALSF_USE_SNE_ONELOOP_SCHEME

int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int fb_loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0;
    int j, k, n;
    double u,r2,h2;
    double v_ejecta_max,kernel_zero,wk,dM,dP;
    double E_coupled,dP_sum,dP_boost_sum;
    
    struct kernel_addFB kernel;
    struct addFBdata_in local;
    struct addFBdata_out out;
    memset(&out, 0, sizeof(struct addFBdata_out));
    
    v_ejecta_max = 5000.0 * 1.0e5/ All.UnitVelocity_in_cm_per_s;
    // 'speed limit' to prevent numerically problematic kicks at low resolution //
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1);
    
    /* Load the data for the particle injecting feedback */
    if(mode == 0)
    particle2in_addFB(&local, target, fb_loop_iteration);
    else
    local = AddFBDataGet[target];
    
    if(local.Msne<=0) return 0; // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml;
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    
    // some units (just used below, but handy to define for clarity) //
    double unitlength_in_kpc=All.UnitLength_in_cm/All.HubbleParam/3.086e21*All.cf_atime;
    double density_to_n=All.cf_a3inv*All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / PROTONMASS;
    double unit_egy_SNe = 1.0e51/(All.UnitEnergy_in_cgs/All.HubbleParam);
    
    
    // now define quantities that will be used below //
    double Esne51;
    Esne51 = 0.5*local.SNe_v_ejecta*local.SNe_v_ejecta*local.Msne / unit_egy_SNe;
    double RsneKPC, RsneKPC_0;//, RsneMAX;
    RsneKPC=0.; //RsneMAX=local.Hsml;
    RsneKPC_0=(0.0284/unitlength_in_kpc) * pow(1+Esne51,0.286); //Cioffi: weak external pressure
    double r2max_phys = 2.0/unitlength_in_kpc; // no super-long-range effects allowed! (of course this is arbitrary in code units) //
    r2max_phys *= r2max_phys;
    
    
    
    /* Now start the actual FB computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = AddFBDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            
            if(numngb_inbox < 0)
            return -1;
            
            E_coupled = dP_sum = dP_boost_sum = 0;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) continue; // same particle //
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) continue; // outside kernel (in both 'directions') //
                if(r2 > r2max_phys) continue; // outside long-range cutoff //
                // calculate kernel quantities //
                kernel.r = sqrt(r2); if(kernel.r <= 0) continue;
                u = kernel.r * kernel.hinv;
                double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j;
                double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = P[j].Mass / SphP[j].Density;
                kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 1);
                kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);
                if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}
                if(V_j<0 || isnan(V_j)) {V_j=0;}
                double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j); // effective face area //
                wk = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r))); // corresponding geometric weight //
                
                if((wk <= 0)||(isnan(wk))) continue; // no point in going further, there's no physical weight here
                double wk_vec[AREA_WEIGHTED_SUM_ELEMENTS] = {0};
                wk_vec[0] = wk;
                if(kernel.dp[0]>0) {wk_vec[1]=wk*kernel.dp[0]/kernel.r; wk_vec[2]=0;} else {wk_vec[1]=0; wk_vec[2]=wk*kernel.dp[0]/kernel.r;}
                if(kernel.dp[1]>0) {wk_vec[3]=wk*kernel.dp[1]/kernel.r; wk_vec[4]=0;} else {wk_vec[3]=0; wk_vec[4]=wk*kernel.dp[1]/kernel.r;}
                if(kernel.dp[2]>0) {wk_vec[5]=wk*kernel.dp[2]/kernel.r; wk_vec[6]=0;} else {wk_vec[5]=0; wk_vec[6]=wk*kernel.dp[2]/kernel.r;}
                
                // if fb_loop_iteration==-1, this is a pre-calc loop to get the relevant weights for coupling //
                if(fb_loop_iteration < 0)
                {
                    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) out.Area_weighted_sum[k] += wk_vec[k];
                    continue;
                }
                // NOW do the actual feedback calculation //
                double wk_norm = 1. / (MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[0])); // normalization for scalar weight sum
                wk *= wk_norm; // this way wk matches the value summed above for the weighting //
                
                if((wk <= 0)||(isnan(wk))) continue;
                
                /* define initial mass and ejecta velocity in this 'cone' */
                double v_bw[3]={0}, e_shock=0;
                double pnorm = 0;
                double pvec[3]={0};
                for(k=0; k<3; k++)
                {
                    double q; q = 0; int i1=2*k+1, i2=i1+1;
                    double q_i1 = fabs(local.Area_weighted_sum[i1]);
                    double q_i2 = fabs(local.Area_weighted_sum[i2]);
                    if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                    {
                        double rr = q_i2/q_i1;
                        double rr2 = rr * rr;
                        if(wk_vec[i1] != 0)
                        {
                            q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                        } else {
                            q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                        }
                    } else {
                        q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                    }
                    pvec[k] = -q;
                    pnorm += pvec[k]*pvec[k];
                }
                pnorm = sqrt(pnorm);
                
                wk = pnorm; // this (vector norm) is the new 'weight function' for our purposes
                dM = wk * local.Msne;
                
                /* now, add contribution from relative star-gas particle motion to shock energy */
                for(k=0;k<3;k++)
                {
                    v_bw[k] = local.SNe_v_ejecta*pvec[k]/pnorm + (local.Vel[k]-P[j].Vel[k])/All.cf_atime;
                    e_shock += v_bw[k]*v_bw[k];
                }
                double mj_preshock, dM_ejecta_in, massratio_ejecta, mu_j;
                mj_preshock = P[j].Mass;
                dM_ejecta_in = dM;
                massratio_ejecta = dM_ejecta_in / (dM_ejecta_in + P[j].Mass);
                mu_j = P[j].Mass / (dM + P[j].Mass);
                e_shock *= pnorm * 0.5*local.Msne * mu_j;
                
                if((wk <= 0)||(isnan(wk))) continue;
                
                RsneKPC = RsneKPC_0;
                double n0 = SphP[j].Density*density_to_n;
                /* this is tedious, but is a fast approximation (essentially a lookup table) for the -0.429 power above */
                if(n0 < 1.e-3) {RsneKPC *= 19.4;} else {
                    if(n0 < 1.e-2) {RsneKPC *= 1.9 + 23./(1.+333.*n0);} else {
                        if(n0 < 1.e-1) {RsneKPC *= 0.7 + 8.4/(1.+33.3*n0);} else {
                            if(n0 < 1) {RsneKPC *= 0.08 + 3.1/(1.+2.5*n0);} else {
                                if(n0 < 10) {RsneKPC *= 0.1 + 1.14/(1.+0.333*n0);} else {
                                    if(n0 < 100) {RsneKPC *= 0.035 + 0.43/(1.+0.0333*n0);} else {
                                        if(n0 < 1000) {RsneKPC *= 0.017 + 0.154/(1.+0.00333*n0);} else {
                                            if(n0 < 1.e4) {RsneKPC *= 0.006 + 0.057/(1.+0.000333*n0);} else {
                                                RsneKPC *= pow(n0, -0.429); }}}}}}}}
                
                
                /* below expression is again just as good a fit to the simulations, and much faster to evaluate */
                double z0 = P[j].Metallicity[0]/All.SolarAbundances[0];
                if(z0 < 0.01)
                {
                    RsneKPC *= 2.0;
                } else {
                    if(z0 < 1)
                    {
                        RsneKPC *= 0.93 + 0.0615 / (0.05 + 0.8*z0);
                    } else {
                        RsneKPC *= 0.8 + 0.4 / (1 + z0);
                    }
                }
                /* calculates cooling radius given density and metallicity in this annulus into which the ejecta propagate */
                
                /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium:
                 from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool */
                double r_eff_ij = sqrt(r2) - Get_Particle_Size(j);
                if(r_eff_ij > RsneKPC) {e_shock *= RsneKPC*RsneKPC*RsneKPC/(r_eff_ij*r_eff_ij*r_eff_ij);}
                
                /* now we have the proper energy to couple */
                E_coupled += e_shock;
                
                /* inject actual mass from mass return */
                if(P[j].Hsml<=0) {if(SphP[j].Density>0){SphP[j].Density*=(1+dM_ejecta_in/P[j].Mass);} else {SphP[j].Density=dM_ejecta_in*kernel.hinv3;}} else {SphP[j].Density+=kernel_zero*dM_ejecta_in*hinv3_j;}
                SphP[j].Density *= 1 + dM_ejecta_in/P[j].Mass; // inject mass at constant particle volume //
                P[j].Mass += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[j].MassTrue += dM_ejecta_in;
#endif
#ifdef METALS
                /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k]=(1-massratio_ejecta)*P[j].Metallicity[k] + massratio_ejecta*local.yields[k];}
#endif
                
                /* inject the post-shock energy and momentum (convert to specific units as needed first) */
                e_shock *= 1 / P[j].Mass;
                SphP[j].InternalEnergy += e_shock;
                SphP[j].InternalEnergyPred += e_shock;
                /* inject momentum */
                double m_ej_input = pnorm * local.Msne;
                /* appropriate factor for the ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                double m_cooling = 4.18879*pnorm*SphP[j].Density*RsneKPC*RsneKPC*RsneKPC;
                /* apply limiter for energy conservation */
                double mom_boost_fac = 1 + sqrt(DMIN(mj_preshock , m_cooling) / m_ej_input);
                if(fb_loop_iteration > 0) {mom_boost_fac=1;}
                
                /* save summation values for outputs */
                dP = local.unit_mom_SNe / P[j].Mass * pnorm;
                dP_sum += dP;
                dP_boost_sum += dP * mom_boost_fac;
                
                /* actually do the injection */
                double q0 = All.cf_atime * (pnorm*local.Msne/P[j].Mass) * mom_boost_fac;
                for(k=0; k<3; k++)
                {
                    double q = q0 * v_bw[k];
                    P[j].Vel[k] += q;
                    SphP[j].VelPred[k] += q;
                }
                
#ifdef PM_HIRES_REGION_CLIPPING
                dP=0; for(k=0;k<3;k++) dP+=P[j].Vel[k]*P[j].Vel[k]; dP=sqrt(dP);
                if(dP>5.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) P[j].Mass=0;
                if(dP>1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) for(k=0;k<3;k++) P[j].Vel[k]*=(1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s)/dP;
#endif
                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = AddFBDataGet[target].NodeList[listindex];
                if(startnode >= 0)
                startnode = Nodes[startnode].u.d.nextnode;    /* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0)
    out2particle_addFB(&out, target, 0, fb_loop_iteration);
    else
    AddFBDataResult[target] = out;
    
    return 0;
} // int addFB_evaluate



#else // un-protected [updated, more fixed energy-injecting SNe scheme]


int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int fb_loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,kernel_zero,wk,dM_ejecta_in,dP,E_coupled,dP_sum,dP_boost_sum;
    struct kernel_addFB kernel;
    struct addFBdata_in local;
    struct addFBdata_out out;
    memset(&out, 0, sizeof(struct addFBdata_out));
    
    /* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_addFB(&local, target, fb_loop_iteration);} else {local = AddFBDataGet[target];}
    if(local.Msne<=0) {return 0;} // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) {return 0;} // zero-extent kernel, no particles //
    
    // some units (just used below, but handy to define for clarity) //
    double h2 = local.Hsml*local.Hsml;
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1); // define the kernel zero-point value, needed to prevent some nasty behavior when no neighbors found
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4); // define kernel quantities
    double unitlength_in_kpc=All.UnitLength_in_cm/All.HubbleParam/3.086e21*All.cf_atime;
    double density_to_n=All.cf_a3inv*All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam / PROTONMASS;
    double unit_egy_SNe = 1.0e51/(All.UnitEnergy_in_cgs/All.HubbleParam);
    double v_ejecta_max = 5000.0 * 1.0e5/ All.UnitVelocity_in_cm_per_s; // 'speed limit' to prevent numerically problematic kicks at low resolution //
    
    // now define quantities that will be used below //
    double psi_cool=1, psi_egycon=1, v_ejecta_eff=local.SNe_v_ejecta;
    double wk_norm = 1. / (MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[0])); // normalization for scalar weight sum
    double pnorm_sum = 1./(MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[10])); // re-normalization after second pass for normalized "pnorm" (should be close to ~1)
    if((local.Area_weighted_sum[0] > MIN_REAL_NUMBER) && (fb_loop_iteration >= 0))
    {
        double vba_2_eff = wk_norm * local.Area_weighted_sum[7]; // phi term for energy: weighted mass-deposited KE for ejecta neighbors
        v_ejecta_eff = sqrt(local.SNe_v_ejecta*local.SNe_v_ejecta + vba_2_eff); // account for all terms to get the revised KE term here
        double beta_egycon = sqrt(pnorm_sum / local.Msne) * (1./v_ejecta_eff) * local.Area_weighted_sum[8]; // beta term for re-normalization for energy [can be positive or negative]
        double beta_cool = pnorm_sum * local.Area_weighted_sum[9]; // beta term if all particles in terminal-momentum-limit
        if(All.ComovingIntegrationOn) {if(fabs(beta_cool) < fabs(beta_egycon)) {beta_egycon = beta_cool;}}
        psi_egycon = sqrt(1. + beta_egycon*beta_egycon) - beta_egycon; // exact solution for energy equation for constant psi
        if(beta_egycon > 20.) {psi_egycon = 1./(2.*beta_egycon);} // replace with series expansion to avoid roundoff error at high beta
        if(beta_cool > 0.5) {psi_cool = 1./(2.*beta_cool);} // for cooling limit, only need upper limit to psi, all else will use less energy
    }
    
    
    double Energy_injected_codeunits = 0.5 * local.Msne * v_ejecta_eff * v_ejecta_eff;
    double Esne51 = Energy_injected_codeunits / unit_egy_SNe;
    double RsneKPC = 0., RsneKPC_3 = 0., m_cooling = 0., v_cooling = 2.1e7 / All.UnitVelocity_in_cm_per_s;
    double RsneKPC_0 = (0.0284/unitlength_in_kpc);
    int feedback_type_is_SNe = 0;
    if(fb_loop_iteration == 0) {feedback_type_is_SNe = 1;} // assume, for now, that loop 0 represents SNe, for purposes of energy-momentum switch below //
    if(feedback_type_is_SNe == 1) // check for SNe specifically
    {
        RsneKPC_0 *= pow(1+Esne51,0.286); //SNe: using scaling from Cioffi with weak external pressure
    } else {
        RsneKPC_0 *= pow(Esne51,0.286); // ensures smooth conservation for winds and tracers as mass-loading goes to vanishingly small values
    }
    double r2max_phys = 2.0/unitlength_in_kpc; // no super-long-range effects allowed! (of course this is arbitrary in code units) //
    if(local.Hsml >= r2max_phys) {psi_egycon=DMIN(psi_egycon,1); psi_cool=DMIN(psi_cool,1);}
    r2max_phys *= r2max_phys;
    
    
    /* Now start the actual FB computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;    /* root node */
    }
    else
    {
        startnode = AddFBDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;    /* open it */
    }
    
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -1;}
            
            E_coupled = dP_sum = dP_boost_sum = 0;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) {continue;} // require a gas particle //
                if(P[j].Mass <= 0) {continue;} // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef BOX_PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) continue; // same particle //
                
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) continue; // outside kernel (in both 'directions') //
                if(r2 > r2max_phys) continue; // outside long-range cutoff //
                // calculate kernel quantities //
                kernel.r = sqrt(r2);
                if(kernel.r <= 0) continue;
                u = kernel.r * kernel.hinv;
                double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j;
                double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = P[j].Mass / SphP[j].Density;
                kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 1);
                kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);
                if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}
                if(V_j<0 || isnan(V_j)) {V_j=0;}
                double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j); // effective face area //
                wk = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r))); // corresponding geometric weight //
                
                if((wk <= 0)||(isnan(wk))) continue; // no point in going further, there's no physical weight here
                
                double wk_vec[AREA_WEIGHTED_SUM_ELEMENTS] = {0};
                wk_vec[0] = wk;
                if(kernel.dp[0]>0) {wk_vec[1]=wk*kernel.dp[0]/kernel.r; wk_vec[2]=0;} else {wk_vec[1]=0; wk_vec[2]=wk*kernel.dp[0]/kernel.r;}
                if(kernel.dp[1]>0) {wk_vec[3]=wk*kernel.dp[1]/kernel.r; wk_vec[4]=0;} else {wk_vec[3]=0; wk_vec[4]=wk*kernel.dp[1]/kernel.r;}
                if(kernel.dp[2]>0) {wk_vec[5]=wk*kernel.dp[2]/kernel.r; wk_vec[6]=0;} else {wk_vec[5]=0; wk_vec[6]=wk*kernel.dp[2]/kernel.r;}
                
                RsneKPC = RsneKPC_0;
                /* calculate cooling radius given density and metallicity in this annulus into which the ejecta propagate */
                if(fb_loop_iteration < 2)
                {
                    double e0 = Esne51;
                    if(fb_loop_iteration < 0) {e0=1;}
                    if(feedback_type_is_SNe == 1) {e0+=1;}
                    double n0 = SphP[j].Density*density_to_n;
                    if(n0 < 0.001) {n0=0.001;}
                    double z0 = P[j].Metallicity[0]/All.SolarAbundances[0], z0_term = 1.;
                    if(z0 < 0.01) {z0 = 0.01;}
                    if(z0 < 1.) {z0_term = z0*sqrt(z0);} else {z0_term = z0;}
                    double nz_dep  = pow(n0 * z0_term , 0.14);;
                    v_cooling = 2.10e7 * DMAX(nz_dep,0.5) / All.UnitVelocity_in_cm_per_s;
                    m_cooling = 4.56e36 * e0 / (nz_dep*nz_dep * All.UnitMass_in_g/All.HubbleParam);
                    RsneKPC = pow( 0.238732 * m_cooling/SphP[j].Density , 1./3. );
                }
                RsneKPC_3 = RsneKPC*RsneKPC*RsneKPC;
                
                // if fb_loop_iteration==-1, this is a pre-calc loop to get the relevant weights for coupling //
                if(fb_loop_iteration < 0)
                {
                    if(fb_loop_iteration==-1) // the Area_weighted_sum quantities are computed on loop=-2; these quantities must be computed on loop=-1 (after Area_weighted_sums are computed)
                    {
                        /* calculate the corrected momentum vectors that we will actually use in the coupling proper */
                        double pnorm=0, pvec[3]={0}, vel_ba_2=0, cos_vel_ba_pcoupled=0;
                        for(k=0;k<3;k++)
                        {
                            double q = 0; int i1=2*k+1, i2=i1+1;
                            double q_i1 = fabs(local.Area_weighted_sum[i1]);
                            double q_i2 = fabs(local.Area_weighted_sum[i2]);
                            if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                            {
                                double rr = q_i2/q_i1;
                                double rr2 = rr * rr;
                                if(wk_vec[i1] != 0)
                                {
                                    q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                                } else {
                                    q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                                }
                            } else {
                                q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                            }
                            pvec[k] = -q;
                            pnorm += pvec[k]*pvec[k];
                        }
                        pnorm = sqrt(pnorm);
                        /* now calculate the additional weights that are needed for energy terms */
                        for(k=0;k<3;k++)
                        {
                            double v_ba = (P[j].Vel[k] - local.Vel[k]) / All.cf_atime; // relative gas-star velocity //
                            vel_ba_2 += v_ba*v_ba; // magnitude of velocity vector (for corrected post-shock energies to distribute)
                            cos_vel_ba_pcoupled += v_ba * pvec[k]/pnorm; // direction of ejecta [after correction loop]
                        }
                        wk_vec[7] = wk * vel_ba_2; // phi_0 term : residual KE term from mass-coupling for {small, second-order} energy correction
                        wk_vec[8] = sqrt(pnorm * P[j].Mass) * cos_vel_ba_pcoupled; // beta_0 term : cross-term for momentum coupling effect on energy-coupling
                        wk_vec[9] = pnorm * cos_vel_ba_pcoupled / v_cooling; // calculate the beta term as if all particles hit terminal: more accurate result in that limit
                        wk_vec[10] = pnorm; // normalization (so that we can divide by its sum to properly normalize the beta_egy and beta_cool quantities)
                    }
                    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {out.Area_weighted_sum[k] += wk_vec[k];}
                    continue;
                }
                // NOW do the actual feedback calculation //
                wk *= wk_norm; // this way wk matches the value summed above for the weighting //
                
                if((wk <= 0)||(isnan(wk))) continue;
                
                /* define initial mass and ejecta velocity in this 'cone' */
                double pnorm = 0, pvec[3] = {0};
                for(k=0; k<3; k++)
                {
                    double q = 0; int i1=2*k+1, i2=i1+1;
                    double q_i1 = fabs(local.Area_weighted_sum[i1]);
                    double q_i2 = fabs(local.Area_weighted_sum[i2]);
                    if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                    {
                        double rr = q_i2/q_i1;
                        double rr2 = rr * rr;
                        if(wk_vec[i1] != 0)
                        {
                            q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                        } else {
                            q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                        }
                    } else {
                        q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                    }
                    pvec[k] = -q;
                    pnorm += pvec[k]*pvec[k];
                }
                pnorm = sqrt(pnorm); // this (vector norm) is the new 'weight function' for our purposes
                dM_ejecta_in = pnorm * local.Msne;
                double mj_preshock, massratio_ejecta;
                mj_preshock = P[j].Mass;
                massratio_ejecta = dM_ejecta_in / (dM_ejecta_in + P[j].Mass);
                
                /* inject actual mass from mass return */
                if(P[j].Hsml<=0) {if(SphP[j].Density>0){SphP[j].Density*=(1+dM_ejecta_in/P[j].Mass);} else {SphP[j].Density=dM_ejecta_in*kernel.hinv3;}} else {SphP[j].Density+=kernel_zero*dM_ejecta_in*hinv3_j;}
                SphP[j].Density *= 1 + dM_ejecta_in/P[j].Mass; // inject mass at constant particle volume //
                P[j].Mass += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                SphP[j].MassTrue += dM_ejecta_in;
#endif
#ifdef METALS
                /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES;k++) {P[j].Metallicity[k]=(1-massratio_ejecta)*P[j].Metallicity[k] + massratio_ejecta*local.yields[k];}
#endif
                /* inject momentum: account for ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                double wk_m_cooling = pnorm * m_cooling; // effective cooling mass for this particle
                double boost_max = sqrt(1 + wk_m_cooling / dM_ejecta_in); // terminal momentum boost-factor
                double boost_egycon = sqrt(1 + mj_preshock / dM_ejecta_in); // energy-conserving limit for coupling through neighbors
                double mom_boost_fac = 1;
                if(feedback_type_is_SNe == 1)
                {
                    double psi0 = 1; // factor to use below for velocity-limiter
                    boost_max *= psi_cool; // appropriately re-weight boost to avoid energy conservation errors [cooling-limit]
                    boost_egycon *= psi_egycon; // appropriately re-weight boost to avoid energy conservation errors [energy-conserving-limit]
                    if((wk_m_cooling < mj_preshock) || (boost_max < boost_egycon)) {mom_boost_fac=boost_max; psi0=DMAX(psi0,psi_cool);} else {mom_boost_fac=boost_egycon; psi0=DMAX(psi0,psi_egycon);} // limit to cooling case if egy-conserving exceeds terminal boost, or coupled mass short of cooling mass
                    if(mom_boost_fac < 1) {mom_boost_fac=1;} // impose lower limit of initial ejecta momentum
                    // finally account for simple physical limiter: if particle moving away faster than cooling terminal velocity, can't reach that velocity //
                    double vcool = DMIN(v_cooling/psi0 , v_ejecta_eff/mom_boost_fac); // effective velocity at stalling/cooling radius
                    double dv_dp_phys = 0; for(k=0;k<3;k++) {dv_dp_phys += (1-massratio_ejecta) * (kernel.dp[k]/kernel.r) * ((local.Vel[k] - P[j].Vel[k])/All.cf_atime);} // recession velocity of particle from SNe
                    double v_cooling_lim = DMAX( vcool , dv_dp_phys ); // cooling vel can't be smaller than actual vel (note: negative dvdp here automatically returns vcool, as desired)
                    double boostfac_max = DMIN(1000. , v_ejecta_eff/v_cooling_lim); // boost factor cant exceed velocity limiter - if recession vel large, limits boost
                    if(mom_boost_fac > boostfac_max) {mom_boost_fac = boostfac_max;} // apply limiter
                } else {
                    mom_boost_fac = DMIN(boost_egycon , boost_max); // simply take minimum - nothing fancy for winds
                }
                
                /* save summation values for outputs */
                dP = local.unit_mom_SNe / P[j].Mass * pnorm;
                dP_sum += dP; dP_boost_sum += dP * mom_boost_fac;
                
                /* actually do the injection */
                double mom_prefactor =  mom_boost_fac * massratio_ejecta * (All.cf_atime*v_ejecta_eff) / pnorm; // this gives the appropriately-normalized tap-able momentum from the energy-conserving solution
                double KE_initial = 0, KE_final = 0;
                for(k=0; k<3; k++)
                {
                    double d_vel = mom_prefactor * pvec[k] + massratio_ejecta*(local.Vel[k] - P[j].Vel[k]); // local.Vel term from extra momentum of moving star, P[j].Vel term from going from momentum to velocity boost with added mass
                    KE_initial += P[j].Vel[k]*P[j].Vel[k]; P[j].Vel[k] += d_vel; SphP[j].VelPred[k] += d_vel; KE_final += P[j].Vel[k]*P[j].Vel[k];
                }
                /* now calculate the residual energy and add it as thermal */
                KE_initial *= 0.5 * mj_preshock * All.cf_a2inv;
                KE_final *= 0.5 * P[j].Mass * All.cf_a2inv;
                double E_sne_initial = pnorm * Energy_injected_codeunits;
                double d_Egy_internal = KE_initial + E_sne_initial - KE_final;
                if(d_Egy_internal < 0.5*E_sne_initial) {d_Egy_internal = 0.5*E_sne_initial;}
                /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium: from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool */
                double r_eff_ij = kernel.r - Get_Particle_Size(j);
                if(r_eff_ij > RsneKPC) {d_Egy_internal *= RsneKPC_3 / (r_eff_ij*r_eff_ij*r_eff_ij);}
                d_Egy_internal /= P[j].Mass; // convert to specific internal energy, finally //
                if(d_Egy_internal > 0) {SphP[j].InternalEnergy += d_Egy_internal; SphP[j].InternalEnergyPred += d_Egy_internal; E_coupled += d_Egy_internal;}
                
#ifdef PM_HIRES_REGION_CLIPPING
                double dP=0; for(k=0;k<3;k++) dP+=P[j].Vel[k]*P[j].Vel[k]; dP=sqrt(dP);
                if(dP>5.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) P[j].Mass=0;
                if(dP>1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s) for(k=0;k<3;k++) P[j].Vel[k]*=(1.e9*All.cf_atime/All.UnitVelocity_in_cm_per_s)/dP;
#endif
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = AddFBDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}    /* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addFB(&out, target, 0, fb_loop_iteration);} else {AddFBDataResult[target] = out;}
    
    return 0;
} // int addFB_evaluate

#endif // GALSF_USE_SNE_ONELOOP_SCHEME else


int addFB_evaluate_active_check(int i, int fb_loop_iteration);
int addFB_evaluate_active_check(int i, int fb_loop_iteration)
{
    if(P[i].Type <= 1) return 0;
    if(P[i].Mass <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(PPP[i].NumNgb <= 0) return 0;
    if(P[i].SNe_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==0) return 1;}
    return 0;
}


void *addFB_evaluate_primary(void *p, int fb_loop_iteration)
{
#define CONDITION_FOR_EVALUATION if(addFB_evaluate_active_check(i,fb_loop_iteration)==1)
#define EVALUATION_CALL addFB_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist, fb_loop_iteration)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *addFB_evaluate_secondary(void *p, int fb_loop_iteration)
{
#define EVALUATION_CALL addFB_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist, fb_loop_iteration);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}




void determine_where_SNe_occur(void)
{
    if(All.Time<=0) return;
    int i; double dt,star_age,npossible,nhosttotal,ntotal,ptotal,dtmean,rmean;
    npossible=nhosttotal=ntotal=ptotal=dtmean=rmean=0;
    double mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean;
    mpi_npossible=mpi_nhosttotal=mpi_ntotal=mpi_ptotal=mpi_dtmean=mpi_rmean=0;
    // loop over particles //
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        P[i].SNe_ThisTimeStep=0;
        if(All.ComovingIntegrationOn) {if(P[i].Type != 4) {continue;}} // in cosmological simulations, 'stars' have particle type=4
        if(All.ComovingIntegrationOn==0) {if((P[i].Type<2)||(P[i].Type>4)) {continue;}} // in non-cosmological sims, types 2,3,4 are valid 'stars'
        if(P[i].Mass<=0) {continue;}
#ifndef WAKEUP
        dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; // dloga to dt_physical
#else
        dt = P[i].dt_step * All.Timebase_interval / All.cf_hubble_a; // get particle timestep //
#endif
        if(dt<=0) {continue;} // no time, no events
        star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age<=0) {continue;} // unphysical age, no events
        // now use a calculation of mechanical event rates to determine where/when the events actually occur //
        npossible++;
        double RSNe = mechanical_fb_calculate_eventrates(i,dt);
        rmean += RSNe; ptotal += RSNe * dt * P[i].Mass * (All.UnitTime_in_Megayears/All.HubbleParam) * (All.UnitMass_in_g/All.HubbleParam)/SOLAR_MASS;
#ifdef GALSF_SFR_IMF_SAMPLING
        if(P[i].IMF_NumMassiveStars>0) {P[i].IMF_NumMassiveStars=DMAX(0,P[i].IMF_NumMassiveStars-P[i].SNe_ThisTimeStep);} // lose an O-star for every SNe //
#endif
        if(P[i].SNe_ThisTimeStep>0) {ntotal+=P[i].SNe_ThisTimeStep; nhosttotal++;}
        dtmean += dt;
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
    
    MPI_Reduce(&dtmean, &mpi_dtmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rmean, &mpi_rmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ptotal, &mpi_ptotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nhosttotal, &mpi_nhosttotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ntotal, &mpi_ntotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&npossible, &mpi_npossible, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
#ifdef IO_REDUCED_MODE
        if(mpi_ntotal > 0 && mpi_nhosttotal > 0 && mpi_dtmean > 0)
#endif
        if(mpi_npossible>0)
        {
            mpi_dtmean /= mpi_npossible; mpi_rmean /= mpi_npossible;
            fprintf(FdSneIIHeating, "%lg %g %g %g %g %g %g \n", All.Time,mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean);
        }
#ifdef IO_REDUCED_MODE
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
        {fflush(FdSneIIHeating);}
    } // if(ThisTask == 0) //
    
} // void determine_where_SNe_occur() //


#endif /* GALSF_FB_MECHANICAL */

