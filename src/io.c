#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "allvars.h"
#include "proto.h"

/*! \file io.c
 *  \brief Output of a snapshot file to disk.
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly to
 * write out new/modified quantities, as needed)
 */

static int n_type[6];
static long long ntot_type_all[6];

static int n_info;

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savepositions(int num)
{
    size_t bytes;
    char buf[500];
    int n, filenr, gr, ngroups, masterTask, lastTask;
    
    CPU_Step[CPU_MISC] += measure_time();
    
    rearrange_particle_sequence();
    /* ensures that new tree will be constructed */
    All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TreeDomainUpdateFrequency * All.TotNumPart);
    
    
    if(DumpFlag == 1)
    {
        if(ThisTask == 0)
            printf("\nwriting snapshot file #%d... \n", num);
        
        size_t MyBufferSize = All.BufferSize;
        if(!(CommBuffer = mymalloc("CommBuffer", bytes = MyBufferSize * 1024 * 1024)))
        {
            printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
            endrun(2);
        }
        
        
        if(NTask < All.NumFilesPerSnapshot)
        {
            if(ThisTask == 0)
                printf
                ("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
            endrun(0);
        }
        if(All.SnapFormat < 1 || All.SnapFormat > 3)
        {
            if(ThisTask == 0)
                printf("Unsupported File-Format\n");
            endrun(0);
        }
#ifndef  HAVE_HDF5
        if(All.SnapFormat == 3)
        {
            if(ThisTask == 0)
                printf("Code wasn't compiled with HDF5 support enabled!\n");
            endrun(0);
        }
#endif
        
        
        /* determine global and local particle numbers */
        for(n = 0; n < 6; n++)
            n_type[n] = 0;
        
        for(n = 0; n < NumPart; n++)
            n_type[P[n].Type]++;
        
        sumup_large_ints(6, n_type, ntot_type_all);
        
        /* assign processors to output files */
        distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);
        
        if(All.NumFilesPerSnapshot > 1)
        {
            if(ThisTask == 0)
            {
                sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
                mkdir(buf, 02755);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        if(All.NumFilesPerSnapshot > 1)
            sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
        else
            sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);
        
        
        ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
        if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
            ngroups++;
        
        for(gr = 0; gr < ngroups; gr++)
        {
            if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
            {
                write_file(buf, masterTask, lastTask);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        myfree(CommBuffer);
        
        if(ThisTask == 0)
            printf("done with snapshot.\n");
        
        All.Ti_lastoutput = All.Ti_Current;
        
        CPU_Step[CPU_SNAPSHOT] += measure_time();
        
    }
    
#ifdef FOF
    if(RestartFlag != 4)
    {
        if(ThisTask == 0)
            printf("\ncomputing group catalogue...\n");
        
        fof_fof(num);
        
        if(ThisTask == 0)
            printf("done with group catalogue.\n");
        
        CPU_Step[CPU_FOF] += measure_time();
    }
#endif
    
#ifdef OUTPUT_POWERSPEC
    if(RestartFlag != 4)
    {
        if(ThisTask == 0)
            printf("\ncomputing power spectra...\n");
        
        calculate_power_spectra(num, &ntot_type_all[0]);
        
        if(ThisTask == 0)
            printf("done with power spectra.\n");
        
        CPU_Step[CPU_MISC] += measure_time();
    }
#endif
    
}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
    int n, k, pindex;
    MyOutputFloat *fp;
    MyOutputPosFloat *fp_pos;
    MyIDType *ip;
    int *ip_int;
    float *fp_single;
    integertime dt_step;
    
#ifdef BOX_PERIODIC
    MyFloat boxSize;
#endif
#ifdef PMGRID
    double dt_gravkick_pm = 0;
#endif
    double dt_gravkick, dt_hydrokick;
    
#ifdef OUTPUT_COOLRATE
    double tcool, u;
#endif
    
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE))
    MyBigFloat half_kick_add[6][6];
    int l;
#endif
    
    
#ifdef MAGNETIC
    double a2_inv = All.cf_a2inv;
    double gizmo2gauss = sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam) / All.UnitMagneticField_in_gauss;
    /* NOTE: we will always work -internally- in code units where MU_0 = 1; hence the 4pi here;
     [much simpler, but be sure of your conversions!] */
#endif
    
    
#ifdef PMGRID
    if(All.ComovingIntegrationOn)
        dt_gravkick_pm =
        get_gravkick_factor(All.PM_Ti_begstep,
                            All.Ti_Current) -
        get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
    else
        dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif
    
    
    fp = (MyOutputFloat *) CommBuffer;
    fp_single = (float *) CommBuffer;
    fp_pos = (MyOutputPosFloat *) CommBuffer;
    ip = (MyIDType *) CommBuffer;
    ip_int = (int *) CommBuffer;
    
    pindex = *startindex;
    
    switch (blocknr)
    {
        case IO_POS:		/* positions */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                    {
                        fp_pos[k] = P[pindex].Pos[k];
#ifdef BOX_PERIODIC
                        boxSize = All.BoxSize;
                        if(k==0) {boxSize = boxSize_X;}
                        if(k==1) {boxSize = boxSize_Y;}
                        if(k==2) {boxSize = boxSize_Z;}
                        while(fp_pos[k] < 0) {fp_pos[k] += boxSize;}
                        while(fp_pos[k] >= boxSize) {fp_pos[k] -= boxSize;}
#endif
                    }
                    n++;
                    fp_pos += 3;
                }
            break;
            
        case IO_VEL:		/* velocities */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    dt_step = (P[pindex].TimeBin ? (((integertime) 1) << P[pindex].TimeBin) : 0);
                    
                    if(All.ComovingIntegrationOn)
                    {
                        dt_gravkick =
                        get_gravkick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
                        get_gravkick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
                        dt_hydrokick =
                        (All.Ti_Current - (P[pindex].Ti_begstep + dt_step / 2)) * All.Timebase_interval / All.cf_hubble_a;
                    }
                    else
                    {
                        dt_gravkick = dt_hydrokick =
                        (All.Ti_Current - (P[pindex].Ti_begstep + dt_step / 2)) * All.Timebase_interval;
                    }
                    
                    for(k = 0; k < 3; k++)
                    {
                        fp[k] = P[pindex].Vel[k] + P[pindex].GravAccel[k] * dt_gravkick;
                        if(P[pindex].Type == 0)
                        {
                            fp[k] += SphP[pindex].HydroAccel[k] * dt_hydrokick * All.cf_atime;
                        }
                    }
#ifdef PMGRID
                    for(k = 0; k < 3; k++)
                        fp[k] += P[pindex].GravPM[k] * dt_gravkick_pm;
#endif
                    for(k = 0; k < 3; k++)
                        fp[k] *= sqrt(All.cf_a3inv);
                    
                    n++;
                    fp += 3;
                }
            break;
            
        case IO_ID:		/* particle ID */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = P[pindex].ID;
                    n++;
                }
            break;

        case IO_CHILD_ID:		/* particle 'child' ID (for splits/mergers) */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = P[pindex].ID_child_number;
                    n++;
                }
            break;

        case IO_GENERATION_ID:	/* particle ID generation (for splits/mergers) */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip++ = P[pindex].ID_generation;
                    n++;
                }
            break;

        case IO_MASS:		/* particle mass */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].Mass;
                    n++;
                }
            break;
            
        case IO_U:			/* internal energy */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if !defined(DMANNIHILATION_DM) || defined(DMAF_INJECT_AFTER_DM_STEP)
                    *fp++ = DMAX(All.MinEgySpec, SphP[pindex].InternalEnergyPred);
#else
// Estimate the energy due to DMAF corresponding to the snapshot time
					*fp++ = DMAX(All.MinEgySpec, SphP[pindex].InternalEnergyPred + drift_DMAF(pindex, P[pindex].Ti_begstep, All.Ti_Current));				
#endif
                    n++;
                }
            break;
            
        case IO_RHO:		/* density */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Density;
                    n++;
                }
            break;
            
        case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(RT_CHEM_PHOTOION)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Ne;
                    n++;
                }
#endif
            break;
            
        case IO_NH:		/* neutral hydrogen fraction */
#if defined(COOLING) || defined(RT_CHEM_PHOTOION)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#if defined(RT_CHEM_PHOTOION)
                    *fp++ = SphP[pindex].HI;
#elif (COOL_GRACKLE_CHEMISTRY > 0)
                    *fp++ = SphP[pindex].grHI;
#else
                    double u, ne, nh0 = 0, mu = 1, temp, nHeII, nhp, nHe0, nHepp;
                    ne = SphP[pindex].Ne;
                    u = DMAX(All.MinEgySpec, SphP[pindex].InternalEnergy); // needs to be in code units
                    temp = ThermalProperties(u, SphP[pindex].Density * All.cf_a3inv, pindex, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);
                    *fp++ = (MyOutputFloat) nh0;
#endif
                    n++;
                }
#endif // #if defined(COOLING) || defined(RT_CHEM_PHOTOION)
            break;
            
        case IO_HII:		/* ionized hydrogen abundance */
#if defined(RT_CHEM_PHOTOION)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].HII;
                    n++;
                }
#endif
            break;
            
        case IO_HeI:		/* neutral Helium */
#if defined(RT_CHEM_PHOTOION_HE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].HeI;
                    n++;
                }
#endif
            break;
            
        case IO_HeII:		/* ionized Helium */
#if defined(RT_CHEM_PHOTOION_HE)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].HeII;
                    n++;
                }
#endif
            break;
            
            
        case IO_HeIII:
            break;
        case IO_H2I:
            break;
        case IO_H2II:
            break;
        case IO_HM:
            break;
        case IO_HD:
            break;
        case IO_DI:
            break;
        case IO_DII:
            break;
        case IO_HeHII:
            break;
            
            
        case IO_HSML:		/* gas kernel length */
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = PPP[pindex].Hsml;
                    n++;
                }
            break;
            
        case IO_SFR:		/* star formation rate */
#ifdef GALSF
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    /* units convert to solar masses per yr */
                    *fp++ = get_starformation_rate(pindex) * ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR));
                    n++;
                }
#endif
            break;
            
        case IO_AGE:		/* stellar formation time */
#ifdef GALSF
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].StellarAge;
                    n++;
                }
#endif
            break;
            
        case IO_OSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].IMF_NumMassiveStars;
                    n++;
                  }
#endif
            break;

        case IO_GRAINSIZE:		/* grain size */
#ifdef GRAIN_FLUID
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].Grain_Size;
                    n++;
                }
#endif
            break;
            
            
        case IO_VSTURB_DISS:
#if defined(TURB_DRIVING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].DuDt_diss;
                    n++;
                }
#endif
            break;
            
        case IO_VSTURB_DRIVE:
#if defined(TURB_DRIVING)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].DuDt_drive;
                    n++;
                }
#endif
            break;
            
        case IO_HSMS:		/* kernel length for star particles */
            break;
            
        case IO_Z:			/* gas and star metallicity */
#ifdef METALS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < NUM_METAL_SPECIES; k++)
                    {
                        fp[k] = P[pindex].Metallicity[k];
                    }
                    fp += NUM_METAL_SPECIES;
                    n++;
                }
#endif
            break;
            
        case IO_POT:		/* gravitational potential */
#if defined(OUTPUT_POTENTIAL)  || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_AND_POTENTIAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].Potential;
                    n++;
                }
#endif
            break;
            
        case IO_BH_DIST:
#ifdef BH_CALC_DISTANCES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].min_dist_to_bh;
                    n++;
                }
#endif
            break;
            
        case IO_ACCEL:		/* acceleration */
#ifdef OUTPUT_ACCELERATION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                        fp[k] = All.cf_a2inv * P[pindex].GravAccel[k];
#ifdef PMGRID
                    for(k = 0; k < 3; k++)
                        fp[k] += All.cf_a2inv * P[pindex].GravPM[k];
#endif
                    if(P[pindex].Type == 0)
                        for(k = 0; k < 3; k++)
                            fp[k] += SphP[pindex].HydroAccel[k];
                    fp += 3;
                    n++;
                }
#endif
            break;
            
        case IO_DTENTR:		/* rate of change of internal energy */
#ifdef OUTPUT_CHANGEOFENERGY
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].DtInternalEnergy;
                    n++;
                }
#endif
            break;
            
        case IO_STRESSDIAG:	/* Diagonal components of viscous shear tensor */
            break;
            
        case IO_STRESSOFFDIAG:	/* Offdiagonal components of viscous shear tensor */
            break;
            
        case IO_STRESSBULK:	/* Viscous bulk tensor */
            break;
            
        case IO_DELAYTIME:
#ifdef GALSF_SUBGRID_WINDS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].DelayTime;
                    n++;
                }
#endif
            break;
            
        case IO_SHEARCOEFF:	/* Shear viscosity coefficient */
            break;
            
        case IO_TSTP:		/* timestep  */
#ifdef OUTPUT_TIMESTEP
            
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (P[pindex].TimeBin ? (((integertime) 1) << P[pindex].TimeBin) : 0) * All.Timebase_interval;
                    n++;
                }
#endif
            break;
            
        case IO_BFLD:		/* magnetic field  */
#ifdef MAGNETIC
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                        *fp++ = (Get_Particle_BField(pindex,k) * a2_inv * gizmo2gauss);
                    n++;
                }
#endif
            break;
            
        case IO_DBDT:		/* rate of change of magnetic field  */
            break;
            
        case IO_VRMS:		/* turbulent velocity around mean */
            break;
            
        case IO_VBULK:		/* mean velocity in kernel */
            break;
            
        case IO_VTAN:		/* turbulent velocity around mean, tangential part */
            break;
            
        case IO_VRAD:		/* turbulent velocity around mean, radial part */
            break;
            
            
        case IO_TRUENGB:		/* True Number of Neighbours  */
            break;
            
            
        case IO_VDIV:		/* Divergence of Vel */
            break;
            
        case IO_VROT:		/* Velocity Curl */
            break;
            
        case IO_VORT:		/* Vorticity */
#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Vorticity[0];
                    *fp++ = SphP[pindex].Vorticity[1];
                    *fp++ = SphP[pindex].Vorticity[2];
                    n++;
                }
#endif
            break;
            
            
        case IO_IMF:		/* parameters describing the IMF  */
#ifdef GALSF_SFR_IMF_VARIATION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < N_IMF_FORMPROPS; k++)
                        fp[k] = P[pindex].IMF_FormProps[k];
                    fp += N_IMF_FORMPROPS;
                    n++;
                }
#endif
            break;
            
        case IO_COSMICRAY_ENERGY:	/* energy in the cosmic ray field  */
            break;

        case IO_COSMICRAY_KAPPA:    /* local CR diffusion constant */
            break;

        case IO_COSMICRAY_ALFVEN:    /* energy in the resonant (~gyro-radii) Alfven modes field, in the +/- (with respect to B) fields  */
            break;

            
        case IO_DIVB:		/* divergence of magnetic field  */
#ifdef MAGNETIC
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    /* divB is saved in physical units */
                    *fp++ = (SphP[pindex].divB * gizmo2gauss * (SphP[pindex].Density*All.cf_a3inv / P[pindex].Mass));
                    n++;
                }
#endif
            break;
            
        case IO_ABVC:		/* artificial viscosity of particle  */
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].alpha * SphP[pindex].alpha_limiter;
                    n++;
                }
#endif
            break;
            
            
        case IO_AMDC:		/* artificial magnetic dissipation of particle  */
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Balpha;
                    n++;
                }
#endif
            break;
            
        case IO_PHI:		/* divBcleaning fuction of particle  */
#ifdef DIVBCLEANING_DEDNER
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = (Get_Particle_PhiField(pindex) * All.cf_a3inv * gizmo2gauss);
                    n++;
                }
#endif
            break;
            
        case IO_GRADPHI:		/* divBcleaning fuction of particle  */
#ifdef DIVBCLEANING_DEDNER
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                        *fp++ = (SphP[pindex].Gradients.Phi[k] * a2_inv*a2_inv * gizmo2gauss);
                    n++;
                }
#endif
            break;
            
        case IO_ROTB:		/* rot of magnetic field  */
            break;
            
        case IO_COOLRATE:		/* current cooling rate of particle  */
#ifdef OUTPUT_COOLRATE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    double ne = SphP[pindex].Ne;
                    /* get cooling time */
                    u = SphP[pindex].InternalEnergyPred;
                    tcool = GetCoolingTime(u, SphP[pindex].Density * All.cf_a3inv, ne, pindex);
                    /* convert cooling time with current thermal energy to du/dt */
                    if(tcool != 0)
                        *fp++ = u / tcool;
                    else
                        *fp++ = 0;
                    n++;
                }
#endif // OUTPUT_COOLRATE
            break;
            
        case IO_CONDRATE:		/* current heating/cooling due to thermal conduction  */
            break;
            
        case IO_DENN:		/* density normalization factor */
            break;
            
            
        case IO_BHMASS:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = BPP(pindex).BH_Mass;
                    n++;
                }
#endif
            break;
            
        case IO_BHMASSALPHA:
#ifdef BH_ALPHADISK_ACCRETION
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = BPP(pindex).BH_Mass_AlphaDisk;
                    n++;
                }
#endif
            break;
            
        case IO_BHMDOT:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = BPP(pindex).BH_Mdot;
                    n++;
                }
#endif
            break;
            
        case IO_BHPROGS:
#ifdef BH_COUNTPROGS
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *ip_int++ = BPP(pindex).BH_CountProgs;
                    n++;
                }
#endif
            break;
            
        case IO_ACRB:
#ifdef BLACK_HOLES
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].Hsml;
                    n++;
                }
#endif
            break;
            
            
        case IO_TIDALTENSORPS:
            /* 3x3 configuration-space tidal tensor that is driving the GDE */
            break;
            
            
        case IO_GDE_DISTORTIONTENSOR:
            /* full 6D phase-space distortion tensor from GDE integration */
            break;
            
        case IO_CAUSTIC_COUNTER:
            /* caustic counter */
            break;
            
        case IO_FLOW_DETERMINANT:
            /* physical NON-CUTOFF corrected stream determinant = 1.0/normed stream density * 1.0/initial stream density */
            break;
            
        case IO_STREAM_DENSITY:
            /* physical CUTOFF corrected stream density = normed stream density * initial stream density */
            break;
            
        case IO_PHASE_SPACE_DETERMINANT:
            /* determinant of phase-space distortion tensor -> should be 1 due to Liouville theorem */
            break;
            
        case IO_ANNIHILATION_RADIATION:
            /* time integrated stream density in physical units */
            break;
            
        case IO_LAST_CAUSTIC:
            /* extensive information on the last caustic the particle has passed */
            break;
            
        case IO_SHEET_ORIENTATION:
            /* initial orientation of the CDM sheet where the particle started */
            break;
            
        case IO_INIT_DENSITY:
            /* initial stream density in physical units  */
            break;
            
        case IO_SECONDORDERMASS:
            break;

        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Abar;
                    n++;
                }
#endif
            break;

        case IO_EOSYE:
#ifdef EOS_CARRIES_YE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Ye;
                    n++;
                }
#endif
            break;

            
        case IO_EOSTEMP:
#ifdef EOS_CARRIES_TEMPERATURE
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Temperature;
                    n++;
                }
#endif
            break;
            
        case IO_PRESSURE:
#if defined(EOS_GENERAL)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = SphP[pindex].Pressure;
                    n++;
                }
#endif
            break;

            case IO_EOSCS:
#if defined(EOS_GENERAL)
            for(n = 0; n < pc; pindex++)
            if(P[pindex].Type == type)
        {
            *fp++ = SphP[pindex].SoundSpeed;
            n++;
        }
#endif
            break;

        case IO_CBE_MOMENTS:
            break;

            
        case IO_EOS_STRESS_TENSOR:
#if defined(EOS_ELASTIC)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                    {
                        int kf;
                        for(kf = 0; kf < 3; kf++)
                            fp[3*k + kf] = SphP[pindex].Elastic_Stress_Tensor[kf][k];
                    }
                    n++;
                    fp += 9;
                }
#endif
            break;

            
            case IO_EOSCOMP:
#ifdef EOS_TILLOTSON
            for(n = 0; n < pc; pindex++)
            if(P[pindex].Type == type)
        {
            *ip_int++ = SphP[pindex].CompositionType;
            n++;
        }
#endif
            break;

        case IO_PARTVEL:
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                            fp[k] = SphP[pindex].ParticleVel[k];

                    n++;
                    fp += 3;
                }
#endif
            break;

            
        case IO_RADGAMMA:
#ifdef RADTRANSFER
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < N_RT_FREQ_BINS; k++)
                        fp[k] = SphP[pindex].E_gamma[k];
                    
                    n++;
                    fp += N_RT_FREQ_BINS;
                }
#endif
            break;
            
        case IO_RAD_ACCEL:
#ifdef RT_RAD_PRESSURE_OUTPUT
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++) {fp[k] = SphP[pindex].RadAccel[k];}                    
                    n++;
                    fp += 3;
                }
#endif
            break;
            
        case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 6; k++)
                    {
                        int kf;
                        for(kf = 0; kf < N_RT_FREQ_BINS; kf++)
                            fp[N_RT_FREQ_BINS*k + kf] = SphP[pindex].ET[kf][k];
                    }
                    n++;
                    fp += 6*N_RT_FREQ_BINS;
                }
#endif
            break;
            
        case IO_DMHSML:
            break;
            
        case IO_DMDENSITY:
            break;
            
        case IO_DMVELDISP:
            break;
            
        case IO_DMHSML_V:
            break;
            
        case IO_DMDENSITY_V:
            break;
            
        case IO_CHEM:
            break;
        case IO_AGS_SOFT:		/* Adaptive Gravitational Softening: softening */
#if defined(ADAPTIVE_GRAVSOFT_FORALL) && defined(AGS_OUTPUTGRAVSOFT)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = PPP[pindex].AGS_Hsml;
                    n++;
                }
#endif
            break;
        case IO_AGS_RHO:        /* Adaptive Gravitational Softening: density */
            break;
        case IO_AGS_QPT:        /* quantum potential (Q) */
            break;
        case IO_AGS_ZETA:		/* Adaptive Gravitational Softening: zeta */
#if defined(ADAPTIVE_GRAVSOFT_FORALL) && defined(AGS_OUTPUTZETA)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = PPPZ[pindex].AGS_zeta;
                    n++;
                }
#endif
            break;
        case IO_AGS_OMEGA:		/* Adaptive Gravitational Softening: omega */
            break;
        case IO_AGS_NGBS:		/* Adaptive Gravitational Softening: neighbours */
            break;
        case IO_AGS_CORR:		/* Adaptive Gravitational Softening: correction */
            break;
            
        case IO_MG_PHI:      /* scalar field fluctuations */
            break;
            
        case IO_MG_ACCEL:      /* modified gravity acceleration */
            break;
            
            
        case IO_grHI:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHI;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grHII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHII;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grHM:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHM;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grHeI:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHeI;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grHeII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHeII;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grHeIII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHeIII;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grH2I:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grH2I;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grH2II:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grH2II;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grDI;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grDII:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grDII;
                    n++;
                }
            }
#endif
            break;
            
        case IO_grHDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            for(n = 0; n < pc; pindex++){
                if(P[pindex].Type == type){
                    *fp++ = SphP[pindex].grHDI;
                    n++;
                }
            }
#endif
            break;
            
		case IO_DENSDM:
#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].DensityDM;
                    n++;
                }
#endif
            break;

		case IO_DENSGAS:
#ifdef DMANNIHILATION_DM
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    *fp++ = P[pindex].DensAroundDM;
                    n++;
                }
#endif
            break;

		case IO_DMAF_Dtu:
#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
#ifdef DMANNIHILATION
					*fp++ = SphP[pindex].DMAF_Dtu;
#else
					*fp++ = SphP[pindex].DMAF_Dtu_tot;
#endif
                    n++;
                }
#endif
            break;


		case IO_VELSHEAR:
#ifdef VWEB
            for(n = 0; n < pc; pindex++)
                if(P[pindex].Type == type)
                {
                    for(k = 0; k < 3; k++)
                    {
                        int kf;
                        for(kf = 0; kf < 3; kf++)
                            fp[3*k + kf] = P[pindex].VelShear[kf][k];
                    }
                    n++;
                    fp += 9;
                }
#endif
			break;
            
        case IO_LASTENTRY:
            endrun(213);
            break;
    }
    
    *startindex = pindex;
}




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
    int bytes_per_blockelement = 0;
    
    switch (blocknr)
    {
        case IO_POS:
            if(mode)
                bytes_per_blockelement = 3 * sizeof(MyInputPosFloat);
            else
                bytes_per_blockelement = 3 * sizeof(MyOutputPosFloat);
            break;
            
        case IO_VEL:
        case IO_PARTVEL:
        case IO_ACCEL:
        case IO_VBULK:
        case IO_BFLD:
        case IO_GRADPHI:
        case IO_DBDT:
        case IO_ROTB:
        case IO_STRESSDIAG:
        case IO_STRESSOFFDIAG:
        case IO_RAD_ACCEL:
        case IO_VORT:
        case IO_MG_ACCEL:
            if(mode)
                bytes_per_blockelement = 3 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
            break;

        case IO_COSMICRAY_ALFVEN:
            if(mode)
                bytes_per_blockelement = 2 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 2 * sizeof(MyOutputFloat);
            break;

        case IO_ID:
            bytes_per_blockelement = sizeof(MyIDType);
            break;

        case IO_CHILD_ID:
            bytes_per_blockelement = sizeof(MyIDType);
            break;

        case IO_EOSCOMP:
        case IO_GENERATION_ID:
        case IO_BHPROGS:
        case IO_TRUENGB:
        case IO_AGS_NGBS:
            bytes_per_blockelement = sizeof(int);
            break;
            
        case IO_MASS:
        case IO_BH_DIST:
        case IO_SECONDORDERMASS:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_HeIII:
        case IO_H2I:
        case IO_H2II:
        case IO_HM:
        case IO_HD:
        case IO_DI:
        case IO_DII:
        case IO_HeHII:
        case IO_HSML:
        case IO_SFR:
        case IO_AGE:
        case IO_OSTAR:
        case IO_GRAINSIZE:
        case IO_DELAYTIME:
        case IO_HSMS:
        case IO_POT:
        case IO_DTENTR:
        case IO_STRESSBULK:
        case IO_SHEARCOEFF:
        case IO_TSTP:
        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_KAPPA:
        case IO_DIVB:
        case IO_VRMS:
        case IO_VRAD:
        case IO_VTAN:
        case IO_VDIV:
        case IO_VROT:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_COOLRATE:
        case IO_CONDRATE:
        case IO_DENN:
        case IO_BHMASS:
        case IO_BHMASSALPHA:
        case IO_ACRB:
        case IO_BHMDOT:
        case IO_CAUSTIC_COUNTER:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_EOSTEMP:
        case IO_EOSABAR:
        case IO_EOSYE:
        case IO_EOSCS:
        case IO_PRESSURE:
        case IO_INIT_DENSITY:
        case IO_AGS_SOFT:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_ZETA:
        case IO_AGS_OMEGA:
        case IO_AGS_CORR:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_MG_PHI:
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
        case IO_grH2I:
        case IO_grH2II:
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
		case IO_DENSDM:
		case IO_DENSGAS:
		case IO_DMAF_Dtu:
            if(mode)
                bytes_per_blockelement = sizeof(MyInputFloat);
            else
                bytes_per_blockelement = sizeof(MyOutputFloat);
            break;
            
            
        case IO_DMHSML:
        case IO_DMDENSITY:
        case IO_DMVELDISP:
        case IO_DMHSML_V:
        case IO_DMDENSITY_V:
            bytes_per_blockelement = sizeof(float);
            break;

            
        case IO_IMF:
#ifdef GALSF_SFR_IMF_VARIATION
            if(mode)
                bytes_per_blockelement = (N_IMF_FORMPROPS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (N_IMF_FORMPROPS) * sizeof(MyOutputFloat);
#endif
            break;
            
        case IO_RADGAMMA:
#ifdef RADTRANSFER
            if(mode)
                bytes_per_blockelement = N_RT_FREQ_BINS * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = N_RT_FREQ_BINS * sizeof(MyOutputFloat);
#endif
            break;
            
        case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
            if(mode)
                bytes_per_blockelement = (6*N_RT_FREQ_BINS) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (6*N_RT_FREQ_BINS) * sizeof(MyOutputFloat);
#endif
            break;

            
        case IO_Z:
#ifdef METALS
            if(mode)
                bytes_per_blockelement = (NUM_METAL_SPECIES) * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = (NUM_METAL_SPECIES) * sizeof(MyOutputFloat);
#endif
            break;
            
            
        case IO_EOS_STRESS_TENSOR:
            if(mode)
                bytes_per_blockelement = 9 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
            break;

            
        case IO_CBE_MOMENTS:

            
        case IO_TIDALTENSORPS:
            if(mode)
                bytes_per_blockelement = 9 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
            break;
            
        case IO_GDE_DISTORTIONTENSOR:
            if(mode)
                bytes_per_blockelement = 36 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 36 * sizeof(MyOutputFloat);
            break;
            
        case IO_ANNIHILATION_RADIATION:
            if(mode)
                bytes_per_blockelement = 3 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
            break;
            
        case IO_LAST_CAUSTIC:
            if(mode)
                bytes_per_blockelement = 20 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 20 * sizeof(MyOutputFloat);
            break;
            
        case IO_SHEET_ORIENTATION:
            if(mode)
                bytes_per_blockelement = 9 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
            break;
            
        case IO_CHEM:
            bytes_per_blockelement = 0;
            break;
            
        case IO_VELSHEAR:
            if(mode)
                bytes_per_blockelement = 9 * sizeof(MyInputFloat);
            else
                bytes_per_blockelement = 9 * sizeof(MyOutputFloat);
            break;
            
        case IO_LASTENTRY:
            endrun(214);
            break;
    }
    
    return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
    int typekey;
    
    switch (blocknr)
    {
#if defined(OUTPUT_POSITIONS_IN_DOUBLE) || defined(INPUT_POSITIONS_IN_DOUBLE)
        case IO_POS:
            typekey = 3; /* pos outputs in HDF5 are double automatically, to prevent overlaps */
            break;
#endif
            
        case IO_ID:
#ifdef LONGIDS
            typekey = 2;		/* native long long */
#else
            typekey = 0;		/* native int */
#endif
            break;

        case IO_CHILD_ID:
#ifdef LONGIDS
            typekey = 2;		/* native long long */
#else
            typekey = 0;		/* native int */
#endif
            break;

        case IO_GENERATION_ID:
        case IO_TRUENGB:
        case IO_BHPROGS:
            typekey = 0;		/* native int */
            break;
            
        default:
            typekey = 1;		/* native MyOutputFloat */
            break;
    }
    
    return typekey;
}



int get_values_per_blockelement(enum iofields blocknr)
{
    int values = 0;
    
    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_PARTVEL:
        case IO_ACCEL:
        case IO_BFLD:
        case IO_GRADPHI:
        case IO_DBDT:
        case IO_ROTB:
        case IO_STRESSDIAG:
        case IO_STRESSOFFDIAG:
        case IO_VBULK:
        case IO_RAD_ACCEL:
        case IO_VORT:
        case IO_MG_ACCEL:
            values = 3;
            break;
            
        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_MASS:
        case IO_BH_DIST:
        case IO_SECONDORDERMASS:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_HeIII:
        case IO_H2I:
        case IO_H2II:
        case IO_HM:
        case IO_HD:
        case IO_DI:
        case IO_DII:
        case IO_HeHII:
        case IO_HSML:
        case IO_SFR:
        case IO_AGE:
        case IO_OSTAR:
        case IO_GRAINSIZE:
        case IO_DELAYTIME:
        case IO_HSMS:
        case IO_POT:
        case IO_DTENTR:
        case IO_STRESSBULK:
        case IO_SHEARCOEFF:
        case IO_TSTP:
        case IO_VRMS:
        case IO_TRUENGB:
        case IO_VTAN:
        case IO_VRAD:
        case IO_VDIV:
        case IO_VROT:
        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_KAPPA:
        case IO_DIVB:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_COOLRATE:
        case IO_CONDRATE:
        case IO_DENN:
        case IO_BHMASS:
        case IO_BHMASSALPHA:
        case IO_ACRB:
        case IO_BHMDOT:
        case IO_BHPROGS:
        case IO_CAUSTIC_COUNTER:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_EOSTEMP:
        case IO_EOSABAR:
        case IO_EOSCS:
        case IO_EOSCOMP:
        case IO_EOSYE:
        case IO_PRESSURE:
        case IO_INIT_DENSITY:
        case IO_DMHSML:
        case IO_DMDENSITY:
        case IO_DMVELDISP:
        case IO_DMHSML_V:
        case IO_DMDENSITY_V:
        case IO_AGS_SOFT:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_ZETA:
        case IO_AGS_OMEGA:
        case IO_AGS_CORR:
        case IO_AGS_NGBS:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_MG_PHI:
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
        case IO_grH2I:
        case IO_grH2II:
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
		case IO_DENSDM:
		case IO_DENSGAS:
		case IO_DMAF_Dtu:
            values = 1;
            break;

        case IO_COSMICRAY_ALFVEN:
            values = 2;
            break;

        case IO_EOS_STRESS_TENSOR:
		case IO_VELSHEAR:
            values = 9;
            break;

        case IO_CBE_MOMENTS:
            values = 0;
            break;

        case IO_EDDINGTON_TENSOR:
#ifdef RADTRANSFER
            values = (6*N_RT_FREQ_BINS);
#else
            values = 0;
#endif
            break;
            
        case IO_RADGAMMA:
#ifdef RADTRANSFER
            values = N_RT_FREQ_BINS;
#else
            values = 0;
#endif
            break;
            
        case IO_Z:
#ifdef METALS
            values = NUM_METAL_SPECIES;
#else
            values = 0;
#endif
            break;
            
            
        case IO_IMF:
#ifdef GALSF_SFR_IMF_VARIATION
            values = N_IMF_FORMPROPS;
#else
            values = 0;
#endif
            break;

            
        case IO_TIDALTENSORPS:
            values = 9;
            break;
        case IO_GDE_DISTORTIONTENSOR:
            values = 36;
            break;
        case IO_ANNIHILATION_RADIATION:
            values = 3;
            break;
        case IO_LAST_CAUSTIC:
            values = 20;
            break;
        case IO_SHEET_ORIENTATION:
            values = 9;
            break;
            
        case IO_CHEM:
            values = 0;
            break;
            
        case IO_LASTENTRY:
            endrun(215);
            break;
    }
    return values;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
long get_particles_in_block(enum iofields blocknr, int *typelist)
{
    long i, nall, nsel, ntot_withmasses, ngas, nstars, nngb;
    
    nall = 0;
    nsel = 0;
    ntot_withmasses = 0;
    
    for(i = 0; i < 6; i++)
    {
        typelist[i] = 0;
        
        if(header.npart[i] > 0)
        {
            nall += header.npart[i];
            typelist[i] = 1;
        }
        
        if(All.MassTable[i] == 0)
            ntot_withmasses += header.npart[i];
    }
    
    ngas = header.npart[0];
    nstars = header.npart[4];
    
    
    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ACCEL:
        case IO_TSTP:
        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_POT:
        case IO_SECONDORDERMASS:
        case IO_AGS_SOFT:
        case IO_AGS_RHO:
        case IO_AGS_QPT:
        case IO_AGS_ZETA:
        case IO_AGS_OMEGA:
        case IO_AGS_CORR:
        case IO_AGS_NGBS:
        case IO_MG_PHI:
        case IO_BH_DIST:
        case IO_CBE_MOMENTS:
            return nall;
            break;
            
        case IO_MASS:
            for(i = 0; i < 6; i++)
            {
                typelist[i] = 0;
                if(All.MassTable[i] == 0 && header.npart[i] > 0)
                    typelist[i] = 1;
            }
            return ntot_withmasses;
            break;
            
        case IO_PARTVEL:
        case IO_RAD_ACCEL:
        case IO_RADGAMMA:
        case IO_EDDINGTON_TENSOR:
        case IO_U:
        case IO_RHO:
        case IO_NE:
        case IO_NH:
        case IO_HII:
        case IO_HeI:
        case IO_HeII:
        case IO_HeIII:
        case IO_H2I:
        case IO_H2II:
        case IO_HM:
        case IO_HD:
        case IO_DI:
        case IO_DII:
        case IO_HeHII:
        case IO_HSML:
        case IO_DELAYTIME:
        case IO_SFR:
        case IO_DTENTR:
        case IO_STRESSDIAG:
        case IO_STRESSOFFDIAG:
        case IO_STRESSBULK:
        case IO_SHEARCOEFF:
        case IO_BFLD:
        case IO_DBDT:
        case IO_VRMS:
        case IO_VBULK:
        case IO_VTAN:
        case IO_VRAD:
        case IO_VDIV:
        case IO_VROT:
        case IO_VORT:
        case IO_COSMICRAY_ENERGY:
        case IO_COSMICRAY_KAPPA:
        case IO_COSMICRAY_ALFVEN:
        case IO_DIVB:
        case IO_ABVC:
        case IO_AMDC:
        case IO_PHI:
        case IO_GRADPHI:
        case IO_ROTB:
        case IO_COOLRATE:
        case IO_CONDRATE:
        case IO_DENN:
        case IO_EOSTEMP:
        case IO_EOSABAR:
        case IO_EOSYE:
        case IO_EOSCS:
        case IO_EOS_STRESS_TENSOR:
        case IO_EOSCOMP:
        case IO_PRESSURE:
        case IO_CHEM:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
        case IO_grH2I:
        case IO_grH2II:
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
            for(i = 1; i < 6; i++)
                typelist[i] = 0;
            return ngas;
            break;
            
        case IO_AGE:
            for(i = 0; i < 6; i++)
#ifdef BLACK_HOLES
                if(i != 4 && i != 5)
                    typelist[i] = 0;
            return nstars + header.npart[5];
#else
            if(i != 4)
                typelist[i] = 0;
            return nstars;
#endif
            break;
            
        case IO_OSTAR:
            for(i = 0; i < 6; i++)
                if(i != 4)
                    typelist[i] = 0;
            return nstars;
            break;
             
        case IO_GRAINSIZE:
            nngb = 0;
            for(i = 0; i < 6; i++)
            {
                if(i == 0) typelist[i] = 0;
                if(i == 4) typelist[i] = 0;
                nngb += header.npart[i];
            }
            return nngb;
            break;


        case IO_IMF:
            for(i = 1; i < 6; i++) {if(i != 4 && i != 5) {typelist[i] = 0;}}
            return nstars + header.npart[5];
            break;
            
            
        case IO_TRUENGB:
            nngb = ngas;
            for(i = 1; i < 4; i++)
                typelist[i] = 0;
            typelist[4] = 0;
#ifdef BLACK_HOLES
            nngb += header.npart[5];
#else
            typelist[5] = 0;
#endif
            return nngb;
            break;
            
        case IO_HSMS:
            for(i = 0; i < 6; i++)
                if(i != 4)
                    typelist[i] = 0;
            return nstars;
            break;
            
        case IO_Z:
            for(i = 0; i < 6; i++)
                if(i != 0 && i != 4)
                    typelist[i] = 0;
            return ngas + nstars;
            break;
            
        case IO_BHMASS:
        case IO_BHMASSALPHA:
        case IO_ACRB:
        case IO_BHMDOT:
        case IO_BHPROGS:
            for(i = 0; i < 6; i++)
                if(i != 5)
                    typelist[i] = 0;
            return header.npart[5];
            break;
            
        case IO_TIDALTENSORPS:
        case IO_GDE_DISTORTIONTENSOR:
        case IO_CAUSTIC_COUNTER:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_ANNIHILATION_RADIATION:
        case IO_LAST_CAUSTIC:
        case IO_SHEET_ORIENTATION:
        case IO_INIT_DENSITY:
            for(i = 0; i < 6; i++)
                if(((1 << i) & (GDE_TYPES)))
                    nsel += header.npart[i];
                else
                    typelist[i] = 0;
            return nsel;
            break;
            
            
        case IO_DMHSML:
        case IO_DMDENSITY:
        case IO_DMVELDISP:
        case IO_DMHSML_V:
        case IO_DMDENSITY_V:
            for(i = 0; i < 6; i++)
                if(((1 << i) & (FOF_PRIMARY_LINK_TYPES)))
                    nsel += header.npart[i];
                else
                    typelist[i] = 0;
            return nsel;
            break;
            
        case IO_MG_ACCEL:
            for(i = 2; i < 6; i++)
                typelist[i] = 0;
            return ngas + header.npart[1];
            break;
            
		case IO_DENSDM:
			for(i = 0; i < 6; i++)
                typelist[i] = 0;
#ifdef DMANNIHILATION
			typelist[0] = 1;
			return ngas;
#elif defined(DMANNIHILATION_DM)
			typelist[1] = 1;
			return header.npart[1];
#endif			
			break; 

		case IO_DENSGAS:
			for(i = 0; i < 6; i++)
                typelist[i] = 0;
			typelist[1] = 1;
			return header.npart[1];			
			break; 

		case IO_DMAF_Dtu:
			for(i = 1; i < 6; i++)
                typelist[i] = 0;		
			return ngas;			
			break;
	
		case IO_VELSHEAR:
#ifdef VWEB
			for(i = 1; i < 6; i++)
                typelist[i] = 0;		
			return ngas;			
			break;
#endif
            
        case IO_LASTENTRY:
            endrun(216);
            break;
    }
    
    endrun(212);
    return 0;
}



/*! This function tells whether a block in the output file is present or not.
 */
int blockpresent(enum iofields blocknr)
{
    switch (blocknr)
    {
        case IO_POS:
        case IO_VEL:
        case IO_ID:
        case IO_CHILD_ID:
        case IO_GENERATION_ID:
        case IO_MASS:
        case IO_U:
        case IO_RHO:
        case IO_HSML:
            return 1;			/* always present */
            
        case IO_NE:
        case IO_NH:
#if defined(COOLING) || defined(RADTRANSFER)
            return 1;
#endif
            return 0;
            break;
            
        case IO_RADGAMMA:
#if defined(RADTRANSFER)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_RAD_ACCEL:
#if defined(RT_RAD_PRESSURE_OUTPUT)
            return 1;
#else
            return 0;
#endif
            
        case IO_HSMS:
            return 0;
            break;
            
        case IO_SFR:
        case IO_AGE:
#ifdef GALSF
            if(blocknr == IO_SFR)
                return 1;
            if(blocknr == IO_AGE)
                return 1;
#endif
            return 0;
            break;
            
        case IO_GRAINSIZE:
#ifdef GRAIN_FLUID
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_Z:
#ifdef METALS
            if(blocknr == IO_Z)
                return 1;
#endif
            return 0;
            break;
            
            
        case IO_DELAYTIME:
#ifdef GALSF_SUBGRID_WINDS
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_HeI:
        case IO_HeII:
#if defined(RT_CHEM_PHOTOION_HE)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_HII:
        case IO_HeIII:
            return 0;
            break;
            
            
            
        case IO_H2I:
        case IO_H2II:
        case IO_HM:
            return 0;
            break;
            
        case IO_HD:
        case IO_DI:
        case IO_DII:
        case IO_HeHII:
            return 0;
            break;
            
        case IO_POT:
#if defined(OUTPUT_POTENTIAL) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_AND_POTENTIAL)
            return 1;
#else
            return 0;
#endif
            
            
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
#if defined(TURB_DRIVING)
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_ACCEL:
#ifdef OUTPUT_ACCELERATION
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_DTENTR:
#ifdef OUTPUT_CHANGEOFENERGY
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_STRESSDIAG:
            return 0;
            break;
            
        case IO_STRESSOFFDIAG:
            return 0;
            break;
            
        case IO_STRESSBULK:
            return 0;
            break;
            
        case IO_SHEARCOEFF:
            return 0;
            break;
            
        case IO_TSTP:
#ifdef OUTPUT_TIMESTEP
            return 1;
#else
            return 0;
#endif
            
            
        case IO_BFLD:
#ifdef MAGNETIC
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_DBDT:
            return 0;
            break;
            
        case IO_VRMS:
        case IO_VBULK:
        case IO_TRUENGB:
            return 0;
            break;
            
        case IO_VRAD:
        case IO_VTAN:
            return 0;
            break;
            
        case IO_VDIV:
        case IO_VROT:
            return 0;
            break;
            
        case IO_VORT:
#if defined(TURB_DRIVING) || defined(OUTPUT_VORTICITY)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_IMF:
#ifdef GALSF_SFR_IMF_VARIATION
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_OSTAR:
#ifdef GALSF_SFR_IMF_SAMPLING
            return 1;
#else
            return 0;
#endif
            break;

        case IO_COSMICRAY_ENERGY:
            return 0;
            break;

        case IO_COSMICRAY_KAPPA:
            return 0;
            break;

        case IO_COSMICRAY_ALFVEN:
            return 0;
            break;

            
        case IO_DIVB:
#ifdef MAGNETIC
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_ABVC:
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_AMDC:
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_PHI:
#ifdef DIVBCLEANING_DEDNER
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_GRADPHI:
#ifdef DIVBCLEANING_DEDNER
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_ROTB:
            return 0;
            break;
            
            
        case IO_COOLRATE:
#ifdef OUTPUT_COOLRATE
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_CONDRATE:
#ifdef OUTPUTCONDRATE
            return 1;
#else
            return 0;
#endif
            break;
            
            
        case IO_DENN:
            return 0;
            break;
            
            
        case IO_ACRB:
        case IO_BHMASS:
        case IO_BHMASSALPHA:
        case IO_BHMDOT:
#ifdef BLACK_HOLES
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_BH_DIST:
#ifdef BH_CALC_DISTANCES
            return 1;
#else
            return 0;
#endif
            
        case IO_BHPROGS:
#ifdef BH_COUNTPROGS
            return 1;
#else
            return 0;
#endif
            break;
                        
            
        case IO_TIDALTENSORPS:
            return 0;
        case IO_GDE_DISTORTIONTENSOR:
            return 0;
            
        case IO_CAUSTIC_COUNTER:
            return 0;
            
        case IO_FLOW_DETERMINANT:
            return 0;
            
        case IO_STREAM_DENSITY:
            return 0;
            
        case IO_PHASE_SPACE_DETERMINANT:
            return 0;
            
        case IO_ANNIHILATION_RADIATION:
            return 0;
            
        case IO_LAST_CAUSTIC:
            return 0;
            
        case IO_SHEET_ORIENTATION:
            return 0;
            
        case IO_INIT_DENSITY:
            return 0;
            
            
        case IO_SECONDORDERMASS:
            if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
                return 1;
            else
                return 0;
            
        case IO_EOSTEMP:
#ifdef EOS_CARRIES_TEMPERATURE
            return 1;
#else
            return 0;
#endif
        case IO_PRESSURE:
#ifdef EOS_GENERAL
            return 1;
#else
            return 0;
#endif
        case IO_EOSCS:
#ifdef EOS_GENERAL
            return 1;
#else
            return 0;
#endif
        case IO_EOS_STRESS_TENSOR:
#ifdef EOS_ELASTIC
            return 1;
#else
            return 0;
#endif
        case IO_EOSCOMP:
#ifdef EOS_TILLOTSON
            return 1;
#else
            return 0;
#endif
        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            return 1;
#else
            return 0;
#endif
        case IO_EOSYE:
#ifdef EOS_CARRIES_YE
            return 1;
#else
            return 0;
#endif
            
        case IO_CBE_MOMENTS:
            return 0;

        case IO_PARTVEL:
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            return 1;
#else
            return 0;
#endif

            
        case IO_EDDINGTON_TENSOR:
#if defined(RADTRANSFER)
            return 1;
#else
            return 0;
#endif
            
        case IO_DMHSML:
        case IO_DMDENSITY:
        case IO_DMVELDISP:
            return 0;
            break;
            
        case IO_DMHSML_V:
        case IO_DMDENSITY_V:
            return 0;
            break;
            
        case IO_CHEM:
            return 0;
            break;
            
        case IO_AGS_SOFT:
#if defined (ADAPTIVE_GRAVSOFT_FORALL) && defined(AGS_OUTPUTGRAVSOFT)
            return 1;
#else
            return 0;
#endif
            break;

        case IO_AGS_RHO:
            return 0;
            break;

        case IO_AGS_QPT:
            return 0;
            break;

        case IO_AGS_ZETA:
#if defined (ADAPTIVE_GRAVSOFT_FORALL) && defined(AGS_OUTPUTZETA)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_AGS_OMEGA:
            return 0;
            break;
            
        case IO_AGS_CORR:
            return 0;
            break;
            
        case IO_AGS_NGBS:
            return 0;
            break;
            
        case IO_MG_PHI:
            return 0;
            break;
            
        case IO_MG_ACCEL:
            return 0;
            break;
            
            
        case IO_grHI:
        case IO_grHII:
        case IO_grHM:
        case IO_grHeI:
        case IO_grHeII:
        case IO_grHeIII:
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_grH2I:
        case IO_grH2II:
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            return 1;
#else
            return 0;
#endif
            break;
            
        case IO_grDI:
        case IO_grDII:
        case IO_grHDI:
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            return 1;
#else
            return 0;
#endif
            break;

		case IO_DENSDM:
#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)
			return 1;
#else
			return 0;
#endif
			break;

		case IO_DENSGAS:
#ifdef DMANNIHILATION_DM
			return 1;
#else
			return 0;
#endif  

		case IO_DMAF_Dtu:
#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)
			return 1;
#else
			return 0;
#endif          

		case IO_VELSHEAR:
#ifdef VWEB
			return 1;
#else
			return 0;
#endif 
            
        case IO_LASTENTRY:
            return 0;			/* will not occur */
    }
    
    
    return 0;			/* default: not present */
}




/*! This function associates a short 4-character block name with each block number.
 *  This is stored in front of each block for snapshot FileFormat=2.
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
    switch (blocknr)
    {
        case IO_POS:
            strncpy(label, "POS ", 4);
            break;
        case IO_VEL:
            strncpy(label, "VEL ", 4);
            break;
        case IO_CHILD_ID:
            strncpy(label, "IDch", 4);
            break;
        case IO_GENERATION_ID:
            strncpy(label, "IDgn", 4);
            break;
        case IO_ID:
            strncpy(label, "ID  ", 4);
            break;
        case IO_MASS:
            strncpy(label, "MASS", 4);
            break;
        case IO_U:
            strncpy(label, "U   ", 4);
            break;
        case IO_RHO:
            strncpy(label, "RHO ", 4);
            break;
        case IO_NE:
            strncpy(label, "NE  ", 4);
            break;
        case IO_NH:
            strncpy(label, "NH  ", 4);
            break;
        case IO_HII:
            strncpy(label, "HII ", 4);
            break;
        case IO_HeI:
            strncpy(label, "HeI ", 4);
            break;
        case IO_HeII:
            strncpy(label, "HeII", 4);
            break;
        case IO_HeIII:
            strncpy(label, "He3 ", 4);
            break;
        case IO_H2I:
            strncpy(label, "H2I ", 4);
            break;
        case IO_H2II:
            strncpy(label, "H2II", 4);
            break;
        case IO_HM:
            strncpy(label, "HM  ", 4);
            break;
        case IO_HD:
            strncpy(label, "HD  ", 4);
            break;
        case IO_DI:
            strncpy(label, "DI  ", 4);
            break;
        case IO_DII:
            strncpy(label, "DII ", 4);
            break;
        case IO_HeHII:
            strncpy(label, "HeHp", 4);
            break;
        case IO_HSML:
            strncpy(label, "HSML", 4);
            break;
        case IO_SFR:
            strncpy(label, "SFR ", 4);
            break;
        case IO_AGE:
            strncpy(label, "AGE ", 4);
            break;
        case IO_GRAINSIZE:
            strncpy(label, "GRSZ", 4);
            break;
        case IO_DELAYTIME:
            strncpy(label, "DETI", 4);
            break;
        case IO_HSMS:
            strncpy(label, "HSMS", 4);
            break;
        case IO_Z:
            strncpy(label, "Z   ", 4);
            break;
        case IO_POT:
            strncpy(label, "POT ", 4);
            break;
        case IO_ACCEL:
            strncpy(label, "ACCE", 4);
            break;
        case IO_DTENTR:
            strncpy(label, "ENDT", 4);
            break;
        case IO_STRESSDIAG:
            strncpy(label, "STRD", 4);
            break;
        case IO_STRESSOFFDIAG:
            strncpy(label, "STRO", 4);
            break;
        case IO_STRESSBULK:
            strncpy(label, "STRB", 4);
            break;
        case IO_SHEARCOEFF:
            strncpy(label, "SHCO", 4);
            break;
        case IO_TSTP:
            strncpy(label, "TSTP", 4);
            break;
        case IO_BFLD:
            strncpy(label, "BFLD", 4);
            break;
        case IO_DBDT:
            strncpy(label, "DBDT", 4);
            break;
        case IO_VBULK:
            strncpy(label, "VBLK", 4);
            break;
        case IO_VRMS:
            strncpy(label, "VRMS", 4);
            break;
        case IO_VTAN:
            strncpy(label, "VTAN", 4);
            break;
        case IO_VRAD:
            strncpy(label, "VRAD", 4);
            break;
        case IO_TRUENGB:
            strncpy(label, "TNGB", 4);
            break;
        case IO_VDIV:
            strncpy(label, "VDIV", 4);
            break;
        case IO_VROT:
            strncpy(label, "VROT", 4);
            break;
        case IO_VORT:
            strncpy(label, "VORT", 4);
            break;
        case IO_IMF:
            strncpy(label, "IMF ", 4);
            break;
        case IO_OSTAR:
            strncpy(label, "IMF ", 4);
            break;    
        case IO_COSMICRAY_ENERGY:
            strncpy(label, "CREG ", 4);
            break;
        case IO_COSMICRAY_KAPPA:
            strncpy(label, "CRK ", 4);
            break;
        case IO_COSMICRAY_ALFVEN:
            strncpy(label, "CRAV ", 4);
            break;
        case IO_DIVB:
            strncpy(label, "DIVB", 4);
            break;
        case IO_ABVC:
            strncpy(label, "ABVC", 4);
            break;
        case IO_AMDC:
            strncpy(label, "AMDC", 4);
            break;
        case IO_PHI:
            strncpy(label, "PHI ", 4);
            break;
        case IO_GRADPHI:
            strncpy(label, "GPHI", 4);
            break;
        case IO_ROTB:
            strncpy(label, "ROTB", 4);
            break;
        case IO_COOLRATE:
            strncpy(label, "COOR", 4);
            break;
        case IO_CONDRATE:
            strncpy(label, "CONR", 4);
            break;
        case IO_DENN:
            strncpy(label, "DENN", 4);
            break;
        case IO_BHMASS:
            strncpy(label, "BHMA", 4);
            break;
        case IO_BH_DIST:
            strncpy(label, "BHR ", 4);
            break;
        case IO_BHMASSALPHA:
            strncpy(label, "BHMa", 4);
            break;
        case IO_ACRB:
            strncpy(label, "ACRB", 4);
            break;
        case IO_BHMDOT:
            strncpy(label, "BHMD", 4);
            break;
        case IO_BHPROGS:
            strncpy(label, "BHPC", 4);
            break;
        case IO_TIDALTENSORPS:
            strncpy(label, "TIPS", 4);
            break;
        case IO_GDE_DISTORTIONTENSOR:
            strncpy(label, "DIPS", 4);
            break;
        case IO_CAUSTIC_COUNTER:
            strncpy(label, "CACO", 4);
            break;
        case IO_FLOW_DETERMINANT:
            strncpy(label, "FLDE", 4);
            break;
        case IO_STREAM_DENSITY:
            strncpy(label, "STDE", 4);
            break;
        case IO_SECONDORDERMASS:
            strncpy(label, "SOMA", 4);
            break;
        case IO_PHASE_SPACE_DETERMINANT:
            strncpy(label, "PSDE", 4);
            break;
        case IO_ANNIHILATION_RADIATION:
            strncpy(label, "ANRA", 4);
            break;
        case IO_LAST_CAUSTIC:
            strncpy(label, "LACA", 4);
            break;
        case IO_SHEET_ORIENTATION:
            strncpy(label, "SHOR", 4);
            break;
        case IO_INIT_DENSITY:
            strncpy(label, "INDE", 4);
            break;
        case IO_EOSTEMP:
            strncpy(label, "TEMP", 4);
            break;
        case IO_EOSABAR:
            strncpy(label, "ABAR", 4);
            break;
        case IO_EOSCS:
            strncpy(label, "EQCS", 4);
            break;
        case IO_EOS_STRESS_TENSOR:
            strncpy(label, "ESTT", 4);
            break;
        case IO_CBE_MOMENTS:
            strncpy(label, "VMOM", 4);
            break;
        case IO_EOSCOMP:
            strncpy(label, "COMP", 4);
            break;
        case IO_PARTVEL:
            strncpy(label, "PVEL", 4);
            break;
        case IO_EOSYE:
            strncpy(label, "YE  ", 4);
            break;
        case IO_PRESSURE:
            strncpy(label, "P   ", 4);
            break;
        case IO_RADGAMMA:
            strncpy(label, "RADG", 4);
            break;
        case IO_RAD_ACCEL:
            strncpy(label, "RADA", 4);
            break;
        case IO_EDDINGTON_TENSOR:
            strncpy(label, "ET", 4);
            break;
        case IO_DMHSML:
            strncpy(label, "DMHS", 4);
            break;
        case IO_DMDENSITY:
            strncpy(label, "DMDE", 4);
            break;
        case IO_DMVELDISP:
            strncpy(label, "DMVD", 4);
            break;
        case IO_DMHSML_V:
            strncpy(label, "VDMH", 4);
            break;
        case IO_DMDENSITY_V:
            strncpy(label, "VDMD", 4);
            break;
        case IO_CHEM:
            strncpy(label, "CHEM", 4);
            break;
        case IO_AGS_SOFT:
            strncpy(label, "AGSH", 4);
            break;
        case IO_AGS_RHO:
            strncpy(label, "ARHO", 4);
            break;
        case IO_AGS_QPT:
            strncpy(label, "AQPT", 4);
            break;
        case IO_AGS_ZETA:
            strncpy(label, "AGSZ", 4);
            break;
        case IO_AGS_OMEGA:
            strncpy(label, "AGSO", 4);
            break;
        case IO_AGS_CORR:
            strncpy(label, "AGSC", 4);
            break;
        case IO_AGS_NGBS:
            strncpy(label, "AGSN", 4);
            break;
        case IO_VSTURB_DISS:
            strncpy(label, "VSDI", 4);
            break;
        case IO_VSTURB_DRIVE:
            strncpy(label, "VSDR", 4);
            break;
        case IO_MG_PHI:
            strncpy(label, "MGPH", 4);
            break;
        case IO_MG_ACCEL:
            strncpy(label, "MGAC", 4);
            break;
        case IO_grHI:
            strncpy(label, "gHI", 4);
            break;
        case IO_grHII:
            strncpy(label, "gHII", 4);
            break;
        case IO_grHM:
            strncpy(label, "gHM", 4);
            break;
        case IO_grHeI:
            strncpy(label, "gHeI", 4);
            break;
        case IO_grHeII:
            strncpy(label, "gHe2", 4);
            break;
        case IO_grHeIII:
            strncpy(label, "gHe3", 4);
            break;
        case IO_grH2I:
            strncpy(label, "gH2I", 4);
            break;
        case IO_grH2II:
            strncpy(label, "H2II", 4);
            break;
        case IO_grDI:
            strncpy(label, "gDI", 4);
            break;
        case IO_grDII:
            strncpy(label, "gDII", 4);
            break;
        case IO_grHDI:
            strncpy(label, "gHDI", 4);
            break;
		case IO_DENSDM:
            strncpy(label, "dDM", 4);
            break;
		case IO_DENSGAS:
            strncpy(label, "dGAS", 4);
            break;
		case IO_DMAF_Dtu:
            strncpy(label, "DMAF", 4);
            break;
		case IO_VELSHEAR:
			strncpy(label, "VWEB", 4);
			break;
            
        case IO_LASTENTRY:
            endrun(217);
            break;
    }
}


void get_dataset_name(enum iofields blocknr, char *buf)
{
    strcpy(buf, "default");
    
    switch (blocknr)
    {
        case IO_POS:
            strcpy(buf, "Coordinates");
            break;
        case IO_VEL:
            strcpy(buf, "Velocities");
            break;
        case IO_ID:
            strcpy(buf, "ParticleIDs");
            break;
        case IO_CHILD_ID:
            strcpy(buf, "ParticleChildIDsNumber");
            break;
        case IO_GENERATION_ID:
            strcpy(buf, "ParticleIDGenerationNumber");
            break;
        case IO_MASS:
            strcpy(buf, "Masses");
            break;
        case IO_U:
            strcpy(buf, "InternalEnergy");
            break;
        case IO_RHO:
            strcpy(buf, "Density");
            break;
        case IO_NE:
            strcpy(buf, "ElectronAbundance");
            break;
        case IO_NH:
            strcpy(buf, "NeutralHydrogenAbundance");
            break;
        case IO_RADGAMMA:
            strcpy(buf, "PhotonEnergy");
            break;
        case IO_RAD_ACCEL:
            strcpy(buf, "RadiativeAcceleration");
            break;
        case IO_HII:
            strcpy(buf, "HII");
            break;
        case IO_HeI:
            strcpy(buf, "HeI");
            break;
        case IO_HeII:
            strcpy(buf, "HeII");
            break;
        case IO_HeIII:
            strcpy(buf, "HeIII");
            break;
        case IO_H2I:
            strcpy(buf, "H2I");
            break;
        case IO_H2II:
            strcpy(buf, "H2II");
            break;
        case IO_HM:
            strcpy(buf, "HM");
            break;
        case IO_HD:
            strcpy(buf, "HD  ");
            break;
        case IO_DI:
            strcpy(buf, "DI  ");
            break;
        case IO_DII:
            strcpy(buf, "DII ");
            break;
        case IO_HeHII:
            strcpy(buf, "HeHp");
            break;
        case IO_DELAYTIME:
            strcpy(buf, "DelayTime");
            break;
        case IO_HSML:
            strcpy(buf, "SmoothingLength");
            break;
        case IO_SFR:
            strcpy(buf, "StarFormationRate");
            break;
        case IO_AGE:
            strcpy(buf, "StellarFormationTime");
            break;
        case IO_GRAINSIZE:
            strcpy(buf, "GrainSize");
            break;
        case IO_HSMS:
            strcpy(buf, "StellarSmoothingLength");
            break;
        case IO_Z:
            strcpy(buf, "Metallicity");
            break;
        case IO_POT:
            strcpy(buf, "Potential");
            break;
        case IO_ACCEL:
            strcpy(buf, "Acceleration");
            break;
        case IO_DTENTR:
            strcpy(buf, "RateOfChangeOfInternalEnergy");
            break;
        case IO_STRESSDIAG:
            strcpy(buf, "DiagonalStressTensor");
            break;
        case IO_STRESSOFFDIAG:
            strcpy(buf, "OffDiagonalStressTensor");
            break;
        case IO_STRESSBULK:
            strcpy(buf, "BulkStressTensor");
            break;
        case IO_SHEARCOEFF:
            strcpy(buf, "ShearCoefficient");
            break;
        case IO_TSTP:
            strcpy(buf, "TimeStep");
            break;
        case IO_BFLD:
            strcpy(buf, "MagneticField");
            break;
        case IO_DBDT:
            strcpy(buf, "RateOfChangeOfMagneticField");
            break;
        case IO_VRMS:
            strcpy(buf, "RMSVelocity");
            break;
        case IO_VBULK:
            strcpy(buf, "BulkVelocity");
            break;
        case IO_VRAD:
            strcpy(buf, "RMSRadialVelocity");
            break;
        case IO_VTAN:
            strcpy(buf, "RMSTangentialVelocity");
            break;
        case IO_TRUENGB:
            strcpy(buf, "TrueNumberOfNeighbours");
            break;
        case IO_VDIV:
            strcpy(buf, "VelocityDivergence");
            break;
        case IO_VROT:
            strcpy(buf, "VelocityCurl");
            break;
        case IO_VORT:
            strcpy(buf, "Vorticity");
            break;
        case IO_IMF:
            strcpy(buf, "IMFFormationProperties");
            break;
        case IO_OSTAR:
            strcpy(buf, "OStarNumber");
            break;    
        case IO_COSMICRAY_ENERGY:
            strcpy(buf, "CosmicRayEnergy");
            break;
        case IO_COSMICRAY_KAPPA:
            strcpy(buf, "CosmicRayDiffusivity");
            break;
        case IO_COSMICRAY_ALFVEN:
            strcpy(buf, "CosmicRayAlfvenEnergyPM");
            break;
        case IO_DIVB:
            strcpy(buf, "DivergenceOfMagneticField");
            break;
        case IO_ABVC:
            strcpy(buf, "ArtificialViscosity");
            break;
        case IO_AMDC:
            strcpy(buf, "ArtMagneticDissipation");
            break;
        case IO_PHI:
            strcpy(buf, "DivBcleaningFunctionPhi");
            break;
        case IO_GRADPHI:
            strcpy(buf, "DivBcleaningFunctionGradPhi");
            break;
        case IO_ROTB:
            strcpy(buf, "RotationB");
            break;
        case IO_COOLRATE:
            strcpy(buf, "CoolingRate");
            break;
        case IO_CONDRATE:
            strcpy(buf, "ConductionRate");
            break;
        case IO_DENN:
            strcpy(buf, "Denn");
            break;
        case IO_BHMASS:
            strcpy(buf, "BH_Mass");
            break;
        case IO_BH_DIST:
            strcpy(buf, "BH_Dist");
            break;
        case IO_BHMASSALPHA:
            strcpy(buf, "BH_Mass_AlphaDisk");
            break;
        case IO_ACRB:
            strcpy(buf, "BH_AccretionLength");
            break;
        case IO_BHMDOT:
            strcpy(buf, "BH_Mdot");
            break;
        case IO_BHPROGS:
            strcpy(buf, "BH_NProgs");
            break;
        case IO_TIDALTENSORPS:
            strcpy(buf, "TidalTensorPS");
            break;
        case IO_GDE_DISTORTIONTENSOR:
            strcpy(buf, "DistortionTensorPS");
            break;
        case IO_CAUSTIC_COUNTER:
            strcpy(buf, "CausticCounter");
            break;
        case IO_FLOW_DETERMINANT:
            strcpy(buf, "FlowDeterminant");
            break;
        case IO_STREAM_DENSITY:
            strcpy(buf, "StreamDensity");
            break;
        case IO_SECONDORDERMASS:
            strcpy(buf, "2lpt-mass");
            break;
        case IO_PHASE_SPACE_DETERMINANT:
            strcpy(buf, "PhaseSpaceDensity");
            break;
        case IO_ANNIHILATION_RADIATION:
            strcpy(buf, "AnnihilationRadiation");
            break;
        case IO_LAST_CAUSTIC:
            strcpy(buf, "LastCaustic");
            break;
        case IO_SHEET_ORIENTATION:
            strcpy(buf, "SheetOrientation");
            break;
        case IO_INIT_DENSITY:
            strcpy(buf, "InitDensity");
            break;
        case IO_EOSTEMP:
            strcpy(buf, "Temperature");
            break;
        case IO_EOSABAR:
            strcpy(buf, "Abar");
            break;
        case IO_EOSCS:
            strcpy(buf, "SoundSpeed");
            break;
        case IO_EOS_STRESS_TENSOR:
            strcpy(buf, "StressTensor");
            break;
        case IO_CBE_MOMENTS:
            strcpy(buf, "VlasovMoments");
            break;
        case IO_EOSCOMP:
            strcpy(buf, "CompositionType");
            break;
        case IO_PARTVEL:
            strcpy(buf, "ParticleVelocities");
            break;
        case IO_EOSYE:
            strcpy(buf, "Ye");
            break;
        case IO_PRESSURE:
            strcpy(buf, "Pressure");
            break;
        case IO_EDDINGTON_TENSOR:
            strcpy(buf, "EddingtonTensor");
            break;
        case IO_DMHSML:
            strcpy(buf, "DM Hsml");
            break;
        case IO_DMDENSITY:
            strcpy(buf, "DM Density");
            break;
        case IO_DMVELDISP:
            strcpy(buf, "DM Velocity Dispersion");
            break;
        case IO_DMHSML_V:
            strcpy(buf, "DM Hsml Voronoi");
            break;
        case IO_DMDENSITY_V:
            strcpy(buf, "DM Density Voronoi");
            break;
        case IO_CHEM:
            strcpy(buf, "ChemicalAbundances");
            break;
        case IO_AGS_SOFT:
            strcpy(buf, "AGS-Softening");
            break;
        case IO_AGS_RHO:
            strcpy(buf, "AGS-Density");
            break;
        case IO_AGS_QPT:
            strcpy(buf, "AGS-QuantumPotentialQ");
            break;
        case IO_AGS_ZETA:
            strcpy(buf, "AGS-Zeta");
            break;
        case IO_AGS_OMEGA:
            strcpy(buf, "AGS-Omega");
            break;
        case IO_AGS_CORR:
            strcpy(buf, "AGS-Correction");
            break;
        case IO_AGS_NGBS:
            strcpy(buf, "AGS-Neighbours");
            break;
        case IO_VSTURB_DISS:
            strcpy(buf, "TurbulenceDissipation");
            break;
        case IO_VSTURB_DRIVE:
            strcpy(buf, "TurbulenceDriving");
            break;
        case IO_MG_PHI:
            strcpy(buf, "ModifiedGravityPhi");
            break;
        case IO_MG_ACCEL:
            strcpy(buf, "ModifiedGravityAcceleration");
            break;
        case IO_grHI:
            strcpy(buf, "GrackleHI");
            break;
        case IO_grHII:
            strcpy(buf, "GrackleHII");
            break;
        case IO_grHM:
            strcpy(buf, "GrackleHM");
            break;
        case IO_grHeI:
            strcpy(buf, "GrackleHeI");
            break;
        case IO_grHeII:
            strcpy(buf, "GrackleHeII");
            break;
        case IO_grHeIII:
            strcpy(buf, "GrackleHeIII");
            break;
        case IO_grH2I:
            strcpy(buf, "GrackleH2I");
            break;
        case IO_grH2II:
            strcpy(buf, "GrackleH2II");
            break;
        case IO_grDI:
            strcpy(buf, "GrackleDI");
            break;
        case IO_grDII:
            strcpy(buf, "GrackleDII");
            break;
        case IO_grHDI:
            strcpy(buf, "GrackleHDI");
            break;
		case IO_DENSDM:
            strcpy(buf, "DensityDM");
            break;
		case IO_DENSGAS:
            strcpy(buf, "GasDensAroundDM");
            break;
		case IO_DMAF_Dtu:
            strcpy(buf, "DMAF_Dtu");
            break;
		case IO_VELSHEAR:
			strcpy(buf, "VelocityShear");
			break;
            
        case IO_LASTENTRY:
            endrun(218);
            break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the gas kernel
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 */
void write_file(char *fname, int writeTask, int lastTask)
{
    int type, bytes_per_blockelement, npart, nextblock, typelist[6];
    int n_for_this_task, n, p, pc, offset = 0, task, i;
    size_t blockmaxlen;
    int ntot_type[6], nn[6];
    enum iofields blocknr;
    char label[8];
    int bnr;
    int blksize;
    MPI_Status status;
    FILE *fd = 0;
    
#ifdef HAVE_HDF5
    hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
    hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
    herr_t hdf5_status;
    hsize_t dims[2], count[2], start[2];
    int rank = 0, pcsum = 0;
    char buf[500];
#endif
    
#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}
    
    
    /* determine particle numbers of each type in file */
    
    if(ThisTask == writeTask)
    {
        for(n = 0; n < 6; n++)
            ntot_type[n] = n_type[n];
        
        for(task = writeTask + 1; task <= lastTask; task++)
        {
            MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
            for(n = 0; n < 6; n++)
                ntot_type[n] += nn[n];
        }
        
        for(task = writeTask + 1; task <= lastTask; task++)
            MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
        MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }
    
    /* fill file header */
    
    for(n = 0; n < 6; n++)
    {
        header.npart[n] = (int) ntot_type[n];
        header.npartTotal[n] = (unsigned int) ntot_type_all[n];
        header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }
    
    if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
        header.flag_ic_info = FLAG_EVOLVED_2LPT;
    
    if(header.flag_ic_info == FLAG_ZELDOVICH_ICS)
        header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;
    
    if(header.flag_ic_info == FLAG_NORMALICS_2LPT)
        header.flag_ic_info = FLAG_EVOLVED_2LPT;
    
    if(header.flag_ic_info == 0 && All.ComovingIntegrationOn != 0)
        header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;
    
    for(n = 0; n < 6; n++)
        header.mass[n] = All.MassTable[n];
    
    
    header.time = All.Time;
    
    if(All.ComovingIntegrationOn)
        header.redshift = 1.0 / All.Time - 1;
    else
        header.redshift = 0;
    
    header.flag_sfr = 0;
    header.flag_feedback = 0;
    header.flag_cooling = 0;
    header.flag_stellarage = 0;
    header.flag_metals = 0;
    
#ifdef COOLING
    header.flag_cooling = 1;
#endif
    
#ifdef GALSF
    header.flag_sfr = 1;
    header.flag_feedback = 1;
    header.flag_stellarage = 1;
#endif
    
#ifdef METALS
    header.flag_metals = NUM_METAL_SPECIES;
#endif
    
    header.num_files = All.NumFilesPerSnapshot;
    header.BoxSize = All.BoxSize;
    header.Omega0 = All.Omega0;
    header.OmegaLambda = All.OmegaLambda;
    header.HubbleParam = All.HubbleParam;
    
#ifdef OUTPUT_IN_DOUBLEPRECISION
    header.flag_doubleprecision = 1;
#else
    header.flag_doubleprecision = 0;
#endif
    
    /* open file and write header */
    
    if(ThisTask == writeTask)
    {
        if(All.SnapFormat == 3)
        {
#ifdef HAVE_HDF5
            sprintf(buf, "%s.hdf5", fname);
            hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            
            hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
            
            for(type = 0; type < 6; type++)
            {
                if(header.npart[type] > 0)
                {
                    sprintf(buf, "/PartType%d", type);
                    hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
                }
            }
            
            write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
        }
        else
        {
            if(!(fd = fopen(fname, "w")))
            {
                printf("can't open file `%s' for writing snapshot.\n", fname);
                endrun(123);
            }
            
            if(All.SnapFormat == 2)
            {
                blksize = sizeof(int) + 4 * sizeof(char);
                SKIP;
                my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
                nextblock = sizeof(header) + 2 * sizeof(int);
                my_fwrite(&nextblock, sizeof(int), 1, fd);
                SKIP;
            }
            
            blksize = sizeof(header);
            SKIP;
            my_fwrite(&header, sizeof(header), 1, fd);
            SKIP;
        }
    }
    
    if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
        n_info = 0;
        InfoBlock = (struct info_block *) mymalloc("InfoBlock", sizeof(struct info_block) * 1000);
        
        for(bnr = 0; bnr < 1000; bnr++)
        {
            blocknr = (enum iofields) bnr;
            
            if(blocknr == IO_SECONDORDERMASS)
                continue;
            
            if(blocknr == IO_LASTENTRY)
                break;
            
            if(blockpresent(blocknr))
            {
                bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
                npart = get_particles_in_block(blocknr, &typelist[0]);
                
                if(npart > 0)
                {
                    for(type = 0; type < 6; type++)
                        InfoBlock[n_info].is_present[type] = typelist[type];
                    InfoBlock[n_info].ndim = get_values_per_blockelement(blocknr);
                    get_Tab_IO_Label(blocknr, label);
                    for(i = 0; i < 4; i++)
                        InfoBlock[n_info].label[i] = label[i];
                    switch (get_datatype_in_block(blocknr))
                    {
                        case 0:
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "LONG    ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "LONGN   ", 8);
                            break;
                        case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
#else
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
#endif
                            break;
                        case 2:
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "LLONG   ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "LLONGN  ", 8);
                            break;
                        case 3:
#ifdef OUTPUT_POSITIONS_IN_DOUBLE
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
#else
                            if(InfoBlock[n_info].ndim <= 1)
                                strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
                            else
                                strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
#endif
                            break;
                    }
                    n_info++;
                }
            }
        }
    }
    
    for(bnr = 0; bnr < 1000; bnr++)
    {
        blocknr = (enum iofields) bnr;
        
        if(blocknr == IO_SECONDORDERMASS)
            continue;
        
        if(blocknr == IO_LASTENTRY)
            break;
        
        if(blockpresent(blocknr))
        {
            bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
            
            size_t MyBufferSize = All.BufferSize;
            blockmaxlen = (size_t) ((MyBufferSize * 1024 * 1024) / bytes_per_blockelement);
            
            npart = get_particles_in_block(blocknr, &typelist[0]);
            
            if(npart > 0)
            {
                if(ThisTask == 0)
                {
                    char buf[1000];
                    
                    get_dataset_name(blocknr, buf);
                    printf("writing block %d (%s)...\n", bnr, buf);
                }
                
                if(ThisTask == writeTask)
                {
                    
                    if(All.SnapFormat == 1 || All.SnapFormat == 2)
                    {
                        if(All.SnapFormat == 2)
                        {
                            blksize = sizeof(int) + 4 * sizeof(char);
                            SKIP;
                            get_Tab_IO_Label(blocknr, label);
                            my_fwrite(label, sizeof(char), 4, fd);
                            nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                            my_fwrite(&nextblock, sizeof(int), 1, fd);
                            SKIP;
                        }
                        blksize = npart * bytes_per_blockelement;
                        SKIP;
                    }
                }
                
                for(type = 0; type < 6; type++)
                {
                    if(typelist[type])
                    {
#ifdef HAVE_HDF5
                        if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                        {
                            switch (get_datatype_in_block(blocknr))
                            {
                                case 0:
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
                                    break;
                                case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                    break;
                                case 2:
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
                                    break;
                                case 3:
#ifdef OUTPUT_POSITIONS_IN_DOUBLE
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else 
                                    hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                    break;
                            }
                            
                            dims[0] = header.npart[type];
                            dims[1] = get_values_per_blockelement(blocknr);
                            if(dims[1] == 1)
                                rank = 1;
                            else
                                rank = 2;
                            
                            get_dataset_name(blocknr, buf);
                            
                            hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
#ifndef IO_COMPRESS_HDF5
                            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
#else
                            if(dims[0] > 10)
			    {
                            	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
                            	hsize_t cdims[2]; cdims[0] = (hsize_t) (dims[0] / 10); cdims[1] = dims[1];
                            	hdf5_status = H5Pset_chunk (plist_id, rank, cdims);
                            	hdf5_status = H5Pset_deflate (plist_id, 4);
                            	hdf5_dataset = H5Dcreate2(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT, plist_id, H5P_DEFAULT);
			    } else {
                            	hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
			    }                      
#endif
                            pcsum = 0;
                        }
#endif
                        
                        for(task = writeTask, offset = 0; task <= lastTask; task++)
                        {
                            if(task == ThisTask)
                            {
                                n_for_this_task = n_type[type];
                                
                                for(p = writeTask; p <= lastTask; p++)
                                    if(p != ThisTask)
                                        MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
                            }
                            else
                                MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
                                         &status);
                            
                            while(n_for_this_task > 0)
                            {
                                pc = n_for_this_task;
                                
                                if(pc > (int)blockmaxlen)
                                    pc = blockmaxlen;
                                
                                if(ThisTask == task)
                                    fill_write_buffer(blocknr, &offset, pc, type);
                                
                                if(ThisTask == writeTask && task != writeTask)
                                    MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
                                             TAG_PDATA, MPI_COMM_WORLD, &status);
                                
                                if(ThisTask != writeTask && task == ThisTask)
                                    MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
                                              TAG_PDATA, MPI_COMM_WORLD);
                                
                                if(ThisTask == writeTask)
                                {
                                    if(All.SnapFormat == 3)
                                    {
#ifdef HAVE_HDF5
                                        start[0] = pcsum;
                                        start[1] = 0;
                                        
                                        count[0] = pc;
                                        count[1] = get_values_per_blockelement(blocknr);
                                        pcsum += pc;
                                        
                                        H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
                                                            start, NULL, count, NULL);
                                        
                                        dims[0] = pc;
                                        dims[1] = get_values_per_blockelement(blocknr);
                                        hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);
                                        
                                        hdf5_status =
                                        H5Dwrite(hdf5_dataset, hdf5_datatype,
                                                 hdf5_dataspace_memory,
                                                 hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);
                                        
                                        H5Sclose(hdf5_dataspace_memory);
#endif
                                    }
                                    else
                                    {
                                        my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                                    }
                                }
                                
                                n_for_this_task -= pc;
                            }
                        }
                        
#ifdef HAVE_HDF5
                        if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
                        {
                            if(All.SnapFormat == 3)
                            {
                                H5Dclose(hdf5_dataset);
                                H5Sclose(hdf5_dataspace_in_file);
                                H5Tclose(hdf5_datatype);
                            }
                        }
#endif
                    }
                }
                
                if(ThisTask == writeTask)
                {
                    if(All.SnapFormat == 1 || All.SnapFormat == 2)
                    {
                        SKIP;
                    }
                }
            }
        }
    }
    
    if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
        myfree(InfoBlock);
    }
    
    if(ThisTask == writeTask)
    {
        if(All.SnapFormat == 3)
        {
#ifdef HAVE_HDF5
            for(type = 5; type >= 0; type--)
                if(header.npart[type] > 0)
                    H5Gclose(hdf5_grp[type]);
            H5Gclose(hdf5_headergrp);
            H5Fclose(hdf5_file);
#endif
        }
        else
        {
            fclose(fd);
        }
    }
}




#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle)
{
    hsize_t adim[1] = { 6 };
    
    hid_t hdf5_dataspace, hdf5_attribute;
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, header.npart);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);	
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Flag_IC_Info", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_ic_info);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
}
#endif





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
    size_t nwritten;
    
    if(size * nmemb > 0)
    {
        if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
        {
            printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
            fflush(stdout);
            endrun(777);
        }
    }
    else
        nwritten = 0;
    
    return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
    size_t nread;
    
    if(size * nmemb == 0)
        return 0;
    
    if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
        if(feof(stream))
            printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
        else
            printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
        fflush(stdout);
        endrun(778);
    }
    return nread;
}

/** This is so we don't have to preface every printf call with the if
 statement to only make it print once. */
void mpi_printf(const char *fmt, ...)
{
    if(ThisTask == 0)
    {
        va_list l;
        va_start(l, fmt);
        vprintf(fmt, l);
#ifndef IO_REDUCED_MODE
        fflush(stdout);
#endif
        va_end(l);
    }
}

#if defined(FLAG_NOT_IN_PUBLIC_CODE_READ_FOF) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
int io_compare_P_ID(const void *a, const void *b)
{
    if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID))
        return -1;
    
    if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID))
        return +1;
    
    return 0;
}

int io_compare_P_GrNr_SubNr(const void *a, const void *b)
{
    if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
        return -1;
    
    if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
        return +1;
    
    if(((struct particle_data *) a)->SubNr < (((struct particle_data *) b)->SubNr))
        return -1;
    
    if(((struct particle_data *) a)->SubNr > (((struct particle_data *) b)->SubNr))
        return +1;
    
    return 0;
}

int io_compare_P_GrNr_ID(const void *a, const void *b)
{
    if(((struct particle_data *) a)->GrNr < (((struct particle_data *) b)->GrNr))
        return -1;
    
    if(((struct particle_data *) a)->GrNr > (((struct particle_data *) b)->GrNr))
        return +1;
    
    if(((struct particle_data *) a)->ID < (((struct particle_data *) b)->ID))
        return -1;
    
    if(((struct particle_data *) a)->ID > (((struct particle_data *) b)->ID))
        return +1;
    
    return 0;
}

#endif

