#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly initializing 
 * new/modified variables, as needed)
 */

/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the initial gas kernel lengths are determined.
 */
void init(void)
{
    int i, j;
    double a3, atime;
    
#ifdef MAGNETIC
    double a2_fac;
    double gauss2gizmo = All.UnitMagneticField_in_gauss / sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam);
    /* NOTE: we will always work -internally- in code units where MU_0 = 1; hence the 4pi here;
        [much simpler, but be sure of your conversions!] */
#endif
    
#ifdef BLACK_HOLES
    int count_holes = 0;
#endif
    
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();    
    
    if(RestartFlag == 3 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");
        endrun(0);
    }
    
    if(RestartFlag == 4 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf("Need to give the snapshot number if snapshot should be converted\n");
        endrun(0);
    }
    
    if(RestartFlag == 5 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf
            ("Need to give the snapshot number if power spectrum and two-point correlation function should be calculated\n");
        endrun(0);
    }
    
    if(RestartFlag == 6 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf
            ("Need to give the snapshot number if velocity power spectrum for the gas cells should be calculated\n");
        endrun(0);
    }
    
    
    switch (All.ICFormat)
    {
        case 1:
        case 2:
        case 3:
        case 4:
            if(RestartFlag >= 2 && RestartSnapNum >= 0)
            {
                char fname[1000];
                
                if(All.NumFilesPerSnapshot > 1)
                    sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase,
                            RestartSnapNum);
                else
                    sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
                
                read_ic(fname);
                
            }
            else
            {
                read_ic(All.InitCondFile);
            }
            break;
            
        default:
            if(ThisTask == 0)
                printf("ICFormat=%d not supported.\n", All.ICFormat);
            endrun(0);
    }
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();
    
    
    
#ifdef COOLING
    IonizeParams();
#endif
    
    if(All.ComovingIntegrationOn)
    {
        All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
        All.Ti_Current = 0;
        a3 = All.Time * All.Time * All.Time;
        atime = All.Time;
#ifdef MAGNETIC
        a2_fac = (All.Time * All.Time);
#endif
    }
    else
    {
        All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
        All.Ti_Current = 0;
        a3 = 1;
        atime = 1;
#ifdef MAGNETIC
        a2_fac = 1;
#endif
    }
        
    set_softenings();
    
    All.NumCurrentTiStep = 0;	/* setup some counters */
    All.SnapshotFileCount = 0;
    if(RestartFlag == 2)
    {
        if(RestartSnapNum < 0)
        {
            char *underscore = strrchr(All.InitCondFile, '_');
            if(!underscore)
            {
                char buf[1000];
                sprintf(buf, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n",
                        All.InitCondFile);
                terminate(buf);
            }
            else
                All.SnapshotFileCount = atoi(underscore + 1) + 1;
        }
        else
            All.SnapshotFileCount = RestartSnapNum + 1;
    }
    
#ifdef OUTPUT_LINEOFSIGHT
    All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
    if(RestartFlag == 2)
        endrun(78787);
#endif
    
    All.TotNumOfForces = 0;
    All.TopNodeAllocFactor = 0.008; /* this will start from a low value and be iteratively increased until it is well-behaved */
    All.TreeAllocFactor = 0.45; /* this will also iteratively increase to fit the particle distribution */
    /* To construct the BH-tree for N particles, somewhat less than N
     internal tree-nodes are necessary for ‘normal’ particle distributions. 
     TreeAllocFactor sets the number of internal tree-nodes allocated in units of the particle number. 
     By experience, space for ≃ 0.65N internal nodes is usually fully sufficient for typical clustered 
     particle distributions, so a value of 0.7 should put you on the safe side. If the employed particle 
     number per processor is very small (less than a thousand or so), or if there are many particle pairs 
     with identical or nearly identical coordinates, a higher value may be required. Since the number of 
     particles on a given processor may be higher by a factor PartAllocFactor than the average particle 
     number, the total amount of memory requested for the BH tree on a single processor scales proportional 
     to PartAllocFactor*TreeAllocFactor. */
    
    
    
#ifdef BOX_PERIODIC
    if(All.ComovingIntegrationOn) check_omega();
#endif
    
    All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if defined(BLACK_HOLES) || defined(GALSF_SUBGRID_WINDS)
    All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif
    
    for(i = 0; i < GRAVCOSTLEVELS; i++)
        All.LevelToTimeBin[i] = 0;
    
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < GRAVCOSTLEVELS; j++)
            P[i].GravCost[j] = 0;
    
   
    if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
        for(i = 0; i < NumPart; i++)
            for(j = 0; j < 3; j++)
	    {
                P[i].Vel[j] *= sqrt(All.Time) * All.Time;
	    }
    }
    
#ifdef DM_SIDM
    init_self_interactions();
#endif

// INJECT ENERGY
#ifdef INJECT_ENERGY
	// All.EnergyRad = 0.02; /* initialize EnergyRad */		
#endif
    
    for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
        for(j = 0; j < 3; j++)
            P[i].GravAccel[j] = 0;
        
        /* DISTORTION PARTICLE SETUP */
        
#ifdef KEEP_DM_HSML_AS_GUESS
        if(RestartFlag != 1)
            P[i].DM_Hsml = -1;
#endif
        
#ifdef PMGRID
        for(j = 0; j < 3; j++)
            P[i].GravPM[j] = 0;
#endif
        P[i].Ti_begstep = 0;
        P[i].Ti_current = 0;
        P[i].TimeBin = 0;
        
        if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
            P[i].OldAcc = 0;	/* Do not zero in 2lpt case as masses are stored here */
        
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)    
        P[i].Potential = 0;
#endif
#ifdef GALSF
        if(RestartFlag == 0)
        {
            P[i].StellarAge = 0;
#ifdef GALSF_SFR_IMF_VARIATION
            P[i].IMF_Mturnover = 2.0; /* gives a solar-type IMF for our calculations in current code */
#endif
#ifdef GALSF_SFR_IMF_SAMPLING
            P[i].IMF_NumMassiveStars = 0;
#endif
        }
#endif
        
        if(RestartFlag != 1)
        {
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES)
            P[i].DensAroundStar = 0;
            P[i].GradRho[0]=0;
            P[i].GradRho[1]=0;
            P[i].GradRho[2]=1;
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
            P[i].SNe_ThisTimeStep = 0;
#endif
#ifdef GALSF_FB_MECHANICAL
            int k; for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {P[i].Area_weighted_sum[k] = 0;}
#endif
        }
        
#if defined(FLAG_NOT_IN_PUBLIC_CODE_X) || defined(GALSF_FB_MECHANICAL) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF_FB_THERMAL)
        if(RestartFlag == 0)
        {
            P[i].StellarAge = -2.0 * All.InitStellarAgeinGyr / (All.UnitTime_in_Megayears*0.001) * get_random_number(P[i].ID + 3);
        }
#endif
        
#ifdef GRAIN_FLUID
        if(RestartFlag == 0)
        {
            /* Change grain mass to change the distribution of sizes.  Grain_Size_Spectrum_Powerlaw parameter sets d\mu/dln(R_d) ~ R_d^Grain_Size_Spectrum_Powerlaw */
            P[i].Grain_Size = All.Grain_Size_Min * exp( gsl_rng_uniform(random_generator) * log(All.Grain_Size_Max/All.Grain_Size_Min) );
            if(P[i].Type==3) {if(All.Grain_Size_Max > All.Grain_Size_Min*1.0001 && fabs(All.Grain_Size_Spectrum_Powerlaw) != 0) {P[i].Mass *= (All.Grain_Size_Spectrum_Powerlaw/(pow(All.Grain_Size_Max/All.Grain_Size_Min,All.Grain_Size_Spectrum_Powerlaw)-1.)) * pow(P[i].Grain_Size/All.Grain_Size_Min,All.Grain_Size_Spectrum_Powerlaw) * log(All.Grain_Size_Max/All.Grain_Size_Min);}}


#ifdef GRAIN_RDI_TESTPROBLEM
	    if(P[i].Type == 3) /* initialize various quantities for test problems from parameters set in the ICs */
	    {
		P[i].Mass *= All.Dust_to_Gas_Mass_Ratio;
		{ 	
			double tS0 = 0.626657 * P[i].Grain_Size * sqrt(GAMMA); /* stopping time [Epstein] for driftvel->0 */
            double a0 = tS0 * All.Vertical_Grain_Accel / (1.+All.Dust_to_Gas_Mass_Ratio); /* acc * tS0 / (1+mu) */
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
			a0 *= All.Grain_Size_Max / P[i].Grain_Size;
#endif
            double ct = cos(All.Vertical_Grain_Accel_Angle * M_PI/180.), st = sin(All.Vertical_Grain_Accel_Angle * M_PI/180.); /* relevant angles */
			int k; double agamma=0.220893; // 9pi/128 //
			double tau2=0, ct2=0, w0=sqrt((sqrt(1.+4.*agamma*a0*a0)-1.)/(2.*agamma)); // exact solution if no Lorentz forces and Epstein drag //
#ifdef GRAIN_LORENTZFORCE
			double tL_i = All.Grain_Charge_Parameter/All.Grain_Size_Max * pow(All.Grain_Size_Max/P[i].Grain_Size,2) * All.BiniZ; // 1/Lorentz in code units
            ct2=ct*ct; double tau2_0=pow(tS0*tL_i,2), f_tau_guess2=0; // variables for below //
			for(k=0;k<20;k++)
			{
			   tau2 = tau2_0 / (1. + agamma*w0*w0); // guess tau [including velocity dependence] //
			   f_tau_guess2 = (1.+tau2*ct2) / (1.+tau2); // what the projection factor (reduction in w from projection) would be //
			   w0 = sqrt((sqrt(1.+4.*agamma*a0*a0*f_tau_guess2)-1.)/(2.*agamma)); // re-calculate w0 with this // 
			}
#endif
		w0 /= sqrt((1.+tau2)*(1.+tau2*ct2)); // ensures normalization to unity with convention below //
        int non_gdir=1; 
        if(GRAV_DIRECTION_RDI==1) {non_gdir=2;}
		P[i].Vel[0] = w0*st; P[i].Vel[non_gdir] = w0*sqrt(tau2)*st; P[i].Vel[GRAV_DIRECTION_RDI] = w0*(1.+tau2)*ct;
        a0 = tS0 * All.Vertical_Gravity_Strength / (1.+All.Dust_to_Gas_Mass_Ratio); w0=sqrt((sqrt(1.+4.*agamma*a0*a0)-1.)/(2.*agamma));
        P[i].Vel[GRAV_DIRECTION_RDI] -= w0;
		}
	    }	    
#endif

            P[i].Gas_Density = P[i].Gas_InternalEnergy = P[i].Gas_Velocity[0]=P[i].Gas_Velocity[1]=P[i].Gas_Velocity[2]=0;
#ifdef GRAIN_LORENTZFORCE
            P[i].Gas_B[0]=P[i].Gas_B[1]=P[i].Gas_B[2]=0;
#endif
        }
#endif
        
        
        
#ifdef METALS
        All.SolarAbundances[0]=0.02;        // all metals (by mass); present photospheric abundances from Asplund et al. 2009 (Z=0.0134, proto-solar=0.0142) in notes;
                                            //   also Anders+Grevesse 1989 (older, but hugely-cited compilation; their Z=0.0201, proto-solar=0.0213)
#ifdef COOL_METAL_LINES_BY_SPECIES
        if (NUM_METAL_SPECIES>=10) {
            All.SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
            All.SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)
            All.SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)
            All.SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)
            All.SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
            All.SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
            All.SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)
            All.SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)
            All.SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
            All.SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)
        }
#endif // COOL_METAL_LINES_BY_SPECIES
        
        if(RestartFlag == 0) {
#if defined(COOL_METAL_LINES_BY_SPECIES) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF_FB_MECHANICAL) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GALSF_FB_THERMAL)
            P[i].Metallicity[0] = All.InitMetallicityinSolar*All.SolarAbundances[0];
#else
            P[i].Metallicity[0] = 0;
#endif
            /* initialize abundance ratios. for now, assume solar */
            for(j=0;j<NUM_METAL_SPECIES;j++) P[i].Metallicity[j]=All.SolarAbundances[j]*P[i].Metallicity[0]/All.SolarAbundances[0];
            /* need to allow for a primordial He abundance */
            if(NUM_METAL_SPECIES>=10) P[i].Metallicity[1]=0.25+(All.SolarAbundances[1]-0.25)*P[i].Metallicity[0]/All.SolarAbundances[0];
        } // if(RestartFlag == 0)
#endif // METALS                        
        
#ifdef BLACK_HOLES
        if(P[i].Type == 5)
        {
            count_holes++;
            if(RestartFlag == 0)
            {
                BPP(i).BH_Mass = All.SeedBlackHoleMass;
#ifdef BH_ALPHADISK_ACCRETION
                BPP(i).BH_Mass_AlphaDisk = All.SeedAlphaDiskMass;
#endif
#ifdef BH_COUNTPROGS
                BPP(i).BH_CountProgs = 1;
#endif
            }
        }
#endif

// DM ANNIHILATION AT DM
#ifdef DMANNIHILATION_DM
		if(RestartFlag == 0)
		{
			if(P[i].Type == 1)
        	{
		    	P[i].HsmlDM = 0;
#ifndef DMAF_DENSITY_WEIGHTED
				P[i].AreaSumDM = 0;
#endif				
		    	P[i].DensityDM = -1;
				P[i].DensAroundDM = -1;		
				P[i].DMAF_MaxTimebin = TIMEBINS;	
#ifndef DMAF_KEEP_DM_MASS
				P[i].DMAF_dM = 0;
#endif				
			}			
		}	
#endif

// VWEB
#ifdef VWEB
		if(RestartFlag == 0)
        {
			int k;
			for(j = 0; j < 3; j++)
				for(k = 0; k < 3; k++)
					P[i].VelShear[j][k] = 0;
		}			
#endif
    }
    
#ifdef BLACK_HOLES
    MPI_Allreduce(&count_holes, &All.TotBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    
    for(i = 0; i < TIMEBINS; i++)
        TimeBinActive[i] = 1;
    
    reconstruct_timebins();
    
#ifdef PMGRID
    All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif
        
    for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        
        for(j = 0; j < 3; j++)
        {
            SphP[i].VelPred[j] = P[i].Vel[j];
            SphP[i].HydroAccel[j] = 0;
            //SphP[i].dMomentum[j] = 0;//manifest-indiv-timestep-debug//
        }
        
        //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
        P[i].Particle_DivVel = 0;
        SphP[i].ConditionNumber = 1;
        SphP[i].DtInternalEnergy = 0;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        SphP[i].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[i].dMass = 0;
        SphP[i].DtMass = 0;
        SphP[i].MassTrue = P[i].Mass;
        for(j=0;j<3;j++) SphP[i].GravWorkTerm[j] = 0;
#endif
        
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[i].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        PPP[i].AGS_Hsml = PPP[i].Hsml;
#endif
#endif
        
#ifdef CONDUCTION
        SphP[i].Kappa_Conduction = 0;
#endif
#ifdef MHD_NON_IDEAL
        SphP[i].Eta_MHD_OhmicResistivity_Coeff = 0;
        SphP[i].Eta_MHD_HallEffect_Coeff = 0;
        SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
#ifdef VISCOSITY
        SphP[i].Eta_ShearViscosity = 0;
        SphP[i].Zeta_BulkViscosity = 0;
#endif
        
        
#ifdef TURB_DIFFUSION
        SphP[i].TD_DiffCoeff = 0;
#endif
        
        if(RestartFlag == 0)
        {
#ifndef INPUT_READ_HSML
            PPP[i].Hsml = 0;
#endif
            SphP[i].Density = -1;

// DM ANNIHILATION at GAS
#ifdef DMANNIHILATION
	        P[i].HsmlDM = 0;
	        P[i].DensityDM = -1;
			SphP[i].DMAF_Dtu = 0;
#endif

// DM ANNIHILATION at DM
#ifdef DMANNIHILATION_DM
			for(j = 0; j < TIMEBINS; j++)
	        	SphP[i].DMAF_Dtu[j] = 0;	
			SphP[i].DMAF_Dtu_tot = 0;
#endif

#ifdef COOLING
            SphP[i].Ne = 1.0;
#endif
        }
#ifdef GALSF_SUBGRID_WINDS
        if(RestartFlag == 0)
            SphP[i].DelayTime = 0;
#if (GALSF_SUBGRID_WIND_SCALING==1)
        SphP[i].HostHaloMass = 0;
#endif
#endif // GALSF_SUBGRID_WINDS //
#ifdef GALSF_FB_TURNOFF_COOLING
        SphP[i].DelayTimeCoolingSNe = 0;
#endif
#ifdef GALSF
        SphP[i].Sfr = 0;
#if (GALSF_SFR_VIRIAL_SF_CRITERION==3)
	SphP[i].AlphaVirial_SF_TimeSmoothed = 0;
#endif
#endif
#ifdef MAGNETIC
#if defined MHD_B_SET_IN_PARAMS
        if(RestartFlag == 0)
        {			/* Set only when starting from ICs */
            SphP[i].B[0]=SphP[i].BPred[0] = All.BiniX;
            SphP[i].B[1]=SphP[i].BPred[1] = All.BiniY;
            SphP[i].B[2]=SphP[i].BPred[2] = All.BiniZ;
        }
#endif /*MHD_B_SET_IN_PARAMS*/
        for(j = 0; j < 3; j++)
        {
            SphP[i].BPred[j] *= a2_fac * gauss2gizmo;
            SphP[i].B[j] = SphP[i].BPred[j];
        }
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
        SphP[i].Balpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
        SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        SphP[i].alpha = 0.0;
#endif
#if defined(BH_THERMALFEEDBACK)
        SphP[i].Injected_BH_Energy = 0;
#endif
    }
    
#ifndef BOX_SHEARING
#if (NUMDIMS==2)
    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[2] = 0;
        //P[i].Vel[2] = 0; // this should be set in the ICs, not here //
        
        P[i].GravAccel[2] = 0;
        
        if(P[i].Type == 0)
        {
            SphP[i].VelPred[2] = 0;
            SphP[i].HydroAccel[2] = 0;
        }
    }
#endif
#endif
    
#if (NUMDIMS==1)
    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[1] = P[i].Pos[2] = 0;
        //P[i].Vel[1] = P[i].Vel[2] = 0; // this should be set in the ICs, not here //
        
        P[i].GravAccel[1] = P[i].GravAccel[2] = 0;
        
        if(P[i].Type == 0)
        {
            SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
            SphP[i].HydroAccel[1] = SphP[i].HydroAccel[2] = 0;
        }
    }
#endif
    
#ifdef ASSIGN_NEW_IDS
    assign_unique_ids();
#endif
    /* assign other ID parameters needed */
    if(RestartFlag==0) {for(i = 0; i < NumPart; i++) {P[i].ID_child_number = 0; P[i].ID_generation = 0;}}
#ifdef NO_CHILD_IDS_IN_ICS
    if(RestartFlag != 1) {for(i = 0; i < NumPart; i++) {P[i].ID_child_number = 0; P[i].ID_generation = 0;}}
#endif
    
    
#ifdef TEST_FOR_IDUNIQUENESS
    test_id_uniqueness();
#endif
    
    Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */
    
    TreeReconstructFlag = 1;
    
    
#ifdef SHIFT_BY_HALF_BOX
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
        {
            double boxtmp = 0;
            if(j==0) {boxtmp = boxSize_X;}
            if(j==1) {boxtmp = boxSize_Y;}
            if(j==2) {boxtmp = boxSize_Z;}
            P[i].Pos[j] += 0.5 * boxtmp;
        }
#endif
    
    
    Gas_split = 0;
#ifdef GALSF
    Stars_converted = 0;
#endif
    domain_Decomposition(0, 0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */
    
    set_softenings();
    
    /* will build tree */
    ngb_treebuild();
    
    All.Ti_Current = 0;
    
    if(RestartFlag != 3 && RestartFlag != 5)
        setup_smoothinglengths();
    
// DM ANNIHILATION
#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)
	if(RestartFlag != 3 && RestartFlag != 5)
		dm_setup_smoothinglengths();
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(RestartFlag != 3 && RestartFlag != 5)
        ags_setup_smoothinglengths();
#endif
    
#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
    if(RestartFlag != 3 && RestartFlag != 5)
        disp_setup_smoothinglengths();
#endif
#endif
    
#if defined GALSF_SFR_IMF_VARIATION
    for(i = 0; i < NumPart; i++)
    {
        P[i].IMF_Mturnover = 2.0; // reset to normal IMF
    }
#endif
    
#if defined(WAKEUP) && defined(ADAPTIVE_GRAVSOFT_FORALL)
    for(i=0;i<NumPart;i++) {P[i].wakeup=0;}
#endif

#if defined(TURB_DRIVING)
    {
        double mass = 0, glob_mass;
        int i;
        for(i=0; i< N_gas; i++)
            mass += P[i].Mass;
        MPI_Allreduce(&mass, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        All.RefDensity = glob_mass / (boxSize_X*boxSize_Y*boxSize_Z);
        All.RefInternalEnergy = All.IsoSoundSpeed*All.IsoSoundSpeed / (GAMMA*GAMMA_MINUS1);
    }
#endif
    
    /* HELLO! This here is where you should insert custom code for hard-wiring the ICs of various test problems */
    
    
    density();
    for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
        int k; k=0;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;

#if defined(TURB_DRIVING) && defined(EOS_ENFORCE_ADIABAT)
        SphP[i].InternalEnergy = All.RefInternalEnergy;
        SphP[i].InternalEnergyPred = All.RefInternalEnergy;
#endif
        // re-match the predicted and initial velocities and B-field values, just to be sure //
        for(j=0;j<3;j++) SphP[i].VelPred[j]=P[i].Vel[j];
        
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==0)
        for(j=0;j<3;j++) {SphP[i].ParticleVel[k] = 0;} // set these to zero and forget them, for the rest of the run //
#endif
        
#ifdef MAGNETIC
        for(j=0;j<3;j++) {SphP[i].B[j] = SphP[i].BPred[j] * P[i].Mass / SphP[i].Density;} // convert to the conserved unit V*B //
        for(j=0;j<3;j++) {SphP[i].BPred[j]=SphP[i].B[j]; SphP[i].DtB[j]=0;}
#endif
#if defined(EOS_ELASTIC)
        if(RestartFlag != 1)
        {
            for(k=0;k<3;k++) {for(j=0;j<3;j++) {SphP[i].Dt_Elastic_Stress_Tensor[j][k] = SphP[i].Elastic_Stress_Tensor_Pred[j][k] = SphP[i].Elastic_Stress_Tensor[j][k] = 0;}}
        } else {
            for(k=0;k<3;k++) {for(j=0;j<3;j++) {SphP[i].Elastic_Stress_Tensor_Pred[j][k] = SphP[i].Elastic_Stress_Tensor[j][k]; SphP[i].Dt_Elastic_Stress_Tensor[j][k] = 0;}}
        }
#endif
        //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
        SphP[i].DtInternalEnergy = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[i].dMass = 0;
        SphP[i].DtMass = 0;
        SphP[i].MassTrue = P[i].Mass;
        for(j=0;j<3;j++) SphP[i].GravWorkTerm[j] = 0;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[i].AGS_zeta = 0;
#endif
#ifdef WAKEUP
        if(RestartFlag!=0) {PPPZ[i].wakeup=0;}
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
        SphP[i].Super_Timestep_Dt_Explicit = 0;
        SphP[i].Super_Timestep_j = 0;
#endif
        
#ifdef COOL_GRACKLE
        if(RestartFlag == 0)
        {
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            SphP[i].grHI    = HYDROGEN_MASSFRAC;
            SphP[i].grHII   = 1.0e-20;
            SphP[i].grHM    = 1.0e-20;
            SphP[i].grHeI   = 1.0 - HYDROGEN_MASSFRAC;
            SphP[i].grHeII  = 1.0e-20;
            SphP[i].grHeIII = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            SphP[i].grH2I   = 1.0e-20;
            SphP[i].grH2II  = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            SphP[i].grDI    = 2.0 * 3.4e-5;
            SphP[i].grDII   = 1.0e-20;
            SphP[i].grHDI   = 1.0e-20;
#endif
        }
#endif
        
    }
    
    
    /* we should define the maximum and minimum particle masses 
        below/above which particles are merged/split */
    if(RestartFlag != 1)
    {
        double mass_min = MAX_REAL_NUMBER;
        double mass_max = -MAX_REAL_NUMBER;
        for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
        {
            if(P[i].Mass > mass_max) mass_max = P[i].Mass;
            if(P[i].Mass < mass_min) mass_min = P[i].Mass;
        }
        /* broadcast this and get the min and max values over all processors */
        double mpi_mass_min,mpi_mass_max;
        MPI_Allreduce(&mass_min, &mpi_mass_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&mass_max, &mpi_mass_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        All.MinMassForParticleMerger = 0.49 * mpi_mass_min;
#ifdef GALSF_GENERATIONS
        All.MinMassForParticleMerger /= (float)GALSF_GENERATIONS;
#endif
        /* All.MaxMassForParticleSplit  = 5.01 * mpi_mass_max; */
        All.MaxMassForParticleSplit  = 3.01 * mpi_mass_max;
#ifdef MERGESPLIT_HARDCODE_MAX_MASS
        All.MaxMassForParticleSplit = MERGESPLIT_HARDCODE_MAX_MASS;
#endif
#ifdef MERGESPLIT_HARDCODE_MIN_MASS
        All.MinMassForParticleMerger = MERGESPLIT_HARDCODE_MIN_MASS;
#endif
    }
    
    
#ifdef PM_HIRES_REGION_CLIPDM
    if(RestartFlag != 1)
    {
        double mpi_m_hires_max, m_hires_max=0.0;
        for(i=0; i<NumPart; i++) {if(P[i].Type==1) {if(P[i].Mass > m_hires_max) {m_hires_max=P[i].Mass;}}}
        MPI_Allreduce(&m_hires_max, &mpi_m_hires_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        All.MassOfClippedDMParticles = mpi_m_hires_max;
    }
#endif
    
    
    if(RestartFlag == 3)
    {
        
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(ThisTask == 0)
            printf("*ADAPTIVE_GRAVSOFT_FORALL* Computation of softening lengths... \n");
        ags_setup_smoothinglengths();
        if(ThisTask == 0)
            printf("*ADAPTIVE_GRAVSOFT_FORALL* Computation of softening lengths done. \n");
#endif
        
#ifdef FOF
        fof_fof(RestartSnapNum);
#endif
        endrun(0);
    }
    
#ifdef OUTPUT_TWOPOINT_ENABLED
    if(RestartFlag == 5)
    {
        /* calculating powerspec and twopoint function */
#ifdef PMGRID
        long_range_init_regionsize();
#ifdef BOX_PERIODIC
        int n, n_type[6];
        long long ntot_type_all[6];
        /* determine global and local particle numbers */
        for(n = 0; n < 6; n++)
            n_type[n] = 0;
        for(n = 0; n < NumPart; n++)
            n_type[P[n].Type]++;
        sumup_large_ints(6, n_type, ntot_type_all);
        
        calculate_power_spectra(RestartSnapNum, ntot_type_all);
#endif
#endif
        force_treebuild(NumPart, NULL);
        twopoint();
        endrun(0);
    }
#endif
    
    
    if(RestartFlag == 4)
    {
        All.Time = All.TimeBegin = header.time;
        sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
        if(ThisTask == 0)
            printf("Start writing file %s\n", All.SnapshotFileBase);
        printf("RestartSnapNum %d\n", RestartSnapNum);
        
        All.TopNodeAllocFactor = 0.008;
        
#ifdef VWEB
		density();
		vweb_calc();
#endif

        savepositions(RestartSnapNum);
        endrun(0);
    }
}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
#ifdef BOX_PERIODIC
void check_omega(void)
{
    double mass = 0, masstot, omega;
    int i;
    
    for(i = 0; i < NumPart; i++)
        mass += P[i].Mass;
    
    MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    omega = masstot / (boxSize_X*boxSize_Y*boxSize_Z) / (3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G));
#ifdef GR_TABULATED_COSMOLOGY_G
    omega *= All.Gini / All.G;
#endif
    
    //if(fabs(omega - All.Omega0) > 1.0e-3)
    // because of how we set up these ICs, allow a little more generous tolerance
    if(fabs(omega - All.Omega0) > 1.0e-2)
    {
        if(ThisTask == 0)
        {
            printf("\n\nI've found something odd!\n");
            printf
            ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
             omega, All.Omega0);
            printf("\nI better stop.\n");
            
            fflush(stdout);
        }
        endrun(1);
    }
}
#endif


/*! This function is used to find an initial kernel length (what used to be called the 
 *  'smoothing length' for SPH, but is just the kernel size for the mesh-free methods) for each gas
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the kernel length is provided to the function density(), which will
 *  then iterate if needed to find the right kernel length.
 */
void setup_smoothinglengths(void)
{
    int i, no, p;
    if((RestartFlag == 0)||(RestartFlag==2)) // best for stability if we re-calc Hsml for snapshot restarts //
    {
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES) || defined(GRAIN_FLUID) || defined(DMANNIHILATION_DM)
        for(i = 0; i < NumPart; i++)
#else
            for(i = 0; i < N_gas; i++)
#endif
            {
                no = Father[i];
                
                while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    
                    if(p < 0)
                        break;
                    
                    no = p;
                }
                
                if((RestartFlag == 0)||(P[i].Type != 0)) // if Restartflag==2, use the saved Hsml of the gas as initial guess //
                {
                    
#ifndef INPUT_READ_HSML
#if NUMDIMS == 3
                    PPP[i].Hsml = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.333333) * Nodes[no].len;
#endif
#if NUMDIMS == 2
                    PPP[i].Hsml = pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.5) * Nodes[no].len;
#endif
#if NUMDIMS == 1
                    PPP[i].Hsml = All.DesNumNgb * (P[i].Mass / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#ifndef SELFGRAVITY_OFF
                    if(All.SofteningTable[0] != 0)
                    {
                        if((PPP[i].Hsml>100.*All.SofteningTable[0])||(PPP[i].Hsml<=0.01*All.SofteningTable[0])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                            PPP[i].Hsml = All.SofteningTable[0];
                    }
#else
                    if((Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) PPP[i].Hsml = 1.0;
#endif
#endif // INPUT_READ_HSML
                } // closes if((RestartFlag == 0)||(P[i].Type != 0))
            }
    }
    
    
#ifdef BLACK_HOLES
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
            if(P[i].Type == 5)
                PPP[i].Hsml = All.SofteningTable[5];
    }
#endif
    
#ifdef GRAIN_FLUID
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
            if(P[i].Type > 0)
                PPP[i].Hsml = All.SofteningTable[P[i].Type];
    }
#endif

    density();    
}

#ifdef DMANNIHILATION_DM
void init_DMAF(void)
{		
	// Calculate DM annihilation
	DMAF_calc(0); // Calculate area sum of neighbouring gas particles and find minimum allowed time steps for gas
	find_timesteps(1); // Update time steps for the gas in such way that gas time step is always smaller or equal the time step of injecting DM particles
	DMAF_calc(1); // Calculate energy injection from DMAF at DM particles into gas particles
	find_timesteps(1); // Update time steps again accounting for the energy that is going to be injected	
}
#endif


void assign_unique_ids(void)
{
    int i, *numpartlist;
    MyIDType idfirst;
    
    numpartlist = (int *) mymalloc("numpartlist", NTask * sizeof(int));
    
    MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MPI_COMM_WORLD);
    
    idfirst = 1;
    
    for(i = 0; i < ThisTask; i++)
        idfirst += numpartlist[i];
    
    for(i = 0; i < NumPart; i++)
    {
        P[i].ID = idfirst;
        idfirst++;
    }
    
    myfree(numpartlist);
}


#ifdef ADAPTIVE_GRAVSOFT_FORALL
void ags_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            P[i].Particle_DivVel = 0;
            PPPZ[i].AGS_zeta = 0;
            if(P[i].Type > 0)
            {
                no = Father[i];
                while(10 * All.AGS_DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0)
                        break;
                    no = p;
                }
                PPP[i].AGS_Hsml = 2. * pow(1.0/NORM_COEFF * All.AGS_DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                if(All.SofteningTable[P[i].Type] != 0)
                {
                    if((PPP[i].AGS_Hsml>ADAPTIVE_GRAVSOFT_FORALL*All.SofteningTable[P[i].Type])||(PPP[i].AGS_Hsml<=0.01*All.SofteningTable[P[i].Type])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                        PPP[i].AGS_Hsml = sqrt(ADAPTIVE_GRAVSOFT_FORALL) * All.SofteningTable[P[i].Type];
                }
            } else {
                PPP[i].AGS_Hsml = PPP[i].Hsml;
            }
        }
    }
    ags_density();
}
#endif // ADAPTIVE_GRAVSOFT_FORALL


#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==2)
void disp_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            if(P[i].Type == 0)
            {
                no = Father[i];
                while(10 * 2.0 * 64 * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0)
                        break;
                    no = p;
                }
                SphP[i].HsmlDM = pow(1.0/NORM_COEFF * 2.0 * 64 * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                if(All.SofteningTable[P[i].Type] != 0)
                {
                    if((SphP[i].HsmlDM >1000.*All.SofteningTable[P[i].Type])||(PPP[i].Hsml<=0.01*All.SofteningTable[P[i].Type])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                        SphP[i].HsmlDM = All.SofteningTable[P[i].Type];
                }
            }
        }
    }
    
    if(ThisTask == 0)
    {
        printf("computing DM Vel_disp around gas particles.\n");
    }
    disp_density();
}
#endif
#endif

// DM ANNIHILATION
#if defined(DMANNIHILATION) || defined(DMANNIHILATION_DM)
void dm_setup_smoothinglengths(void)
{
    int i;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            P[i].HsmlDM = PPP[i].Hsml; // guess that the dm smoothing lengths are initially the same as gas smoothing length
        }
    }
    
    if(ThisTask == 0)
    {
        printf("Initializing smoothing lengths for DM density calculations...\n");
    }
	
	dm_density();
#ifdef DMANNIHILATION_DM
	gas_density_for_DM();  /* computes gas density around DM particles */
#endif
	
}
#endif


void test_id_uniqueness(void)
{
    double t0, t1;
#ifndef BOX_BND_PARTICLES
    int i;
    MyIDType *ids, *ids_first;
#endif
    
    if(ThisTask == 0)
    {
        printf("Testing ID uniqueness...\n");
    }
    
    if(NumPart == 0)
    {
        printf("need at least one particle per cpu\n");
        endrun(8);
    }
    
    t0 = my_second();
    
#ifndef BOX_BND_PARTICLES
    ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
    ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));
    
    for(i = 0; i < NumPart; i++)
        ids[i] = P[i].ID;
    
#ifdef ALTERNATIVE_PSORT
    init_sort_ID(ids, NumPart);
#else
    parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);
#endif
    
    for(i = 1; i < NumPart; i++)
        if(ids[i] == ids[i - 1])
        {
#ifdef LONGIDS
            printf("non-unique ID=%d%09d found on task=%d (i=%d NumPart=%d)\n",
                   (int) (ids[i] / 1000000000), (int) (ids[i] % 1000000000), ThisTask, i, NumPart);
            
#else
            printf("non-unique ID=%d found on task=%d   (i=%d NumPart=%d)\n", (int) ids[i], ThisTask, i, NumPart);
#endif
            endrun(12);
        }
    
    MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);
    
    if(ThisTask < NTask - 1)
        if(ids[NumPart - 1] == ids_first[ThisTask + 1])
        {
            printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
            endrun(13);
        }
    
    myfree(ids_first);
    myfree(ids);
#endif
    
    t1 = my_second();
    
    if(ThisTask == 0)
    {
        printf("success.  took=%g sec\n", timediff(t0, t1));
    }
}

int compare_IDs(const void *a, const void *b)
{
    if(*((MyIDType *) a) < *((MyIDType *) b))
        return -1;
    
    if(*((MyIDType *) a) > *((MyIDType *) b))
        return +1;
    
    return 0;
}
