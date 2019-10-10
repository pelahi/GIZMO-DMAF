#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"

/* Routines for the drift/predict step */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially in detail (although the actual algorithm 
 * structure remains essentially the same) 
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void reconstruct_timebins(void)
{
    int i, bin;
    long long glob_sum;
    
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        TimeBinCount[bin] = 0;
        TimeBinCountSph[bin] = 0;
        FirstInTimeBin[bin] = -1;
        LastInTimeBin[bin] = -1;
#ifdef GALSF
        TimeBinSfr[bin] = 0;
#endif
#ifdef BLACK_HOLES
        TimeBin_BH_mass[bin] = 0;
        TimeBin_BH_dynamicalmass[bin] = 0;
        TimeBin_BH_Mdot[bin] = 0;
        TimeBin_BH_Medd[bin] = 0;
#endif
    }
    
    for(i = 0; i < NumPart; i++)
    {
        bin = P[i].TimeBin;
        
        if(TimeBinCount[bin] > 0)
        {
            PrevInTimeBin[i] = LastInTimeBin[bin];
            NextInTimeBin[i] = -1;
            NextInTimeBin[LastInTimeBin[bin]] = i;
            LastInTimeBin[bin] = i;
        }
        else
        {
            FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
            PrevInTimeBin[i] = NextInTimeBin[i] = -1;
        }
        TimeBinCount[bin]++;
        if(P[i].Type == 0)
            TimeBinCountSph[bin]++;
        
#ifdef GALSF
        if(P[i].Type == 0)
            TimeBinSfr[bin] += SphP[i].Sfr;
#endif
#ifdef BLACK_HOLES
        if(P[i].Type == 5)
        {
            TimeBin_BH_mass[bin] += BPP(i).BH_Mass;
            TimeBin_BH_dynamicalmass[bin] += P[i].Mass;
            TimeBin_BH_Mdot[bin] += BPP(i).BH_Mdot;
            TimeBin_BH_Medd[bin] += BPP(i).BH_Mdot / BPP(i).BH_Mass;
        }
#endif
    }
    
    make_list_of_active_particles();
    
    for(i = FirstActiveParticle, NumForceUpdate = 0; i >= 0; i = NextActiveParticle[i])
    {
        NumForceUpdate++;
        if(i >= NumPart)
        {
            printf("Bummer i=%d\n", i);
            terminate("inconsistent list");
        }
    }
    
    sumup_large_ints(1, &NumForceUpdate, &glob_sum);
    
    GlobNumForceUpdate = glob_sum;
}





void drift_particle(int i, integertime time1)
{
    int j;
    double dt_drift;
    
    integertime time0 = P[i].Ti_current;
    
    if(time1 < time0)
    {
        printf("i=%d time0=%d time1=%d\n", i, (int)time0, (int)time1);
        terminate("no prediction into past allowed");
    }
    
    if(time1 == time0)
        return;
    
    if(All.ComovingIntegrationOn)
        dt_drift = get_drift_factor(time0, time1);
    else
        dt_drift = (time1 - time0) * All.Timebase_interval;
    
    
#if !defined(FREEZE_HYDRO)
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(P[i].Type==0) {advect_mesh_point(i,dt_drift);} else {for(j=0;j<3;j++) {P[i].Pos[j] += P[i].Vel[j] * dt_drift;}}
#else
    for(j=0;j<3;j++) {P[i].Pos[j] += P[i].Vel[j] * dt_drift;}
#endif
#endif
#if (NUMDIMS==1)
    P[i].Pos[1]=P[i].Pos[2]=0; // force zero-ing
#endif
#if (NUMDIMS==2)
    P[i].Pos[2]=0; // force zero-ing
#endif
    
    double divv_fac = P[i].Particle_DivVel * dt_drift;
    double divv_fac_max = 0.3; //1.5; // don't allow Hsml to change too much in predict-step //
    if(divv_fac > +divv_fac_max) divv_fac = +divv_fac_max;
    if(divv_fac < -divv_fac_max) divv_fac = -divv_fac_max;

// DM ANNIHILATION at GAS
#if defined(DMANNIHILATION)
    double divv_facDM = P[i].Particle_DivVelDM * dt_drift;
    double divv_fac_maxDM = 0.3; //1.5; // don't allow Hsml to change too much in predict-step //
    if(divv_facDM > +divv_fac_maxDM) divv_facDM = +divv_fac_maxDM;
    if(divv_facDM < -divv_fac_maxDM) divv_facDM = -divv_fac_maxDM;
#endif

#ifdef GRAIN_FLUID
    if(P[i].Type > 0)
    {
        PPP[i].Hsml *= exp((double)divv_fac / ((double)NUMDIMS));
        if(PPP[i].Hsml < All.MinHsml) {PPP[i].Hsml = All.MinHsml;}
        if(PPP[i].Hsml > All.MaxHsml) {PPP[i].Hsml = All.MaxHsml;}
    }
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(P[i].Type>0)
    {
        if(dt_drift>0)
        {
            double minsoft = ags_return_minsoft(i);
            double maxsoft = ags_return_maxsoft(i);
            PPP[i].AGS_Hsml *= exp((double)divv_fac / ((double)NUMDIMS));
            if(PPP[i].AGS_Hsml < minsoft) {PPP[i].AGS_Hsml = minsoft;}
            if(PPP[i].AGS_Hsml > maxsoft) {PPP[i].AGS_Hsml = maxsoft;}
        }
    }
#endif
    
    
    if((P[i].Type == 0) && (P[i].Mass > 0))
        {
            double dt_gravkick, dt_hydrokick, dt_entr;
            
            if(All.ComovingIntegrationOn)
            {
                dt_entr = dt_hydrokick = (time1 - time0) * All.Timebase_interval / All.cf_hubble_a;
                dt_gravkick = get_gravkick_factor(time0, time1);
            }
            else
            {
                dt_entr = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
            }
            
#ifdef PMGRID
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick +
                    SphP[i].HydroAccel[j]*All.cf_atime * dt_hydrokick; /* make sure v is in code units */
#else
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick +
                    SphP[i].HydroAccel[j]*All.cf_atime * dt_hydrokick; /* make sure v is in code units */
#endif
            
#if defined(TURB_DRIVING)
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += SphP[i].TurbAccel[j] * dt_gravkick;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            for(j = 0; j < 3; j++)
                SphP[i].VelPred[j] += SphP[i].RadAccel[j] * All.cf_atime * dt_hydrokick;
#endif
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            P[i].Mass = DMAX(P[i].Mass + SphP[i].DtMass * dt_entr, 0.5 * SphP[i].MassTrue);
#endif
            
            SphP[i].Density *= exp(-divv_fac);

// DM ANNIHILATION at GAS  
#if defined(DMANNIHILATION)
	    	P[i].DensityDM *= exp(-divv_facDM);
#if (HYDRO_FIX_MESH_MOTION > 0)
	    	P[i].HsmlDM *= exp((double)divv_facDM / ((double)NUMDIMS));
            if(P[i].HsmlDM < All.MinHsml) {P[i].HsmlDM = All.MinHsml;}
            if(P[i].HsmlDM > All.MaxHsml) {P[i].HsmlDM = All.MaxHsml;}
#endif
#endif

            double etmp = SphP[i].InternalEnergyPred + SphP[i].DtInternalEnergy * dt_entr;

            if(etmp<0.5*SphP[i].InternalEnergyPred) {SphP[i].InternalEnergyPred *= 0.5;} else {SphP[i].InternalEnergyPred=etmp;}
            if(SphP[i].InternalEnergyPred<All.MinEgySpec) SphP[i].InternalEnergyPred=All.MinEgySpec;
            
#ifdef HYDRO_PRESSURE_SPH
            SphP[i].EgyWtDensity *= exp(-divv_fac);
#endif
            
#if (HYDRO_FIX_MESH_MOTION > 0)
            PPP[i].Hsml *= exp((double)divv_fac / ((double)NUMDIMS));
            if(PPP[i].Hsml < All.MinHsml) {PPP[i].Hsml = All.MinHsml;}
            if(PPP[i].Hsml > All.MaxHsml) {PPP[i].Hsml = All.MaxHsml;}
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            PPP[i].AGS_Hsml = PPP[i].Hsml;
#endif
#endif
            drift_sph_extra_physics(i, time0, time1, dt_entr);

        
            SphP[i].Pressure = get_pressure(i);
#ifdef EOS_ENFORCE_ADIABAT
            SphP[i].InternalEnergyPred = SphP[i].Pressure / (SphP[i].Density * GAMMA_MINUS1);
#endif
        }
    
    /* check for reflecting boundaries: if so, do the reflection! */
#if defined(BOX_REFLECT_X) || defined(BOX_REFLECT_Y) || defined(BOX_REFLECT_Z)
    double box_upper[3]; box_upper[0]=box_upper[1]=box_upper[2]=1;
#ifdef BOX_PERIODIC
    box_upper[0]=boxSize_X; box_upper[1]=boxSize_Y; box_upper[2]=boxSize_Z;
#endif
    for(j = 0; j < 3; j++)
    {
        /* skip the non-reflecting boundaries */
#ifndef BOX_REFLECT_X
        if(j==0) continue;
#endif
#ifndef BOX_REFLECT_Y
        if(j==1) continue;
#endif
#ifndef BOX_REFLECT_Z
        if(j==2) continue;
#endif
        if(P[i].Pos[j] <= 0)
        {
            if(P[i].Vel[j]<0) {P[i].Vel[j]=-P[i].Vel[j]; if(P[i].Type==0) {SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0;}}
            P[i].Pos[j]=(0+((double)P[i].ID)*1.e-9)*box_upper[j];
        }
        if(P[i].Pos[j] >= box_upper[j])
        {
            if(P[i].Vel[j]>0) {P[i].Vel[j]=-P[i].Vel[j]; if(P[i].Type==0) {SphP[i].VelPred[j]=P[i].Vel[j]; SphP[i].HydroAccel[j]=0;}}
            P[i].Pos[j]=box_upper[j]*(1-((double)P[i].ID)*1.e-9);
        }
    }
#endif

    
    P[i].Ti_current = time1;
}





void move_particles(integertime time1)
{
    int i;
    for(i = 0; i < NumPart; i++)
        drift_particle(i, time1);
}





void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
#ifdef MAGNETIC
    int kB;
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    for(kB=0;kB<3;kB++) {SphP[i].BPred[kB] += SphP[i].DtB[kB] * dt_entr * BphysVolphys_to_BcodeVolCode;} // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi fluxes (otherwise coefficient is 1)
    double dtphi_code = (PhiphysVolphys_to_PhicodeVolCode) * SphP[i].DtPhi;
    SphP[i].PhiPred += dtphi_code  * dt_entr;
    double t_damp = Get_Particle_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)))
    {
        SphP[i].PhiPred *= exp( -dt_entr * t_damp );
    }
#endif
#ifdef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    for(k=0;k<3;k++) {SphP[i].B[k]=SphP[i].BPred[k];}
#ifdef DIVBCLEANING_DEDNER
    SphP[i].Phi=SphP[i].PhiPred;
#endif
#endif
#endif
#ifdef RADTRANSFER
    rt_update_driftkick(i,dt_entr,1);
#endif
#ifdef EOS_ELASTIC
    elastic_body_update_driftkick(i,dt_entr,1);
#endif
}





/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef BOX_PERIODIC
void do_box_wrapping(void)
{
    int i, j;
    double boxsize[3];
    boxsize[0] = boxSize_X;
    boxsize[1] = boxSize_Y;
    boxsize[2] = boxSize_Z;
    
    for(i = 0; i < NumPart; i++)
    {
        for(j = 0; j < 3; j++)
        {
            while(P[i].Pos[j] < 0)
            {
                P[i].Pos[j] += boxsize[j];
#ifdef BOX_SHEARING
                if(j==0)
                {
                    P[i].Vel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
                    if(P[i].Type==0)
                    {
                        SphP[i].VelPred[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                        SphP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#endif
                    }
#if (BOX_SHEARING > 1)
                    /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                    P[i].Pos[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Pos_Offset;
#endif
                }
#endif
            }
            
            while(P[i].Pos[j] >= boxsize[j])
            {
                P[i].Pos[j] -= boxsize[j];
#ifdef BOX_SHEARING
                if(j==0)
                {
                    P[i].Vel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
                    if(P[i].Type==0)
                    {
                        SphP[i].VelPred[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                        SphP[i].ParticleVel[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#endif
                    }
#if (BOX_SHEARING > 1)
                    /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                    P[i].Pos[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Pos_Offset;
#endif
                }
#endif
            }
        }
    }
}
#endif




/* ====================================================================== */
/* ================== Functions for physical information ================ */
/* ====================================================================== */


/* this function returns the effective (grid-equivalent) particle 'size'; useful for things like 
    time-stepping and limiter functions */
double INLINE_FUNC Get_Particle_Size(int i)
{
    /* in previous versions of the code, we took NumNgb^(1/NDIMS) here; however, now we 
        take that when NumNgb is computed (at the end of the density routine), so we 
        don't have to re-compute it each time. That makes this function fast enough to 
        call -inside- of loops (e.g. hydro computations) */
#if (NUMDIMS == 1)
    return 2.00000 * PPP[i].Hsml / PPP[i].NumNgb;
#endif
#if (NUMDIMS == 2)
    return 1.25331 * PPP[i].Hsml / PPP[i].NumNgb; // sqrt(Pi/2)
#endif
#if (NUMDIMS == 3)
    return 1.61199 * PPP[i].Hsml / PPP[i].NumNgb; // (4pi/3)^(1/3)
#endif
}



double INLINE_FUNC Get_Particle_Expected_Area(double h)
{
#if (NUMDIMS == 1)
    return 2;
#endif
#if (NUMDIMS == 2)
    return 2 * M_PI * h;
#endif
#if (NUMDIMS == 3)
    return 4 * M_PI * h * h;
#endif
}


/* return the estimated local column (physical units) from integrating the gradient in the density (separated here for convenience) */
double evaluate_NH_from_GradRho(MyFloat gradrho[3], double hsml, double rho, double numngb_ndim, double include_h)
{
    double gradrho_mag;
    if(rho<=0)
    {
        gradrho_mag = 0;
    } else {
        gradrho_mag = sqrt(gradrho[0]*gradrho[0]+gradrho[1]*gradrho[1]+gradrho[2]*gradrho[2]);
        if(gradrho_mag > 0) {gradrho_mag = rho*rho/gradrho_mag;} else {gradrho_mag=0;}
        if(include_h > 0) if(numngb_ndim > 0) gradrho_mag += include_h * rho * hsml / numngb_ndim; // quick-and-dirty approximation to the effective neighbor number needed here
        //if(include_h > 0) gradrho_mag += include_h * rho * (hsml * (0.124 + 11.45 / (26.55 + All.DesNumNgb))); // quick-and-dirty approximation to the effective neighbor number needed here
        // account for the fact that 'h' is much larger than the inter-particle separation //
    }
    return gradrho_mag * All.cf_a2inv; // (physical units) // *(Z/Zsolar) add metallicity dependence
}







#ifdef MAGNETIC
double INLINE_FUNC Get_Particle_BField(int i_particle_id, int k_vector_component)
{
    return SphP[i_particle_id].BPred[k_vector_component] * SphP[i_particle_id].Density / P[i_particle_id].Mass;
}

/* this function is needed to control volume fluxes of the normal components of B and phi in the 
    -bad- situation where the meshless method 'faces' do not properly close (usually means you are 
    using boundary conditions that you should not) */
double Get_DtB_FaceArea_Limiter(int i)
{
#ifdef HYDRO_SPH
    return 1;
#else
    /* define some variables */
    double dt_entr = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
    int j; double dB[3],dBmag=0,Bmag=0;
    /* check the magnitude of the predicted change in B-fields, vs. B-magnitude */
    for(j=0;j<3;j++)
    {
        dB[j] = SphP[i].DtB[j] * dt_entr / All.cf_atime; /* converts to code units of Vol_code*B_code = Vol_phys*B_phys/a */
        dBmag += dB[j]*dB[j];
        Bmag += SphP[i].BPred[j]*SphP[i].BPred[j];
    }
    dBmag = sqrt(dBmag); Bmag = sqrt(Bmag);
    /* also make sure to check the actual pressure, since if P>>B, we will need to allow larger changes in B per timestep */
    double P_BV_units = sqrt(2.*SphP[i].Pressure*All.cf_a3inv)*P[i].Mass/SphP[i].Density * All.cf_afac3 / All.cf_a2inv;
    /* the above should be in CODE Bcode*Vol_code units! */
    double Bmag_max = DMAX(Bmag, DMIN( P_BV_units, 10.*Bmag ));
    /* now check how accurately the cell is 'closed': the face areas are ideally zero */
    double area_sum = fabs(SphP[i].Face_Area[0])+fabs(SphP[i].Face_Area[1])+fabs(SphP[i].Face_Area[2]);
    /* but this needs to be normalized to the 'expected' area given Hsml */
    double area_norm = Get_Particle_Expected_Area(PPP[i].Hsml * All.cf_atime);
    /* ok, with that in hand, define an error tolerance based on this */
    if(area_norm>0)
    {
        double area_norm_min_threshold = 0.001;
        double area_norm_weight = 200.0;
#ifdef PM_HIRES_REGION_CLIPPING
        area_norm_min_threshold *= 0.01;
        area_norm_weight *= 2.5; // can be as low as 1.0 (PFH) //
#endif
        if(area_sum/area_norm > area_norm_min_threshold)
        {
            double tol = (All.CourantFac/0.2) * DMAX( 0.01, area_norm/(area_norm_weight * area_sum) );
            tol *= Bmag_max; /* give the limiter dimensions */
            if(dBmag > tol) {return tol/dBmag;} /* now actually check if we exceed this */
        }
    }
    return 1;
#endif
}


#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Particle_PhiField(int i_particle_id)
{
    //return SphP[i_particle_id].PhiPred * SphP[i_particle_id].Density / P[i_particle_id].Mass; // volumetric phy-flux (requires extra term compared to mass-based flux)
    return SphP[i_particle_id].PhiPred / P[i_particle_id].Mass; // mass-based phi-flux
}

double INLINE_FUNC Get_Particle_PhiField_DampingTimeInv(int i_particle_id)
{
    /* this timescale should always be returned as a -physical- time */
#ifdef HYDRO_SPH
    /* PFH: add simple damping (-phi/tau) term */
    double damping_tinv = 0.5 * All.DivBcleanParabolicSigma * (SphP[i_particle_id].MaxSignalVel*All.cf_afac3 / (All.cf_atime*Get_Particle_Size(i_particle_id)));
#else
    double damping_tinv;
#ifdef SELFGRAVITY_OFF
    damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveSpeed / Get_Particle_Size(i_particle_id); // fastest wavespeed has units of [vphys]
    //double damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveDecay * All.cf_a2inv; // no improvement over fastestwavespeed; decay has units [vphys/rphys]
#else
    // only see a small performance drop from fastestwavespeed above to maxsignalvel below, despite the fact that below is purely local (so allows more flexible adapting to high dynamic range)
    damping_tinv = 0.0;
    
    if(PPP[i_particle_id].Hsml > 0)
    {
        double h_eff = Get_Particle_Size(i_particle_id);
        double vsig2 = 0.5 * All.cf_afac3 * fabs(SphP[i_particle_id].MaxSignalVel);
        double phi_B_eff = 0.0;
        if(vsig2 > 0) {phi_B_eff = Get_Particle_PhiField(i_particle_id) / (All.cf_atime * vsig2);}
        double vsig1 = 0.0;
        if(SphP[i_particle_id].Density > 0)
        {
            vsig1 = All.cf_afac3 *
            sqrt( Particle_effective_soundspeed_i(i_particle_id)*Particle_effective_soundspeed_i(i_particle_id) +
                 (All.cf_afac1 / All.cf_atime) *
                 (Get_Particle_BField(i_particle_id,0)*Get_Particle_BField(i_particle_id,0) +
                  Get_Particle_BField(i_particle_id,1)*Get_Particle_BField(i_particle_id,1) +
                  Get_Particle_BField(i_particle_id,2)*Get_Particle_BField(i_particle_id,2) +
                  phi_B_eff*phi_B_eff) / SphP[i_particle_id].Density );
        }
        vsig1 = DMAX(vsig1, vsig2);
        vsig2 = 0.0;
        int j,k;
        for(j=0;j<3;j++) for(k=0;k<3;k++) {vsig2 += SphP[i_particle_id].Gradients.Velocity[j][k]*SphP[i_particle_id].Gradients.Velocity[j][k];}
        vsig2 = sqrt(vsig2);
        vsig2 = 3.0 * h_eff * DMAX( vsig2, fabs(P[i_particle_id].Particle_DivVel)) / All.cf_atime;
        double prefac_fastest = 0.1;
        double prefac_tinv = 0.5;
        double area_0 = 0.1;
#ifdef MHD_CONSTRAINED_GRADIENT
        prefac_fastest = 1.0;
        prefac_tinv = 2.0;
        area_0 = 0.05;
        vsig2 *= 5.0;
        if(SphP[i_particle_id].FlagForConstrainedGradients <= 0) prefac_tinv *= 30;
#endif
        prefac_tinv *= sqrt(1. + SphP[i_particle_id].ConditionNumber/100.);
        double area = fabs(SphP[i_particle_id].Face_Area[0]) + fabs(SphP[i_particle_id].Face_Area[1]) + fabs(SphP[i_particle_id].Face_Area[2]);
        area /= Get_Particle_Expected_Area(PPP[i_particle_id].Hsml);
        prefac_tinv *= (1. + area/area_0)*(1. + area/area_0);
        
        double vsig_max = DMAX( DMAX(vsig1,vsig2) , prefac_fastest * All.FastestWaveSpeed );
        damping_tinv = prefac_tinv * All.DivBcleanParabolicSigma * (vsig_max / (All.cf_atime * h_eff));
    }
#endif
#endif
    return damping_tinv;
}

#endif // dedner
#endif // magnetic





/* -------------------------------------------------------------------------------------------------------------------------------------
 ------------------- the following routines are not setting the velocity, but instead are useful routines for computation of
 -------------------  various quantities needed in the mesh motion for different coordinate systems or assumed mesh shapes
 ------------------------------------------------------------------------------------------------------------------------------------- */

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
/* time-step the positions of the mesh points. this is trivial except if we are evolving the mesh points in non-cartesian coordinates
    (cylindrical or spherical) based on assumed fixed initial velocities (if HYDRO_FIX_MESH_MOTION=2 or 3),
    in which case we have to convert back and forth. */
void advect_mesh_point(int i, double dt)
{
    int k;

#if (HYDRO_FIX_MESH_MOTION == 2) || (HYDRO_FIX_MESH_MOTION == 3) // cylindrical or spherical coordinates
    // define the location relative to the origin (needed in these coordinate systems)
    double dp[3], dp_offset[3]={0}; for(k=0;k<3;k++) {dp[k]=P[i].Pos[k];} // assume center is at coordinate origin
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE) // unless we use a BH anchor, to define the center
    for(k=0;k<3;k++) {dp_offset[k] = -P[i].min_xyz_to_bh[k] + P[i].Pos[k];}
#elif defined(BOX_PERIODIC) // or if periodic, the box mid-point is instead the center
#if (NUMDIMS==1)
    dp_offset[0] = -boxHalf_X;
#elif (NUMDIMS==2)
    dp_offset[0] = -boxHalf_X; dp_offset[1] = -boxHalf_Y;
#else
    dp_offset[0] = -boxHalf_X; dp_offset[1] = -boxHalf_Y; dp_offset[2] = -boxHalf_Z;
#endif
#endif
    for(k=0;k<3;k++) {dp[k] += dp_offset[k];}
#if (HYDRO_FIX_MESH_MOTION == 2) // cylindrical
    double r2=dp[0]*dp[0]+dp[1]*dp[1], r=sqrt(r2), c0=dp[0]/r, s0=dp[1]/r, z=dp[2]; // get r, sin/cos theta, z
    double vr=c0*SphP[i].ParticleVel[0] + s0*SphP[i].ParticleVel[1], vt=s0*SphP[i].ParticleVel[0] - c0*SphP[i].ParticleVel[1]; vz=SphP[i].ParticleVel[2]; // velocities in these directions
    double r_n=r+vr*dt, z_n=z+vz*dt, c_n=c0-s0*(vt/r)*dt, s_n=s0+c0*(vt/r)*dt; // updated cylindrical values
    dp[0] = c_n*r_n; dp[1] = s_n*r_n; dp[2] = z_n; // back to coordinates
    SphP[i].ParticleVel[0] = c_n*vr + s_n*vt; // re-set velocities in these coordinates //
    SphP[i].ParticleVel[1] = s_n*vr - c_n*vt
    SphP[i].ParticleVel[2] = vz;
    return;
#elif (HYDRO_FIX_MESH_MOTION == 3) // spherical
    double dp[3],v[3],r2=0; for(k=0;k<3;k++) {r2+=dp[k]*dp[k]; v[k]=SphP[i].ParticleVel[k];} // assume center is at coordinate origin
    double r=sqrt(r2), rxy=sqrt(dp[0]*dp[0]+dp[1]*dp[1]); vr=(dp[0]*v[0] + dp[1]*v[1] + dp[2]*v[2])/r; // updated r is easy
    double ct = 1./sqrt(1.+dp[1]*dp[1]/(dp[0]*dp[0])), st = (dp[1]/dp[0])*ct; // cos and sin theta
    double cp = sqrt(1.-dp[2]*dp[2]/(r*r)), sp = dp[2]/r; // cos and sin phi
    double t_dot = (v[0]*dp[1]-v[1]*dp[0])/(rxy*rxy), p_dot = (dp[2]*(dp[0]*v[0]+dp[1]*v[1])-rxy*rxy*v[2])/(r*r*rxy); // theta, phi derivatives
    double r_n=r+vr*dt, ct_n=ct-st*t_dot, st_n=st+ct*t_dot, cp_n=cp-sp*t_dot, sp_n=sp+cp*t_dot; // updated angles and positions in spherical
    dp[0] = r_n * ct_n * cp_n; dp[1] = r_n * st_n * cp_n; dp[2] = r_n * sp_n; // back to coordinates
    rxy = sqrt(dp[0]*dp[0] + dp[1]*dp[1]); // updated rxy
    SphP[i].ParticleVel[0] = (dp[0]/r_n) * vr + dp[1] * t_dot + dp[0]*dp[2]/rxy * p_dot; // back to cartesian velocities
    SphP[i].ParticleVel[1] = (dp[1]/r_n) * vr - dp[0] * t_dot + dp[1]*dp[2]/rxy * p_dot; // back to cartesian velocities
    SphP[i].ParticleVel[2] = (dp[2]/r_n) * vr - rxy * p_dot; // back to cartesian velocities
    return;
#endif
    // ok now have the updated x/y/z positions relative to the origin, convert these back to the simulation coordinate frame
    for(k=0;k<3;k++) {P[i].Pos[k] = dp[k] - dp_offset[k];}
#endif // ok done with cylindrical/spherical coordinates
    
    
    // ok anything else ('normal' coordinates), does down here
    for(k=0;k<3;k++) {P[i].Pos[k] += SphP[i].ParticleVel[k] * dt;} // for standard grid velocities, this is trivial //
    return;
}




/* routine to calculate the overlapping face area of two cuboids in NDIMS dimensions based on their relative positions */
double calculate_face_area_for_cartesian_mesh(double *dp, double rinv, double l_side, double *Face_Area_Vec)
{
    Face_Area_Vec[0]=Face_Area_Vec[1]=Face_Area_Vec[2]=0; double Face_Area_Norm;
#if (NUMDIMS==1)
    Face_Area_Norm = 1; Face_Area_Vec[0] = Face_Area_Norm * dp[0]/fabs(dp[0]);
#elif (NUMDIMS==2)
    if(fabs(dp[0]) > fabs(dp[1])) {Face_Area_Vec[0] = Face_Area_Norm = DMAX(0,l_side-fabs(dp[1])) * dp[0]/fabs(dp[0]);} else {Face_Area_Vec[1] = Face_Area_Norm = DMAX(0,l_side-fabs(dp[0])) * dp[1]/fabs(dp[1]);}
#else
    double dp_abs[3]; int k,kdir; for(k=0;k<3;k++) {dp_abs[k] = fabs(dp[k]);}
    if((dp_abs[0]>=dp_abs[1])&&(dp_abs[0]>=dp_abs[2])) {kdir=0;} else if ((dp_abs[1]>=dp_abs[0])&&(dp_abs[1]>=dp_abs[2])) {kdir=1;} else {kdir=2;}
    Face_Area_Norm=1; for(k=0;k<3;k++) {if(k!=kdir) {Face_Area_Norm *= DMAX(0,l_side-dp_abs[k]);}}
    Face_Area_Vec[kdir] = Face_Area_Norm * dp[kdir]/fabs(dp[kdir]);
#endif
    return fabs(Face_Area_Norm);
}

#endif
