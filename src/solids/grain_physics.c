#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*
 
 This module contains the self-contained sub-routines needed for
 grain-specific physics in proto-planetary/proto-stellar/planetary cases.
 It's also potentially use-able for GMC and ISM scales, and terrestrial
 turbulence. Anything where aerodynamic particles are interesting
 
 
 This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 
 */



#ifdef GRAIN_FLUID


/* function to apply the drag on the grains from surrounding gas properties */
void apply_grain_dragforce(void)
{
    
    CPU_Step[CPU_MISC] += measure_time();
    
    int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((P[i].Type != 0)&&(P[i].Type != 4))
        {
#ifdef BOX_BND_PARTICLES
            if(P[i].ID > 0)
#endif
            if(P[i].Gas_Density > 0)
            {
                double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                if(dt > 0)
                {
                    double cs = sqrt( GAMMA * GAMMA_MINUS1 * P[i].Gas_InternalEnergy);
                    double R_grain_cgs = P[i].Grain_Size;
                    double R_grain_code = R_grain_cgs / (All.UnitLength_in_cm / All.HubbleParam);
                    double rho_gas = P[i].Gas_Density * All.cf_a3inv;
                    double rho_grain_physical = All.Grain_Internal_Density; // cgs units //
                    double rho_grain_code = rho_grain_physical / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam); // code units //
                    double vgas_mag = 0.0;
                    for(k=0;k<3;k++) {vgas_mag+=(P[i].Gas_Velocity[k]-P[i].Vel[k])*(P[i].Gas_Velocity[k]-P[i].Vel[k]);}
                    
                    
                    if(vgas_mag > 0)
                    {
                        vgas_mag = sqrt(vgas_mag) / All.cf_atime;
                        double x0 = 0.469993*sqrt(GAMMA) * vgas_mag/cs; // (3/8)*sqrt[pi/2]*|vgas-vgrain|/cs //
                        double tstop_inv = 1.59577/sqrt(GAMMA) * rho_gas * cs / (R_grain_code * rho_grain_code); // 2*sqrt[2/pi] * 1/tstop //



#ifdef GRAIN_LORENTZFORCE
                        /* calculate the grain charge following Draine & Sutin */
                        double cs_cgs = cs * All.UnitVelocity_in_cm_per_s;
                        double tau_draine_sutin = R_grain_cgs * (2.3*PROTONMASS) * (cs_cgs*cs_cgs) / (GAMMA * ELECTRONCHARGE*ELECTRONCHARGE);
                        double Z_grain = -DMAX( 1./(1. + sqrt(1.0e-3/tau_draine_sutin)) , 2.5*tau_draine_sutin );
                        if(isnan(Z_grain)||(Z_grain>=0)) {Z_grain=0;}
#endif

#ifdef GRAIN_EPSTEIN_STOKES
                        double mu = 2.3 * PROTONMASS;
                        double temperature = mu * (P[i].Gas_InternalEnergy*All.UnitEnergy_in_cgs*All.HubbleParam/All.UnitMass_in_g) / BOLTZMANN;
                        double cross_section = GRAIN_EPSTEIN_STOKES * 2.0e-15 * (1. + 70./temperature);
                        cross_section /= (All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam*All.HubbleParam));
                        double n_mol = rho_gas / (mu * All.HubbleParam/All.UnitMass_in_g);
                        double mean_free_path = 1 / (n_mol * cross_section); // should be in code units now //
                        double corr_mfp = R_grain_code / ((9./4.) * mean_free_path);
                        if(corr_mfp > 1) {tstop_inv /= corr_mfp;}
#ifdef GRAIN_LORENTZFORCE
                        /* also have charged grains, so we will calculate Coulomb forces as well */
                        double a_Coulomb = sqrt(2.*GAMMA*GAMMA*GAMMA/(9.*M_PI));
                        double tstop_Coulomb_inv = 0.797885/sqrt(GAMMA) * rho_gas * cs / (R_grain_code * rho_grain_code); // base normalization //
                        tstop_Coulomb_inv /= (1. + a_Coulomb *(vgas_mag/cs)*(vgas_mag/cs)*(vgas_mag/cs)) * sqrt(1.+x0*x0); // velocity dependence (force becomes weak when super-sonic)
                        tstop_Coulomb_inv *= (Z_grain/tau_draine_sutin) * (Z_grain/tau_draine_sutin) / 17.; // coulomb attraction terms, assuming ions have charge ~1, and Coulomb logarithm is 17
                        // don't need super-accuration gas ionization states, just need approximate estimate, which we can make based on temperature //
                        double T_Kelvin = (2.3*PROTONMASS) * (cs_cgs*cs_cgs) / (1.3807e-16 * GAMMA), f_ion_to_use = 0; // temperature in K
#ifdef COOLING  // in this case, have the ability to calculate more accurate ionization fraction
                        {
                            double u_tmp, ne_tmp = 1, nh0_tmp = 0, mu_tmp = 1, temp_tmp, nHeII_tmp, nhp_tmp, nHe0_tmp, nHepp_tmp;
                            u_tmp = 1.3807e-16 * T_Kelvin / (2.3*PROTONMASS) * (All.UnitMass_in_g/All.UnitEnergy_in_cgs); // needs to be in code units
                            temp_tmp = ThermalProperties(u_tmp, rho_gas, -1, &mu_tmp, &ne_tmp, &nh0_tmp, &nhp_tmp, &nHe0_tmp, &nHeII_tmp, &nHepp_tmp);
                            f_ion_to_use = DMIN(ne_tmp , 1.);
                        }
#else           // without cooling active, use a simple approximate guess for ionization
                        if(T_Kelvin > 1000.) {f_ion_to_use = exp(-15000./T_Kelvin);}
#endif
                        tstop_Coulomb_inv *= f_ion_to_use; // correct for ionization fraction
                        tstop_inv += tstop_Coulomb_inv; // add both forces
#endif // LORENTZ force (Coulomb computation)
#endif
                        double C1 = (-1-sqrt(1+x0*x0)) / x0;
                        double xf = 0.0;
                        double dt_tinv = dt * tstop_inv;
                        if(dt_tinv < 100.)
                        {
                            double C2 = C1 * exp( dt_tinv );
                            xf = -2 * C2 / (C2*C2 -1);
                        }
                        double slow_fac = 1 - xf / x0;
                        // note that, with an external (gravitational) acceleration, we can still solve this equation for the relevant update //

                        double dv[3]={0}, external_forcing[3]={0};
                        for(k=0;k<3;k++) {external_forcing[k] = 0;}
                        /* this external_forcing parameter includes additional grain-specific forces. note that -anything- which imparts an 
                            identical acceleration onto gas and dust will cancel in the terms in t_stop, and just act like a 'normal' acceleration
                            on the dust. for this reason the gravitational acceleration doesn't need to enter our 'external_forcing' parameter */
#ifdef GRAIN_LORENTZFORCE
                        /* Lorentz force on a grain = Z*e/c * ([v_grain-v_gas] x B) */
                        double v_cross_B[3];
                        v_cross_B[0] = (P[i].Vel[1]-P[i].Gas_Velocity[1])*P[i].Gas_B[2] - (P[i].Vel[2]-P[i].Gas_Velocity[2])*P[i].Gas_B[1];
                        v_cross_B[1] = (P[i].Vel[2]-P[i].Gas_Velocity[2])*P[i].Gas_B[0] - (P[i].Vel[0]-P[i].Gas_Velocity[0])*P[i].Gas_B[2];
                        v_cross_B[2] = (P[i].Vel[0]-P[i].Gas_Velocity[0])*P[i].Gas_B[1] - (P[i].Vel[1]-P[i].Gas_Velocity[1])*P[i].Gas_B[0];

                        double grain_mass = (4.*M_PI/3.) * R_grain_code*R_grain_code*R_grain_code * rho_grain_code; // code units
                        double lorentz_units = sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam); // code B to Gauss
                        lorentz_units *= (ELECTRONCHARGE/C) * All.UnitVelocity_in_cm_per_s / (All.UnitMass_in_g / All.HubbleParam); // converts acceleration to cgs
                        lorentz_units /= All.UnitVelocity_in_cm_per_s / (All.UnitTime_in_s / All.HubbleParam); // converts it to code-units acceleration

                        /* define unit vectors and B for evolving the lorentz force */
                        double bhat[3]={0}, bmag=0, efield[3]={0}, efield_coeff=0;
                        for(k=0;k<3;k++) {bhat[k]=P[i].Gas_B[k]; bmag+=bhat[k]*bhat[k]; dv[k]=P[i].Vel[k]-P[i].Gas_Velocity[k];}
                        if(bmag>0) {bmag=sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}} else {bmag=0;}
                        double grain_charge_cinv = Z_grain / grain_mass * lorentz_units;
#ifdef GRAIN_RDI_TESTPROBLEM
                        if(All.Grain_Charge_Parameter != 0) {grain_charge_cinv = -All.Grain_Charge_Parameter/All.Grain_Size_Max * pow(All.Grain_Size_Max/P[i].Grain_Size,2);} // set charge manually //
#endif
                        /* now apply the boris integrator */
                        double lorentz_coeff = (0.5*dt) * bmag * grain_charge_cinv; // dimensionless half-timestep term for boris integrator //
                        double v_m[3]={0}, v_t[3]={0}, v_p[3]={0}, vcrosst[3]={0};
                        for(k=0;k<3;k++) {v_m[k] = dv[k] + 0.5*efield_coeff*efield[k];} // half-step from E-field
                        /* cross-product for rotation */
                        vcrosst[0] = v_m[1]*bhat[2] - v_m[2]*bhat[1]; vcrosst[1] = v_m[2]*bhat[0] - v_m[0]*bhat[2]; vcrosst[2] = v_m[0]*bhat[1] - v_m[1]*bhat[0];
                        for(k=0;k<3;k++) {v_t[k] = v_m[k] + lorentz_coeff * vcrosst[k];} // first half-rotation
                        vcrosst[0] = v_t[1]*bhat[2] - v_t[2]*bhat[1]; vcrosst[1] = v_t[2]*bhat[0] - v_t[0]*bhat[2]; vcrosst[2] = v_t[0]*bhat[1] - v_t[1]*bhat[0];
                        for(k=0;k<3;k++) {v_p[k] = v_m[k] + (2.*lorentz_coeff/(1.+lorentz_coeff*lorentz_coeff)) * vcrosst[k];} // second half-rotation
                        for(k=0;k<3;k++) {v_p[k] += 0.5*efield_coeff*efield[k];} // half-step from E-field
                        /* calculate effective acceleration from discrete step in velocity */
                        for(k=0;k<3;k++) {external_forcing[k] += (v_p[k] - dv[k]) / dt;} // boris integrator
                        //for(k=0;k<3;k++) {external_forcing[k] += grain_charge_cinv * v_cross_B[k];} // standard explicit integrator



                        /* note: if grains moving super-sonically with respect to gas, and charge equilibration time is much shorter than the 
                            streaming/dynamical timescales, then the charge is slightly reduced, because the ion collision rate is increased while the 
                            electron collision rate is increased less (since electrons are moving much faster, we assume the grain is still sub-sonic 
                            relative to the electron sound speed. in this case, for the large-grain limit, the Draine & Sutin results can be generalized; 
                            the full expressions are messy but can be -approximated- fairly well for Mach numbers ~3-30 by simply 
                            suppressing the equilibrium grain charge by a power ~exp[-0.04*mach]  (weak effect, though can be significant for mach>10) */
#endif
                        
                        double delta_egy = 0;
                        double delta_mom[3];
                        for(k=0; k<3; k++)
                        {
                            /* measure the imparted energy and momentum as if there were no external acceleration */
                            double v_init = P[i].Vel[k];
                            double vel_new = v_init + slow_fac * (P[i].Gas_Velocity[k]-v_init);
                            /* now calculate the updated velocity accounting for any external, non-standard accelerations */
                            double vdrift = 0;
                            if(tstop_inv > 0) {vdrift = external_forcing[k] / (tstop_inv * sqrt(1+x0*x0));}
                            dv[k] = slow_fac * (P[i].Gas_Velocity[k] - v_init + vdrift);
                            if(isnan(vdrift)||isnan(slow_fac)) {dv[k] = 0;}

			    vel_new = v_init + dv[k];
                            delta_mom[k] = P[i].Mass * (vel_new - v_init);
                            delta_egy += 0.5*P[i].Mass * (vel_new*vel_new - v_init*v_init);

                            /* note, we can directly apply this by taking P[i].Vel[k] += dv[k]; but this is not as accurate as our
                                normal leapfrog integration scheme.
                                we can also account for the -gas- acceleration, by including it like vdrift;
                                for a constant t_stop, the gas acceleration term appears as 
                                P[i].Vel[l] += Gas_Accel[k] * dt + slow_fac * (Gas-Accel[k] / tstop_inv) */
                            /* note that we solve the equations with an external acceleration already (external_forcing above): therefore add to forces
                             like gravity that are acting on the gas and dust in the same manner (in terms of acceleration) */
                            P[i].GravAccel[k] += dv[k] / dt;
                            //P[i].Vel[k] += dv[k];
                        }

                    
#ifdef GRAIN_BACKREACTION
                        double dvel, degy, r2nearest, *pos,h,h2,hinv,hinv3,hinv4,r2,rho,u,wk,dwk;
                        int N_MAX_KERNEL,N_MIN_KERNEL,MAXITER_FB,NITER,startnode,dummy,numngb_inbox,jnearest,j,k,n;
                        Ngblist = (int *) mymalloc("Ngblist",NumPart * sizeof(int));
                        
                        /* now add in a loop to find particles in same domain, share back the
                         momentum and energy to them (to be properly conservative) */
                        N_MIN_KERNEL=4;N_MAX_KERNEL=10.*All.DesNumNgb;MAXITER_FB=30;NITER=0;jnearest=0;
                        startnode=All.MaxPart;dummy=0;h=0;numngb_inbox=0;pos=P[i].Pos;
                        h=PPP[i].Hsml; if(h<=0) h=All.SofteningTable[0];
                        do {
                            numngb_inbox = ngb_treefind_variable_targeted(pos,h,-1,&startnode,0,&dummy,&dummy,1); // search for gas: 2^0=1
                            h2=h*h; hinv=1/h; hinv3=hinv*hinv*hinv; hinv4=hinv3*hinv; rho=0;
                            if((numngb_inbox>=N_MIN_KERNEL)&&(numngb_inbox<=N_MAX_KERNEL))
                            {
                                jnearest=0;r2nearest=1.0e10;
                                for(n=0; n<numngb_inbox; n++)
                                {
                                    j = Ngblist[n];
                                    r2=0;for(k=0;k<3;k++) r2+=(P[i].Pos[k]-P[j].Pos[k])*(P[i].Pos[k]-P[j].Pos[k]);
                                    if((r2<r2nearest)&&(P[j].Mass>0)) {
                                        r2nearest=r2;jnearest=j;
                                    }
                                    if ((r2<=h2)&&(P[j].Mass>0)&&(SphP[j].Density>0)) {
                                        u=sqrt(r2)*hinv;
                                        kernel_main(u,hinv3,hinv4,&wk,&dwk,0);
                                        rho += (P[j].Mass*wk);
                                    }
                                } /* for(n=0; n<numngb_inbox; n++) */
                            } /* if(numngb_inbox>0) */
                            else
                            {
                                startnode=All.MaxPart;
                                if(numngb_inbox<N_MIN_KERNEL)
                                {
                                    if(numngb_inbox<=0) {
                                        h*=2.0;
                                    } else {
                                        if(NITER<=5)
                                            h*=pow((float)numngb_inbox/(float)N_MIN_KERNEL,-1/NUMDIMS);
                                        else
                                            h*=1.26; /* iterate until find appropriate > N_MIN # particles */
                                    }
                                }
                                if(numngb_inbox>N_MAX_KERNEL)
                                {
                                    if(NITER<=5)
                                        h*=pow((float)numngb_inbox/(float)N_MAX_KERNEL,-1/NUMDIMS);
                                    else
                                        h/=1.31; /* iterate until find appropriate < N_MAX # particles */
                                }
                            }
                            NITER++;
                        } while((startnode >= 0)&&(NITER<=MAXITER_FB));
                        if(jnearest != 0) {if((P[jnearest].Mass<=0)||(SphP[jnearest].Density<=0)) jnearest=0;}
                        
                        if((jnearest != 0)&&(numngb_inbox>0))
                        {
                            // share the coupled momentum and energy back to the nearby gas //
                            if(rho>0)
                            {
                                for(n=0; n<numngb_inbox; n++)
                                {
                                    j = Ngblist[n];
#ifdef BOX_BND_PARTICLES
                                    if(P[j].ID <= 0) continue;
#endif
                                    r2=0; for(k=0;k<3;k++) r2+=(P[i].Pos[k]-P[j].Pos[k])*(P[i].Pos[k]-P[j].Pos[k]);
                                    if ((r2<=h2)&&(P[j].Mass>0)&&(SphP[j].Density>0))
                                    {
                                        u=sqrt(r2)*hinv;
                                        kernel_main(u,hinv3,hinv4,&wk,&dwk,0);
                                        wk *= P[j].Mass / rho;
                                        degy=0;
                                        MyDouble VelPred_j[3];
                                        for(k=0;k<3;k++) {VelPred_j[k]=P[j].Vel[k];}
#ifdef BOX_SHEARING
                                        if(P[i].Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                                        if(P[i].Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                                        for(k=0; k<3; k++)
                                        {
                                            dvel=-wk*delta_mom[k]/P[j].Mass;
                                            degy-=0.5*P[j].Mass*((VelPred_j[k]+dvel)*(VelPred_j[k]+dvel) - VelPred_j[k]*VelPred_j[k]);
                                            P[j].Vel[k] += dvel;
                                            SphP[j].VelPred[k] += dvel;
                                        }
                                        degy += -wk*delta_egy; // energy 'donated' - kinetic energy change = thermal energy change
                                        // degy is now the change in thermal energy, convert to specific energy change //
                                        degy *= 1 / P[j].Mass;
                                        if(degy<-0.9*SphP[j].InternalEnergy) SphP[j].InternalEnergy*=0.1; else SphP[j].InternalEnergy += degy;
                                    }
                                } // for(n=0; n<numngb_inbox; n++)
                            } else {
                                j=jnearest;
#ifdef BOX_BND_PARTICLES
                                if(P[j].ID <= 0) continue;
#endif
                                wk=1; degy=0;
                                MyDouble VelPred_j[3];
                                for(k=0;k<3;k++) {VelPred_j[k]=P[j].Vel[k];}
#ifdef BOX_SHEARING
                                if(P[i].Pos[0] - P[j].Pos[0] > +boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;}
                                if(P[i].Pos[0] - P[j].Pos[0] < -boxHalf_X) {VelPred_j[BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;}
#endif
                                for(k=0; k<3; k++)
                                {
                                    dvel=-wk*delta_mom[k]/P[j].Mass;
                                    degy-=0.5*P[j].Mass*((VelPred_j[k]+dvel)*(VelPred_j[k]+dvel) - VelPred_j[k]*VelPred_j[k]);
                                    P[j].Vel[k] += dvel;
                                    SphP[j].VelPred[k] += dvel;
                                }
                                degy += -wk*delta_egy; // energy 'donated' - kinetic energy change = thermal energy change
                                // degy is now the change in thermal energy, convert to specific energy change //
                                degy *= 1 / P[j].Mass;
                                if(degy<-0.9*SphP[j].InternalEnergy) SphP[j].InternalEnergy*=0.1; else SphP[j].InternalEnergy += degy;
                            } // if(rho>0) else
                        } // closes if((jnearest != 0)&&(rho>0)) //
                        
                        myfree(Ngblist);
#endif // closes GRAIN_BACKREACTION

                    } // closes check for if(v_mag > 0)
                } // closes check for if(dt > 0)
            } // closes check for if(P[i].Gas_Density > 0)
        } // closes check for if(P[i].Type != 0)
    } // closes main particle loop
    
    CPU_Step[CPU_DRAGFORCE] += measure_time();
    
}











#endif


