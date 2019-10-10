#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"

#ifdef COOL_GRACKLE
#include <grackle.h>
#define ENDRUNVAL 91234

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return pressure
//     4 == calculate and return gamma (only valid when COOL_GRACKLE_CHEMISTRY>0)
//
double CallGrackle(double u_old, double rho, double dt, double ne_guess, int target, int mode)
{
    gr_float returnval = 0.0;
    
    // Set grid dimension and size.
    // grid_start and grid_end are used to ignore ghost zones.
    int field_size = 1;
    int grid_rank = 3;
    int grid_dimension[3], grid_start[3], grid_end[3];
    int i;
    for (i = 0;i < 3;i++) {
        grid_dimension[i] = 1; // the active dimension not including ghost zones.
        grid_start[i]     = 0;
        grid_end[i]       = 0;
    }
    grid_dimension[0] = field_size;
    grid_end[0]       = field_size - 1;
    
    gr_float density, metal_density, energy, velx, vely, velz;
    gr_float cooling_time, temperature, pressure, gamma;
    velx          = SphP[target].VelPred[0];
    vely          = SphP[target].VelPred[1];
    velz          = SphP[target].VelPred[2];
    density       = rho;
    energy        = u_old;
#ifdef METALS
    metal_density = density * P[target].Metallicity[0];
#else
    metal_density = density * 0.02;
#endif
    gamma         = GAMMA;
    
#if (COOL_GRACKLE_CHEMISTRY >  0) // non-tabular
    gr_float ne_density;
    gr_float HI_density, HII_density, HM_density;
    gr_float HeI_density, HeII_density, HeIII_density;
    gr_float H2I_density, H2II_density;
    gr_float DI_density, DII_density, HDI_density;
    gr_float tiny = 1.0e-20;
    
    // Atomic
    ne_density    = density * ne_guess;
    
    HI_density    = density * SphP[target].grHI;  //initialized with HYDROGEN_MASSFRAC
    HII_density   = density * SphP[target].grHII;
    HM_density    = density * SphP[target].grHM;
    
    HeI_density   = density * SphP[target].grHeI;
    HeII_density  = density * SphP[target].grHeII;
    HeIII_density = density * SphP[target].grHeIII;
    
    H2I_density  = density * tiny;
    H2II_density = density * tiny;
    DI_density   = density * tiny;
    DII_density  = density * tiny;
    HDI_density  = density * tiny;
    
#if (COOL_GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
    H2I_density  = density * SphP[target].grH2I;
    H2II_density = density * SphP[target].grH2II;
#endif
    
#if (COOL_GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
    DI_density   = density * SphP[target].grDI;
    DII_density  = density * SphP[target].grDII;
    HDI_density  = density * SphP[target].grHDI;
#endif
    
    switch(mode) {
        case 0:  //solve chemistry & update values
            if(solve_chemistry(&All.GrackleUnits,
                               All.cf_atime, dt,
                               grid_rank, grid_dimension,
                               grid_start, grid_end,
                               &density, &energy,
                               &velx, &vely, &velz,
                               &HI_density, &HII_density, &HM_density,
                               &HeI_density, &HeII_density, &HeIII_density,
                               &H2I_density, &H2II_density,
                               &DI_density, &DII_density, &HDI_density,
                               &ne_density, &metal_density) == 0) {
                fprintf(stderr, "Error in solve_chemistry.\n");
                endrun(ENDRUNVAL);
            }
            
            // Assign variables back
            SphP[target].grHI    = HI_density    / density;
            SphP[target].grHII   = HII_density   / density;
            SphP[target].grHM    = HM_density    / density;
            
            SphP[target].grHeI   = HeI_density   / density;
            SphP[target].grHeII  = HeII_density  / density;
            SphP[target].grHeIII = HeIII_density / density;
            
#if (COOL_GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
            SphP[target].grH2I   = H2I_density   / density;
            SphP[target].grH2II  = H2II_density  / density;
#endif
            
#if (COOL_GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
            SphP[target].grDI    = DI_density    / density;
            SphP[target].grDII   = DII_density   / density;
            SphP[target].grHDI   = HDI_density   / density;
#endif
            returnval = energy;
            break;
            
        case 1:  //cooling time
            if(calculate_cooling_time(&All.GrackleUnits, All.cf_atime,
                                      grid_rank, grid_dimension,
                                      grid_start, grid_end,
                                      &density, &energy,
                                      &velx, &vely, &velz,
                                      &HI_density, &HII_density, &HM_density,
                                      &HeI_density, &HeII_density, &HeIII_density,
                                      &H2I_density, &H2II_density,
                                      &DI_density, &DII_density, &HDI_density,
                                      &ne_density, &metal_density,
                                      &cooling_time) == 0) {
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature
            if(calculate_temperature(&All.GrackleUnits, All.cf_atime,
                                     grid_rank, grid_dimension,
                                     grid_start, grid_end,
                                     &density, &energy,
                                     &HI_density, &HII_density, &HM_density,
                                     &HeI_density, &HeII_density, &HeIII_density,
                                     &H2I_density, &H2II_density,
                                     &DI_density, &DII_density, &HDI_density,
                                     &ne_density, &metal_density,
                                     &temperature) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = temperature;
            break;
        case 3:  //calculate pressure
            if(calculate_pressure(&All.GrackleUnits, All.cf_atime,
                                  grid_rank, grid_dimension,
                                  grid_start, grid_end,
                                  &density, &energy,
                                  &HI_density, &HII_density, &HM_density,
                                  &HeI_density, &HeII_density, &HeIII_density,
                                  &H2I_density, &H2II_density,
                                  &DI_density, &DII_density, &HDI_density,
                                  &ne_density, &metal_density,
                                  &pressure) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = pressure;
            break;
        case 4:  //calculate gamma
            if(calculate_gamma(&All.GrackleUnits, All.cf_atime,
                               grid_rank, grid_dimension,
                               grid_start, grid_end,
                               &density, &energy,
                               &HI_density, &HII_density, &HM_density,
                               &HeI_density, &HeII_density, &HeIII_density,
                               &H2I_density, &H2II_density,
                               &DI_density, &DII_density, &HDI_density,
                               &ne_density, &metal_density,
                               &gamma) == 0) {
                fprintf(stderr, "Error in calculate_gamma.\n");
                endrun(ENDRUNVAL);
            }
            returnval = gamma;
            break;
    } //end switch
    
#else // tabular
    
    switch(mode){
        case 0:  //solve chemistry & update values (table)
            if(solve_chemistry_table(&All.GrackleUnits,
                                     All.cf_atime, dt,
                                     grid_rank, grid_dimension,
                                     grid_start, grid_end,
                                     &density, &energy,
                                     &velx, &vely, &velz,
                                     &metal_density) == 0){
                fprintf(stderr, "Error in solve_chemistry_table.\n");
                endrun(ENDRUNVAL);
            }
            double nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess;
            convert_u_to_temp(energy, rho, target, &ne_guess, &nH0_guess, &nHp_guess, &nHe0_guess, &nHep_guess, &nHepp_guess); //need to update *ne_guess for tabular!!, this may be wrong
#ifdef RT_CHEM_PHOTOION
            if(target >= 0)
            {
                SphP[target].HI = nH0_guess; SphP[target].HII = nHp_guess;
#ifdef RT_CHEM_PHOTOION_HE
                SphP[target].HeI = nHe0_guess; SphP[target].HeII = nHep_guess; SphP[target].HeIII = nHepp_guess;
#endif
            }
#endif
            returnval = energy;
            break;
        case 1:  //cooling time (table)
            if(calculate_cooling_time_table(&All.GrackleUnits,
                                            All.cf_atime,
                                            grid_rank, grid_dimension,
                                            grid_start, grid_end,
                                            &density, &energy,
                                            &velx, &vely, &velz,
                                            &metal_density,
                                            &cooling_time) == 0){
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature (table)
            if(calculate_temperature_table(&All.GrackleUnits,
                                           All.cf_atime,
                                           grid_rank, grid_dimension,
                                           grid_start, grid_end,
                                           &density, &energy,
                                           &metal_density,
                                           &temperature) == 0){
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = temperature;
            break;
        case 3:  //calculate pressure (table)
            if(calculate_pressure_table(&All.GrackleUnits,
                                        All.cf_atime,
                                        grid_rank, grid_dimension,
                                        grid_start, grid_end,
                                        &density, &energy,
                                        &pressure) == 0){
                fprintf(stderr, "Error in calculate_pressure.\n");
                endrun(ENDRUNVAL);
            }
            returnval = pressure;
            break;
    } //end switch
    
#endif // COOL_GRACKLE_CHEMISTRY
    
    return returnval;
}




//Initialize Grackle
void InitGrackle(void)
{
    grackle_verbose = 0;
    // Enable output
    if(ThisTask == 0) grackle_verbose = 1;
    
    // First, set up the units system.
    // These are conversions from code units to cgs.
    All.GrackleUnits.comoving_coordinates = 0; //All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
    All.GrackleUnits.density_units        = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
    All.GrackleUnits.length_units         = All.UnitLength_in_cm / All.HubbleParam;
    All.GrackleUnits.time_units           = All.UnitTime_in_s / All.HubbleParam;
    All.GrackleUnits.velocity_units       = All.UnitVelocity_in_cm_per_s;
    All.GrackleUnits.a_units              = 1.0; // units for the expansion factor
    
    // Second, create a chemistry object for parameters and rate data.
    if (set_default_chemistry_parameters() == 0) {
        fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
        exit(ENDRUNVAL);
    }
    // Third, set parameter values for chemistry & cooling
    
    /* optional flags:: */
    
    // Flag to control which three-body H2 formation rate is used.
    //    0: Abel, Bryan & Norman (2002),
    //    1: Palla, Salpeter & Stahler (1983),
    //    2: Cohen & Westberg (1983),
    //    3: Flower & Harris (2007),
    //    4: Glover (2008).
    //    These are discussed in Turk et. al. (2011). Default: 0.
    grackle_data.three_body_rate        = 0;
    
#ifdef METALS
    // Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    grackle_data.metal_cooling          = 1;                   // metal cooling on
    // Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the dust to gas ratio scales with the metallicity. Default: 0.
    grackle_data.h2_on_dust             = 0;                   // dust cooling/chemistry on
    // Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
    // If photoelectric_heating enabled, photoelectric_heating_rate is the heating rate in units of erg cm-3 s-1. Default: 8.5e-26.
    //   (Caution: this tends to heat gas even at extremely high densities to ~3000 K, when it should be entirely self-shielding)
    grackle_data.photoelectric_heating            = 1;         // photo-electric on [but not adjusted to local background, beware!]
    grackle_data.photoelectric_heating_rate       = 8.5e-26;
#else
    grackle_data.metal_cooling          = 0;                   // metal cooling on
    grackle_data.h2_on_dust             = 0;                   // dust cooling/chemistry off
    grackle_data.photoelectric_heating            = 0;
    grackle_data.photoelectric_heating_rate       = 8.5e-26;
#endif
    
    // Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate. Default: 1.
    grackle_data.cmb_temperature_floor  = 1;
    // Flag to enable a UV background. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    grackle_data.UVbackground           = 1;                  // UV background on
    // Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
    grackle_data.Compton_xray_heating   = 1;
    
    
    // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0.
    grackle_data.cie_cooling                      = 0;
    // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0
    grackle_data.h2_optical_depth_approximation   = 0;
    
    // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
    //    in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
    grackle_data.LWbackground_intensity           = 0;
    // Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
    //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
    grackle_data.LWbackground_sawtooth_suppression = 0;
    
    
    /* fixed flags:: */
    
    // Flag to activate the grackle machinery:
    grackle_data.use_grackle            = 1;                   // grackle on (duh)
    // Path to the data file containing the metal cooling and UV background tables:
    grackle_data.grackle_data_file      = All.GrackleDataFile; // data file
    // Flag to include radiative cooling and actually update the thermal energy during the
    // chemistry solver. If off, the chemistry species will still be updated. The most
    // common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
    grackle_data.with_radiative_cooling = 1;                   // cooling on
    // The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3.
    grackle_data.Gamma                  = GAMMA;              // our eos set in Config.sh
    // Flag to control which primordial chemistry network is used (set by Config file)
#ifndef COOL_GRACKLE_CHEMISTRY
    grackle_data.primordial_chemistry = 0;                     // fully tabulated cooling
#else
    grackle_data.primordial_chemistry = COOL_GRACKLE_CHEMISTRY;
#endif
    
    // Set initial expansion factor (for internal units).
    // Set expansion factor to 1 for non-cosmological simulation.
    double a_value = 1.0;
    if(All.ComovingIntegrationOn) a_value = All.TimeBegin;
    
    // Finally, initialize the chemistry object.
    if (initialize_chemistry_data(&All.GrackleUnits, a_value) == 0) {
        fprintf(stderr, "Error in initialize_chemistry_data.\n");
        exit(ENDRUNVAL);
    }
    
    if(ThisTask == 0)
        printf("Grackle Initialized\n");
}

#endif  //COOL_GRACKLE
