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
	grackle_field_data my_fields;    

    // Set grid dimension and size.
    // grid_start and grid_end are used to ignore ghost zones.
    int field_size = 1;
    my_fields.grid_rank = 3;

	my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int));
    my_fields.grid_start = malloc(field_size * sizeof(int));
    my_fields.grid_end = malloc(field_size * sizeof(int));

    int i;	
    for (i = 0;i < 3;i++) 
	{
        my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
        my_fields.grid_start[i]     = 0;
        my_fields.grid_end[i]       = 0;
    }
    my_fields.grid_dimension[0] = field_size;
    my_fields.grid_end[0]       = field_size - 1;

    my_fields.x_velocity          = &SphP[target].VelPred[0];
    my_fields.y_velocity          = &SphP[target].VelPred[1];
    my_fields.z_velocity          = &SphP[target].VelPred[2];
    my_fields.density       	  = &rho;
    my_fields.internal_energy     = &u_old;
#ifdef METALS
	gr_float metal_density  = *my_fields.density * P[target].Metallicity[0];
    my_fields.metal_density = &metal_density;
#else
	gr_float metal_density  = *my_fields.density * 0.02;
    my_fields.metal_density = &metal_density;
#endif
    gr_float gamma         = GAMMA;
    
// No distinction for the function calls whether tabular or not in Grackle version 3!
    gr_float tiny = 1.0e-20;
    
    // Atomic
	gr_float e_density, HI_density, HII_density, HM_density;
	gr_float HeI_density, HeII_density, HeIII_density;
	gr_float H2I_density, H2II_density, DI_density, DII_density, HDI_density;

	e_density 	  = *my_fields.density * ne_guess;

	HI_density    = *my_fields.density * SphP[target].grHI;  //initialized with HYDROGEN_MASSFRAC
	HII_density   = *my_fields.density * SphP[target].grHII;
	HM_density    = *my_fields.density * SphP[target].grHM;

	HeI_density   = *my_fields.density * SphP[target].grHeI;
	HeII_density  = *my_fields.density * SphP[target].grHeII;
	HeIII_density = *my_fields.density * SphP[target].grHeIII;

	H2I_density   = *my_fields.density * tiny;
	H2II_density  = *my_fields.density * tiny;
	DI_density    = *my_fields.density * tiny;
	DII_density   = *my_fields.density * tiny;
	HDI_density   = *my_fields.density * tiny;

    my_fields.e_density  = &e_density;
    
    my_fields.HI_density = &HI_density;
    my_fields.HII_density = &HII_density;
    my_fields.HM_density = &HM_density;
    
    my_fields.HeI_density = &HeI_density;
    my_fields.HeII_density = &HeII_density;
    my_fields.HeIII_density = &HeIII_density;
    
    my_fields.H2I_density = &H2I_density;
    my_fields.H2II_density = &H2II_density;
    my_fields.DI_density   = &DI_density;
    my_fields.DII_density  = &DII_density;
    my_fields.HDI_density  = &HDI_density;
    
#if (COOL_GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
	H2I_density  = *my_fields.density * SphP[target].grH2I;
	H2II_density = *my_fields.density * SphP[target].grH2II;
    my_fields.H2I_density  = &H2I_density;
    my_fields.H2II_density = &H2II_density;
#endif
    
#if (COOL_GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
	DI_density  = *my_fields.density * SphP[target].grDI;
	DII_density = *my_fields.density * SphP[target].grDII;
	HDI_density = *my_fields.density * SphP[target].grHDI;
    my_fields.DI_density   = &DI_density;
    my_fields.DII_density  = &DII_density;
    my_fields.HDI_density  = &HDI_density;
#endif
    
	gr_float *cooling_time, *temperature, *pressure;
	cooling_time = malloc(field_size * sizeof(gr_float));
	temperature = malloc(field_size * sizeof(gr_float));
	pressure = malloc(field_size * sizeof(gr_float));

	// Set current redshift
	All.GrackleUnits.a_value = All.cf_atime;

    switch(mode) 
	{
        case 0:  //solve chemistry & update values
            if(solve_chemistry(&All.GrackleUnits, &my_fields, dt) == 0) 
			{
                fprintf(stderr, "Error in solve_chemistry.\n");
                endrun(ENDRUNVAL);
            }
            
            // Assign variables back
            SphP[target].grHI    = *my_fields.HI_density    / *my_fields.density;
            SphP[target].grHII   = *my_fields.HII_density   / *my_fields.density;
            SphP[target].grHM    = *my_fields.HM_density    / *my_fields.density;
            
            SphP[target].grHeI   = *my_fields.HeI_density   / *my_fields.density;
            SphP[target].grHeII  = *my_fields.HeII_density  / *my_fields.density;
            SphP[target].grHeIII = *my_fields.HeIII_density / *my_fields.density;
            
#if (COOL_GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
            SphP[target].grH2I   = *my_fields.H2I_density   / *my_fields.density;
            SphP[target].grH2II  = *my_fields.H2II_density  / *my_fields.density;
#endif
            
#if (COOL_GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
            SphP[target].grDI    = *my_fields.DI_density    / *my_fields.density;
            SphP[target].grDII   = *my_fields.DII_density   / *my_fields.density;
            SphP[target].grHDI   = *my_fields.HDI_density   / *my_fields.density;
#endif
            returnval = *my_fields.internal_energy;
            break;
            
        case 1:  //cooling time
            if(calculate_cooling_time(&All.GrackleUnits, &my_fields, cooling_time) == 0) 
			{
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *cooling_time;
            break;
        case 2:  //calculate temperature
            if(calculate_temperature(&All.GrackleUnits, &my_fields, temperature) == 0) 
			{
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *temperature;
            break;
        case 3:  //calculate pressure
            if(calculate_pressure(&All.GrackleUnits, &my_fields, pressure) == 0)
			{
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = *pressure;
            break;
        case 4:  //calculate gamma
            if(calculate_gamma(&All.GrackleUnits, &my_fields, &gamma) == 0)
			{
                fprintf(stderr, "Error in calculate_gamma.\n");
                endrun(ENDRUNVAL);
            }
            returnval = gamma;
            break;
    } //end switch
        

	// Clean up
	free(pressure);
	free(temperature);
	free(cooling_time);
	free(my_fields.grid_end);
	free(my_fields.grid_start);
	free(my_fields.grid_dimension);
	
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

	// Set initial expansion factor (for internal units).
    // Set expansion factor to 1 for non-cosmological simulation.
    All.GrackleUnits.a_value = 1.0;
    if(All.ComovingIntegrationOn) All.GrackleUnits.a_value = All.TimeBegin;
    
    // Second, create a chemistry object for parameters and rate data.
	chemistry_data *my_grackle_data;
	my_grackle_data = malloc(sizeof(chemistry_data));	
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
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
    grackle_data->three_body_rate        = 0;
    
#ifdef METALS
    // Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    grackle_data->metal_cooling          = 1;                   // metal cooling on
    // Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the dust to gas ratio scales with the metallicity. Default: 0.
    grackle_data->h2_on_dust             = 0;                   // dust cooling/chemistry on
    // Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008). Default: 0.
    // If photoelectric_heating enabled, photoelectric_heating_rate is the heating rate in units of erg cm-3 s-1. Default: 8.5e-26.
    //   (Caution: this tends to heat gas even at extremely high densities to ~3000 K, when it should be entirely self-shielding)
    grackle_data->photoelectric_heating            = 1;         // photo-electric on [but not adjusted to local background, beware!]
    grackle_data->photoelectric_heating_rate       = 8.5e-26;
#else
    grackle_data->metal_cooling          = 0;                   // metal cooling on
    grackle_data->h2_on_dust             = 0;                   // dust cooling/chemistry off
    grackle_data->photoelectric_heating            = 0;
    grackle_data->photoelectric_heating_rate       = 8.5e-26;
#endif
    
    // Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate. Default: 1.
    grackle_data->cmb_temperature_floor  = 1;
    // Flag to enable a UV background. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
    grackle_data->UVbackground           = 1;                  // UV background on
    // Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
    grackle_data->Compton_xray_heating   = 1;
    
    
    // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0.
    grackle_data->cie_cooling                      = 0;
    // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0
    grackle_data->h2_optical_depth_approximation   = 0;
    
    // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
    //    in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
    grackle_data->LWbackground_intensity           = 0;
    // Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
    //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
    grackle_data->LWbackground_sawtooth_suppression = 0;
    
    
    /* fixed flags:: */
    
    // Flag to activate the grackle machinery:
    grackle_data->use_grackle            = 1;                   // grackle on (duh)
    // Path to the data file containing the metal cooling and UV background tables:
    grackle_data->grackle_data_file      = All.GrackleDataFile; // data file
    // Flag to include radiative cooling and actually update the thermal energy during the
    // chemistry solver. If off, the chemistry species will still be updated. The most
    // common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
    grackle_data->with_radiative_cooling = 1;                   // cooling on
    // The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3.
    grackle_data->Gamma                  = GAMMA;              // our eos set in Config.sh
    // Flag to control which primordial chemistry network is used (set by Config file)
#ifndef COOL_GRACKLE_CHEMISTRY
    grackle_data->primordial_chemistry = 0;                     // fully tabulated cooling
#else
    grackle_data->primordial_chemistry = COOL_GRACKLE_CHEMISTRY;
#endif
    
	// Set number of OpenMP threads
#ifdef OPENMP
	grackle_data->omp_nthreads = OPENMP;
#endif
	    
    // Finally, initialize the chemistry object.
    if (initialize_chemistry_data(&All.GrackleUnits) == 0) {
        fprintf(stderr, "Error in initialize_chemistry_data.\n");
        exit(ENDRUNVAL);
    }
    
    if(ThisTask == 0)
        printf("Grackle Initialized\n");
}

#endif  //COOL_GRACKLE
