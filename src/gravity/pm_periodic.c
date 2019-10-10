#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*! \file pm_periodic.c
 *  \brief routines for periodic PM-force computation
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * slightly by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef PMGRID
#ifdef BOX_PERIODIC

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

#define  PMGRID2 (2*(PMGRID/2 + 1))

#if (PMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#define d_fftw_real fftw_real

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

static int slab_to_task[PMGRID];
static int *slabs_per_task;
static int *first_slab_of_task;

static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

static int fftsize, maxfftsize;

static fftw_real *rhogrid, *forcegrid, *workspace;
static d_fftw_real *d_rhogrid, *d_forcegrid, *d_workspace;
#ifdef KSPACE_NEUTRINOS
static fftw_complex *Cdata;
#endif



static fftw_complex *fft_of_rhogrid;


static MyFloat to_slab_fac;

void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch);
void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch);
int pm_periodic_compare_sortindex(const void *a, const void *b);


static struct part_slab_data
{
  large_array_offset globalindex;
  int partindex;
  int localindex;
} *part;

static int *part_sortindex;


/*! This routines generates the FFTW-plans to carry out the parallel FFTs
 *  later on. Some auxiliary variables are also initialized.
 */
void pm_init_periodic(void)
{
  int i;
  int slab_to_task_local[PMGRID];

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID; /* note that these routines REQUIRE a uniform (BOX_LONG_X=BOX_LONG_Y=BOX_LONG_Z=1) box, so we can just use 'BoxSize' */
  All.Rcut[0] = RCUT * All.Asmth[0];

  /* Set up the FFTW plan files. */

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  slabs_per_task = (int *) mymalloc("slabs_per_task", NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc("first_slab_of_task", NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = PMGRID / All.BoxSize;

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

#ifdef KSPACE_NEUTRINOS
  kspace_neutrinos_init();
#endif
}


/*! This function allocates the memory neeed to compute the long-range PM
 *  force. Three fields are used, one to hold the density (and its FFT, and
 *  then the real-space potential), one to hold the force field obtained by
 *  finite differencing, and finally a workspace field, which is used both as
 *  workspace for the parallel FFT, and as buffer for the communication
 *  algorithm used in the force computation.
 */
void pm_init_periodic_allocate(void)
{
  double bytes_tot = 0;
  size_t bytes;

  /* allocate the memory to hold the FFT fields */

  if(!(rhogrid = (fftw_real *) mymalloc("rhogrid", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(forcegrid = (fftw_real *) mymalloc("forcegrid", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-forcegrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!
     (part = (struct part_slab_data *) mymalloc("part", bytes = 8 * NumPart * sizeof(struct part_slab_data))))
    {
      printf("failed to allocate memory for `part' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(part_sortindex = (int *) mymalloc("part_sortindex", bytes = 8 * NumPart * sizeof(int))))
    {
      printf("failed to allocate memory for `part_sortindex' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  bytes_tot += bytes;


  if(ThisTask == 0)
    printf("Using %g MByte for periodic FFT computation. (presently allocated=%g MB)\n",
	   bytes_tot / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));

  workspace = forcegrid;

  fft_of_rhogrid = (fftw_complex *) & rhogrid[0];

  d_rhogrid = (d_fftw_real *) rhogrid;
  d_forcegrid = (d_fftw_real *) forcegrid;
  d_workspace = (d_fftw_real *) workspace;
}



/*! This routine frees the space allocated for the parallel FFT algorithm.
 */
void pm_init_periodic_free(void)
{
  /* allocate the memory to hold the FFT fields */
  myfree(part_sortindex);
  myfree(part);
  myfree(forcegrid);
  myfree(rhogrid);
}

#ifdef ALT_QSORT
#define KEY_TYPE int
#define KEY_BASE_TYPE large_array_offset
#define KEY_GETVAL(pk) (part[*(pk)].globalindex)
#define KEY_COPY(pk1,pk2)       \
{                               \
  *(pk2) = *(pk1);      \
}
#define QSORT qsort_pm_periodic
#include "system/myqsort.h"
#endif

/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The potential is finite differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved. Note that the particle distribution is not in the slab
 *  decomposition that is used for the FFT. Instead, overlapping patches
 *  between local domains and FFT slabs are communicated as needed.
 *
 *  For mode=0, normal force calculation, mode=1, only PS calculation.
 */
void pmforce_periodic(int mode, int *typelist)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, acc_dim;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, yl, zl, yr, zr, yll, zll, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  MyDouble pp[3], *pos;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

#ifdef DM_SCALARFIELD_SCREENING
  int phase;
  double kscreening2;

  kscreening2 = pow(All.BoxSize / All.ScalarScreeningLength / (2 * M_PI), 2);
#endif
 
  if(ThisTask == 0)
    {
      printf("Starting periodic PM calculation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
#ifndef IO_REDUCED_MODE
      fflush(stdout);
#endif
    }
    
  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */

#ifdef KSPACE_NEUTRINOS
  double rhocrit = 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
  double kspace_prefac = sqrt(pow(2 * M_PI / All.BoxSize, 3.0)) * All.OmegaNu * rhocrit * pow(All.BoxSize, 3);
#endif

  pm_init_periodic_allocate();

#ifdef DM_SCALARFIELD_SCREENING
  for(phase = 0; phase < 2; phase++)
    {
#endif

      /* determine the cells each particles accesses */
      for(i = 0, num_on_grid = 0; i < NumPart; i++)
	{
	  if(mode)
	    {
	      /* only power spectrum calculation */
	      if(typelist[P[i].Type] == 0)
		continue;
	    }

#ifdef DM_SCALARFIELD_SCREENING
	  if(phase == 1)
	    if(P[i].Type == 0)	/* don't bin baryonic mass in this phase */
	      continue;
#endif

        /* possible bugfix: Y.Feng:  (was if(mode)) */
	  if(mode > -1)
	    {
	      /* make sure that particles are properly box-wrapped */
	      for(j = 0; j < 3; j++)
		{
		  pp[j] = P[i].Pos[j];
            pp[j] = WRAP_POSITION_UNIFORM_BOX(pp[j]);
		}
	      pos = pp;
	    }
	  else
	    pos = P[i].Pos;

	  slab_x = (int) (to_slab_fac * pos[0]);
	  slab_y = (int) (to_slab_fac * pos[1]);
	  slab_z = (int) (to_slab_fac * pos[2]);

        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

	  for(xx = 0; xx < 2; xx++)
	    for(yy = 0; yy < 2; yy++)
	      for(zz = 0; zz < 2; zz++)
		{
		  slab_xx = slab_x + xx;
		  slab_yy = slab_y + yy;
		  slab_zz = slab_z + zz;

		  if(slab_xx >= PMGRID)
		    slab_xx -= PMGRID;
		  if(slab_yy >= PMGRID)
		    slab_yy -= PMGRID;
		  if(slab_zz >= PMGRID)
		    slab_zz -= PMGRID;

		  offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

		  part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
		  part[num_on_grid].globalindex = offset;
		  part_sortindex[num_on_grid] = num_on_grid;
		  num_on_grid++;
		}
	}
      /* note: num_on_grid will be  8 times larger than the particle number,
         but num_field_points will generally be much smaller */

      /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
      mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
      qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

      /* determine the number of unique field points */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  num_field_points++;
	}

      /* allocate the local field */
      localfield_globalindex =
	(large_array_offset *) mymalloc("localfield_globalindex",
					num_field_points * sizeof(large_array_offset));
      localfield_d_data =
	(d_fftw_real *) mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
      localfield_data = (fftw_real *) localfield_d_data;
      localfield_first = (int *) mymalloc("localfield_first", NTask * sizeof(int));
      localfield_count = (int *) mymalloc("localfield_count", NTask * sizeof(int));
      localfield_offset = (int *) mymalloc("localfield_offset", NTask * sizeof(int));
      localfield_togo = (int *) mymalloc("localfield_togo", NTask * NTask * sizeof(int));

      for(i = 0; i < NTask; i++)
	{
	  localfield_first[i] = 0;
	  localfield_count[i] = 0;
	}

      /* establish the cross link between the part[] array and the local list of
         mesh points. Also, count on which CPU how many of the needed field points are stored */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	      num_field_points++;

	  part[part_sortindex[i]].localindex = num_field_points;

	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

	  slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
	  task = slab_to_task[slab];
	  if(localfield_count[task] == 0)
	    localfield_first[task] = num_field_points;
	  localfield_count[task]++;
	}
      num_field_points++;

      for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
	localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

      if(mode != 0)
	{
	  foldonitself(typelist);
	  powerspec(0, typelist);
	}

      /* now bin the local particle data onto the mesh list */

      for(i = 0; i < num_field_points; i++)
	localfield_d_data[i] = 0;

      for(i = 0; i < num_on_grid; i += 8)
	{
	  pindex = (part[i].partindex >> 3);
        if(P[pindex].Mass<=0) continue;

        /* possible bugfix: Y.Feng:  (was if(mode)) */
	  if(mode > -1)
	    {
	      /* make sure that particles are properly box-wrapped */
	      for(j = 0; j < 3; j++)
		{
		  pp[j] = P[pindex].Pos[j];
            pp[j] = WRAP_POSITION_UNIFORM_BOX(pp[j]);
		}
	      pos = pp;
	    }
	  else
	    pos = P[pindex].Pos;

	  slab_x = (int) (to_slab_fac * pos[0]);
	  slab_y = (int) (to_slab_fac * pos[1]);
	  slab_z = (int) (to_slab_fac * pos[2]);
        /* possible bugfix: Y.Feng:  */
        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

	  dx = to_slab_fac * pos[0] - slab_x;
	  dy = to_slab_fac * pos[1] - slab_y;
	  dz = to_slab_fac * pos[2] - slab_z;

	  localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
	  localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
	}

      /* clear local FFT-mesh density field */
      for(i = 0; i < fftsize; i++)
	d_rhogrid[i] = 0;

      /* exchange data and add contributions to the local mesh-path */

      MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(level > 0)
		{
		  import_d_data =
		    (d_fftw_real *) mymalloc("import_d_data", localfield_togo[recvTask * NTask + ThisTask] *
					     sizeof(d_fftw_real));
		  import_globalindex =
		    (large_array_offset *) mymalloc("import_globalindex",
						    localfield_togo[recvTask * NTask +
								    ThisTask] * sizeof(large_array_offset));

		  if(localfield_togo[sendTask * NTask + recvTask] > 0
		     || localfield_togo[recvTask * NTask + sendTask] > 0)
		    {
		      MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_A, import_d_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		    }
		}
	      else
		{
		  import_d_data = localfield_d_data + localfield_offset[ThisTask];
		  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		}

	      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		{
		  /* determine offset in local FFT slab */
		  offset =
		    import_globalindex[i] -
		    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

		  d_rhogrid[offset] += import_d_data[i];
		}

	      if(level > 0)
		{
		  myfree(import_globalindex);
		  myfree(import_d_data);
		}
	    }
	}

      /* Do the FFT of the density field */

      report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC");

      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      if(mode != 0)
	{
	  powerspec(1, typelist);
	}

      if(mode == 0)		/* only carry out this part for the ordinary force calculation */
	{
	  /* multiply with Green's function for the potential */

	  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
	    for(x = 0; x < PMGRID; x++)
	      for(z = 0; z < PMGRID / 2 + 1; z++)
		{
		  if(x > PMGRID / 2)
		    kx = x - PMGRID;
		  else
		    kx = x;
		  if(y > PMGRID / 2)
		    ky = y - PMGRID;
		  else
		    ky = y;
		  if(z > PMGRID / 2)
		    kz = z - PMGRID;
		  else
		    kz = z;

		  k2 = kx * kx + ky * ky + kz * kz;

		  if(k2 > 0)
		    {
#ifdef DM_SCALARFIELD_SCREENING
		      if(phase == 1)
			smth = -All.ScalarBeta * exp(-k2 * asmth2) / (k2 + kscreening2);
		      else
#endif
			smth = -exp(-k2 * asmth2) / k2;

		      /* do deconvolution */

		      fx = fy = fz = 1;
		      if(kx != 0)
			{
			  fx = (M_PI * kx) / PMGRID;
			  fx = sin(fx) / fx;
			}
		      if(ky != 0)
			{
			  fy = (M_PI * ky) / PMGRID;
			  fy = sin(fy) / fy;
			}
		      if(kz != 0)
			{
			  fz = (M_PI * kz) / PMGRID;
			  fz = sin(fz) / fz;
			}
		      ff = 1 / (fx * fy * fz);
		      smth *= ff * ff * ff * ff;

		      /* end deconvolution */

		      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
		      fft_of_rhogrid[ip].re *= smth;
		      fft_of_rhogrid[ip].im *= smth;

#ifdef KSPACE_NEUTRINOS
		      double ampl =
			smth * kspace_prefac *
			sqrt(get_neutrino_powerspec(sqrt(k2) * 2 * M_PI / All.BoxSize, All.Time));

		      fft_of_rhogrid[ip].re += ampl * Cdata[ip].re;
		      fft_of_rhogrid[ip].im += ampl * Cdata[ip].im;
#endif
		    }
		}

	  if(slabstart_y == 0)
	    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

	  /* Do the inverse FFT to get the potential */

	  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

	  /* Now rhogrid holds the potential */

#ifdef EVALPOTENTIAL		/* now read out the potential */
	  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ level;

	      if(recvTask < NTask)
		{
		  if(level > 0)
		    {
		      import_data =
			(fftw_real *) mymalloc("import_data",
					       localfield_togo[recvTask * NTask +
							       ThisTask] * sizeof(fftw_real));
		      import_globalindex =
			(large_array_offset *) mymalloc("import_globalindex",
							localfield_togo[recvTask * NTask +
									ThisTask] *
							sizeof(large_array_offset));

		      if(localfield_togo[sendTask * NTask + recvTask] > 0
			 || localfield_togo[recvTask * NTask + sendTask] > 0)
			{
			  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				       localfield_togo[sendTask * NTask +
						       recvTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_PERIODIC_C, import_globalindex,
				       localfield_togo[recvTask * NTask +
						       sendTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
			}
		    }
		  else
		    {
		      import_data = localfield_data + localfield_offset[ThisTask];
		      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		    }

		  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		    {
		      offset =
			import_globalindex[i] -
			first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
		      import_data[i] = rhogrid[offset];
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(import_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_A,
				   localfield_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		      myfree(import_globalindex);
		      myfree(import_data);
		    }
		}
	    }

	  /* read out the potential values, which all have been assembled in localfield_data */

	  double pot;

	  for(i = 0, j = 0; i < NumPart; i++)
	    {
	      while(j < num_on_grid && (part[j].partindex >> 3) != i)
              j++;

            /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
            /* make sure that particles are properly box-wrapped */
            for(xx = 0; xx < 3; xx++)
            {
                pp[xx] = P[i].Pos[xx];
                pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
            }
            slab_x = (int) (to_slab_fac * pp[0]);
            slab_y = (int) (to_slab_fac * pp[1]);
            slab_z = (int) (to_slab_fac * pp[2]);
            dx = to_slab_fac * pp[0] - slab_x;
            dy = to_slab_fac * pp[1] - slab_y;
            dz = to_slab_fac * pp[2] - slab_z;

            pot =
            +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
            + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
            + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
            + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
            + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
            + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
            + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
            + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	      P[i].PM_Potential += pot * fac * (2 * All.BoxSize / PMGRID);
	      /* compensate the finite differencing factor */ ;
	    }

#endif


	  /* get the force components by finite differencing the potential for each dimension,
	     and send back the results to the right CPUs */

	  for(dim = 2; dim >= 0; dim--)	/* Calculate each component of the force. */
	    {			/* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */
	      if(dim == 0)
		pm_periodic_transposeA(rhogrid, forcegrid);	/* compute the transpose of the potential field */

	      for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
		for(y = 0; y < PMGRID; y++)
		  for(z = 0; z < PMGRID; z++)
		    {
		      x = xx - slabstart_x;

		      yrr = yll = yr = yl = y;
		      zrr = zll = zr = zl = z;

		      switch (dim)
			{
			case 0:	/* note: for the x-direction, we difference the transposed direction (y) */
			case 1:
			  yr = y + 1;
			  yl = y - 1;
			  yrr = y + 2;
			  yll = y - 2;
			  if(yr >= PMGRID)
			    yr -= PMGRID;
			  if(yrr >= PMGRID)
			    yrr -= PMGRID;
			  if(yl < 0)
			    yl += PMGRID;
			  if(yll < 0)
			    yll += PMGRID;
			  break;
			case 2:
			  zr = z + 1;
			  zl = z - 1;
			  zrr = z + 2;
			  zll = z - 2;
			  if(zr >= PMGRID)
			    zr -= PMGRID;
			  if(zrr >= PMGRID)
			    zrr -= PMGRID;
			  if(zl < 0)
			    zl += PMGRID;
			  if(zll < 0)
			    zll += PMGRID;
			  break;
			}

		      if(dim == 0)
			{
			  forcegrid[PMGRID * (x + y * nslab_x) + z]
			    =
			    fac * ((4.0 / 3) *
				   (rhogrid[PMGRID * (x + yl * nslab_x) + zl] -
				    rhogrid[PMGRID * (x + yr * nslab_x) + zr]) -
				   (1.0 / 6) * (rhogrid[PMGRID * (x + yll * nslab_x) + zll] -
						rhogrid[PMGRID * (x + yrr * nslab_x) + zrr]));
			}
		      else
			forcegrid[PMGRID2 * (PMGRID * x + y) + z]
			  =
			  fac * ((4.0 / 3) *
				 (rhogrid[PMGRID2 * (PMGRID * x + yl) + zl] -
				  rhogrid[PMGRID2 * (PMGRID * x + yr) + zr]) -
				 (1.0 / 6) * (rhogrid[PMGRID2 * (PMGRID * x + yll) + zll] -
					      rhogrid[PMGRID2 * (PMGRID * x + yrr) + zrr]));
		    }

	      if(dim == 0)
		pm_periodic_transposeB(forcegrid, rhogrid);	/* compute the transpose of the potential field */

	      /* send the force components to the right processors */

	      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
		{
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ level;

		  if(recvTask < NTask)
		    {
		      if(level > 0)
			{
			  import_data =
			    (fftw_real *) mymalloc("import_data",
						   localfield_togo[recvTask * NTask +
								   ThisTask] * sizeof(fftw_real));
			  import_globalindex =
			    (large_array_offset *) mymalloc("import_globalindex",
							    localfield_togo[recvTask * NTask +
									    ThisTask] *
							    sizeof(large_array_offset));

			  if(localfield_togo[sendTask * NTask + recvTask] > 0
			     || localfield_togo[recvTask * NTask + sendTask] > 0)
			    {
			      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
					   localfield_togo[sendTask * NTask +
							   recvTask] * sizeof(large_array_offset), MPI_BYTE,
					   recvTask, TAG_PERIODIC_C, import_globalindex,
					   localfield_togo[recvTask * NTask +
							   sendTask] * sizeof(large_array_offset), MPI_BYTE,
					   recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
			    }
			}
		      else
			{
			  import_data = localfield_data + localfield_offset[ThisTask];
			  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
			}

		      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
			{
			  /* determine offset in local FFT slab */
			  offset =
			    import_globalindex[i] -
			    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
			  import_data[i] = forcegrid[offset];
			}

		      if(level > 0)
			{
			  MPI_Sendrecv(import_data,
				       localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real),
				       MPI_BYTE, recvTask, TAG_PERIODIC_A,
				       localfield_data + localfield_offset[recvTask],
				       localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real),
				       MPI_BYTE, recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

			  myfree(import_globalindex);
			  myfree(import_data);
			}
		    }
		}

	      /* read out the forces, which all have been assembled in localfield_data */

	      for(i = 0, j = 0; i < NumPart; i++)
		{
#ifdef DM_SCALARFIELD_SCREENING
		  if(phase == 1)
		    if(P[i].Type == 0)	/* baryons don't get an extra scalar force */
		      continue;
#endif
            while(j < num_on_grid && (part[j].partindex >> 3) != i)
                j++;

            /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
            /* make sure that particles are properly box-wrapped */
            for(xx = 0; xx < 3; xx++)
            {
                pp[xx] = P[i].Pos[xx];
                pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
            }
            slab_x = (int) (to_slab_fac * pp[0]);
            slab_y = (int) (to_slab_fac * pp[1]);
            slab_z = (int) (to_slab_fac * pp[2]);
            dx = to_slab_fac * pp[0] - slab_x;
            dy = to_slab_fac * pp[1] - slab_y;
            dz = to_slab_fac * pp[2] - slab_z;

            acc_dim =
		    +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
		    + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
		    + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
		    + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
		    + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
		    + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
		    + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
		    + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;
            
		  P[i].GravPM[dim] += acc_dim;
		}

	    }			/* end of if(mode==0) block */

	}

      /* free locallist */
      myfree(localfield_togo);
      myfree(localfield_offset);
      myfree(localfield_count);
      myfree(localfield_first);
      myfree(localfield_d_data);
      myfree(localfield_globalindex);
#ifdef DM_SCALARFIELD_SCREENING
    }
#endif

  pm_init_periodic_free();
}


/*! Calculates the long-range potential using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 */
void pmpotential_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, pot;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MyDouble pp[3];
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("Starting periodic PM-potential calculation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
      //fflush(stdout);
    }
#endif
  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */

  pm_init_periodic_allocate();


  /* determine the cells each particles accesses */
  for(i = 0, num_on_grid = 0; i < NumPart; i++)
    {
        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

      for(xx = 0; xx < 2; xx++)
	for(yy = 0; yy < 2; yy++)
	  for(zz = 0; zz < 2; zz++)
	    {
	      slab_xx = slab_x + xx;
	      slab_yy = slab_y + yy;
	      slab_zz = slab_z + zz;

	      if(slab_xx >= PMGRID)
		slab_xx -= PMGRID;
	      if(slab_yy >= PMGRID)
		slab_yy -= PMGRID;
	      if(slab_zz >= PMGRID)
		slab_zz -= PMGRID;

	      offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

	      part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
	      part[num_on_grid].globalindex = offset;
	      part_sortindex[num_on_grid] = num_on_grid;
	      num_on_grid++;
	    }
    }

  /* note: num_on_grid will be  8 times larger than the particle number,
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
  qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

  /* determine the number of unique field points */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex =
    (large_array_offset *) mymalloc("localfield_globalindex", num_field_points * sizeof(large_array_offset));
  localfield_d_data = (d_fftw_real *) mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
  localfield_data = (fftw_real *) localfield_d_data;
  localfield_first = (int *) mymalloc("localfield_first", NTask * sizeof(int));
  localfield_count = (int *) mymalloc("localfield_count", NTask * sizeof(int));
  localfield_offset = (int *) mymalloc("localfield_offset", NTask * sizeof(int));
  localfield_togo = (int *) mymalloc("localfield_togo", NTask * NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_count[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of
     mesh points. Also, count on which CPU how many of the needed field points are stored */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	  num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

      slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
      task = slab_to_task[slab];
      if(localfield_count[task] == 0)
	localfield_first[task] = num_field_points;
      localfield_count[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

  /* now bin the local particle data onto the mesh list */

  for(i = 0; i < num_field_points; i++)
    localfield_d_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      pindex = (part[i].partindex >> 3);
        if(P[pindex].Mass<=0) continue;

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[pindex].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        dx = to_slab_fac * pp[0] - slab_x;
        dy = to_slab_fac * pp[1] - slab_y;
        dz = to_slab_fac * pp[2] - slab_z;

      localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
      localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
    }

  /* clear local FFT-mesh density field */
  for(i = 0; i < fftsize; i++)
    d_rhogrid[i] = 0;

  /* exchange data and add contributions to the local mesh-path */

  MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_d_data =
		(d_fftw_real *) mymalloc("import_d_data",
					 localfield_togo[recvTask * NTask + ThisTask] * sizeof(d_fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc("import_globalindex",
						localfield_togo[recvTask * NTask +
								ThisTask] * sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A,
			       import_d_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_B, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_d_data = localfield_d_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

	      d_rhogrid[offset] += import_d_data[i];
	    }

	  if(level > 0)
	    {
	      myfree(import_globalindex);
	      myfree(import_d_data);
	    }
	}
    }

  report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC_POTENTIAL");

  /* Do the FFT of the density field */
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2 * fac;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;

	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      fft_of_rhogrid[ip].re *= smth;
	      fft_of_rhogrid[ip].im *= smth;
	    }
	}

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the inverse FFT to get the potential */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* Now rhogrid holds the potential */


  /* now read out the potential */

  /* send the force components to the right processors */

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_data =
		(fftw_real *) mymalloc("import_data",
				       localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc("import_globalindex",
						localfield_togo[recvTask * NTask +
								ThisTask] * sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_C, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_data = localfield_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
	      import_data[i] = rhogrid[offset];
	    }

	  if(level > 0)
	    {
	      MPI_Sendrecv(import_data,
			   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_PERIODIC_A,
			   localfield_data + localfield_offset[recvTask],
			   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

	      myfree(import_globalindex);
	      myfree(import_data);
	    }
	}
    }

  /* read out the potential values, which all have been assembled in localfield_data */

  for(i = 0, j = 0; i < NumPart; i++)
    {
        while(j < num_on_grid && (part[j].partindex >> 3) != i)
            j++;

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        dx = to_slab_fac * pp[0] - slab_x;
        dy = to_slab_fac * pp[1] - slab_y;
        dz = to_slab_fac * pp[2] - slab_z;

        pot =
        +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
        + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
        + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
        + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
        + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
        + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
        + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
        + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL)
      P[i].Potential += pot;
#endif
    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm_init_periodic_free();
 
#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("done PM-Potential.\n");
      fflush(stdout);
    }
#endif
}



int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *) a].globalindex < part[*(int *) b].globalindex)
    return -1;

  if(part[*(int *) a].globalindex > part[*(int *) b].globalindex)
    return +1;

  return 0;
}

static void msort_pmperiodic_with_tmp(int *b, size_t n, int *t)
{
  int *tmp;
  int *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_pmperiodic_with_tmp(b1, n1, t);
  msort_pmperiodic_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(part[*b1].globalindex <= part[*b2].globalindex)
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(int));

  memcpy(b, t, (n - n2) * sizeof(int));
}

void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
  const size_t size = n * s;

  int *tmp = (int *) mymalloc("int *tmp", size);

  msort_pmperiodic_with_tmp((int *) b, n, tmp);

  myfree(tmp);
}

void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < PMGRID; z++)
	  {
	    scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
			      x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z] =
	      field[PMGRID2 * (PMGRID * x + y) + z];
	  }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif
}



void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }


  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < PMGRID; z++)
	  {
	    field[PMGRID2 * (PMGRID * x + y) + z] =
	      scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
				x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z];
	  }

}



#endif
#endif







#ifdef PMGRID
#ifdef BOX_PERIODIC



/*           Here comes code for the power-sepctrum computation.
 */
#define BINS_PS  2000		/* number of bins for power spectrum computation */
#define POWERSPEC_FOLDFAC 32
static long long CountModes[2][BINS_PS];
static double SumPower[2][BINS_PS];
static double SumPowerUncorrected[2][BINS_PS];	/* without binning correction (as for shot noise) */
static double Power[2][BINS_PS];
static double PowerUncorrected[2][BINS_PS];	/* without binning correction */
static double Delta[2][BINS_PS];
static double DeltaUncorrected[2][BINS_PS];	/* without binning correction */
static double ShotLimit[2][BINS_PS];
static double Kbin[BINS_PS];
static double K0, K1;
static double binfac;
static char power_spec_fname[500];
static long long power_spec_totnumpart;
static double power_spec_totmass;


void powerspec(int flag, int *typeflag)
{
  int i, n, x, y, z, kx, ky, kz, bin, ip, rep, zz;
  double k, k2, po, ponorm, smth, fac;
  double fx, fy, fz, ff, mass;
  double *powerbuf;
  long long *countbuf;
  double tstart, tend;

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("begin power spectrum. (step=%d)  POWERSPEC_FOLDFAC=%d\n", flag, POWERSPEC_FOLDFAC);
      fflush(stdout);
    }
#endif

  tstart = my_second();

  for(i = 0, mass = 0; i < NumPart; i++)
    if(typeflag[P[i].Type] && (P[i].Mass>0))
      mass += P[i].Mass;

  MPI_Allreduce(&mass, &power_spec_totmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  fac = 1.0 / power_spec_totmass;

  K0 = 2 * M_PI / All.BoxSize;	/* minimum k */
  K1 = K0 * All.BoxSize / All.SofteningTable[1];	/* maximum k */
  binfac = BINS_PS / (log(K1) - log(K0));

  if(flag == 0)
    {
      for(rep = 0; rep < 2; rep++)
	for(i = 0; i < BINS_PS; i++)
	  {
	    SumPower[rep][i] = 0;
	    SumPowerUncorrected[rep][i] = 0;
	    CountModes[rep][i] = 0;
	  }
    }

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID; z++)
	{
	  zz = z;
	  if(z >= PMGRID / 2 + 1)
	    zz = PMGRID - z;

	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      if(k2 < (PMGRID / 2.0) * (PMGRID / 2.0))
		{
		  /* do deconvolution */

		  fx = fy = fz = 1;
		  if(kx != 0)
		    {
		      fx = (M_PI * kx) / PMGRID;
		      fx = sin(fx) / fx;
		    }
		  if(ky != 0)
		    {
		      fy = (M_PI * ky) / PMGRID;
		      fy = sin(fy) / fy;
		    }
		  if(kz != 0)
		    {
		      fz = (M_PI * kz) / PMGRID;
		      fz = sin(fz) / fz;
		    }
		  ff = 1 / (fx * fy * fz);
		  smth = ff * ff * ff * ff;

		  /* end deconvolution */

		  ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + zz;

		  po = (fft_of_rhogrid[ip].re * fft_of_rhogrid[ip].re
			+ fft_of_rhogrid[ip].im * fft_of_rhogrid[ip].im);

		  po *= fac * fac * smth;

		  k = sqrt(k2) * 2 * M_PI / All.BoxSize;

		  if(flag == 0)
		    k *= POWERSPEC_FOLDFAC;

		  if(k >= K0 && k < K1)
		    {
		      bin = log(k / K0) * binfac;

		      ponorm = po / PowerSpec_Efstathiou(k);

		      SumPower[flag][bin] += ponorm;
		      SumPowerUncorrected[flag][bin] += po;


		      CountModes[flag][bin] += 1;

#ifdef RAYLEIGH_BINS
		      if(flag == 1)	/* only for course grid */
			{
			  ratio = sqrt(ponorm * RayleighFactor[bin]);

			  if(ratio > RayleighMax[bin])
			    RayleighMax[bin] = ratio;

			  if(ratio >= 10.0)
			    RayleighCountAbove[bin]++;
			  else
			    {
			      i = RAYLEIGH_BINS * ratio / 10.0;
			      RayleighCountModes[bin][i]++;
			    }
			}
#endif
		    }
		}
	    }
	}

  /* Now compute the power spectrum */

  countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));

  MPI_Allgather(CountModes[flag], BINS_PS * sizeof(long long), MPI_BYTE,
		countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      CountModes[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	CountModes[flag][i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPower[flag], BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	SumPower[flag][i] += powerbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPowerUncorrected[flag], BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPowerUncorrected[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	SumPowerUncorrected[flag][i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);


  for(i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[flag][i] > 0)
	{
	  Power[flag][i] = PowerSpec_Efstathiou(Kbin[i]) * SumPower[flag][i] / CountModes[flag][i];
	  PowerUncorrected[flag][i] = SumPowerUncorrected[flag][i] / CountModes[flag][i];
	}
      else
	{
	  Power[flag][i] = 0;
	  PowerUncorrected[flag][i] = 0;
	}

      Delta[flag][i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * Power[flag][i];

      DeltaUncorrected[flag][i] = 4 * M_PI * pow(Kbin[i], 3) /
	pow(2 * M_PI / All.BoxSize, 3) * PowerUncorrected[flag][i];

      ShotLimit[flag][i] = 4 * M_PI * pow(Kbin[i], 3) /
	pow(2 * M_PI / All.BoxSize, 3) * (1.0 / power_spec_totnumpart);
    }

  if(flag == 1)
    {
      powerspec_save();
#ifdef RAYLEIGH_BINS
      rayleigh_save();
#endif
    }

  tend = my_second();

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("end power spectrum. (step=%d) took %g seconds\n", flag, timediff(tstart, tend));
      fflush(stdout);
    }
#endif
}

double PowerSpec_Efstathiou(double k)
{
  double AA, BB, CC, nu, ShapeGamma;

  ShapeGamma = 0.21;
  AA = 6.4 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  BB = 3.0 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  CC = 1.7 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  nu = 1.13;


  return k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}



void calculate_power_spectra(int num, long long *ntot_type_all)
{
  int i, typeflag[6];

  power_spec_totnumpart = 0;

  for(i = 0; i < 6; i++)
    {
      typeflag[i] = 1;
      power_spec_totnumpart += ntot_type_all[i];
    }

  sprintf(power_spec_fname, "%s/powerspec_%03d.txt", All.OutputDir, num);

  pmforce_periodic(1, typeflag);	/* calculate power spectrum for all particle types */

#ifdef OUTPUT_POWERSPEC_EACH_TYPE
  if(ntot_type_all)
    for(i = 0; i < 6; i++)
      {
	if(ntot_type_all[i] > 0)
	  {
	    int j;
	    for(j = 0; j < 6; j++)
	      typeflag[j] = 0;

	    typeflag[i] = 1;
	    power_spec_totnumpart = ntot_type_all[i];

	    sprintf(power_spec_fname, "%s/powerspec_type%d_%03d.txt", All.OutputDir, i, num);

	    pmforce_periodic(1, typeflag);	/* calculate power spectrum for type i */
	  }
      }
#endif
}

void powerspec_save(void)
{
  FILE *fd;
  char buf[500];
  int i, flag;

  if(ThisTask == 0)
    {
      if(!(fd = fopen(power_spec_fname, "w")))
	{
	  sprintf(buf, "can't open file `%s`\n", power_spec_fname);
	  terminate(buf);
	}

      for(flag = 0; flag < 2; flag++)
	{
	  fprintf(fd, "%g\n", All.Time);
	  i = BINS_PS;
	  fprintf(fd, "%d\n", i);
	  fprintf(fd, "%g\n", power_spec_totmass);
	  fprintf(fd, "%d%09d\n", (int) (power_spec_totnumpart / 1000000000),
		  (int) (power_spec_totnumpart % 1000000000));

	  for(i = 0; i < BINS_PS; i++)
	    {
	      fprintf(fd, "%g %g %g %g %g %g %g %g %g %g\n", Kbin[i], Delta[flag][i], ShotLimit[flag][i],
		      Power[flag][i], (double) CountModes[flag][i], DeltaUncorrected[flag][i],
		      PowerUncorrected[flag][i], PowerSpec_Efstathiou(Kbin[i]), SumPower[flag][i],
		      4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3));
	    }
	}
      fclose(fd);
    }
}



void foldonitself(int *typelist)
{
  int i, j, level, sendTask, recvTask, istart, nbuf, n, rest, iter = 0;
  int slab_x, slab_xx, slab_y, slab_yy, slab_z, slab_zz;
  int *nsend_local, *nsend_offset, *nsend, count, buf_capacity;
  double to_slab_fac_folded, dx, dy, dz;
  double tstart0, tstart, tend, t0, t1;
  MyDouble pp[3];
  MyFloat *pos_sendbuf, *pos_recvbuf, *pos;
  MPI_Status status;

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("begin folding for power spectrum estimation...\n");
      fflush(stdout);
    }
#endif
  tstart0 = tstart = my_second();

  nsend_local = (int *) mymalloc("nsend_local", NTask * sizeof(int));
  nsend_offset = (int *) mymalloc("nsend_offset", NTask * sizeof(int));
  nsend = (int *) mymalloc("nsend", NTask * NTask * sizeof(int));

  buf_capacity = (maxfftsize * sizeof(d_fftw_real)) / (4 * sizeof(MyFloat));
  buf_capacity /= 2;

  pos_sendbuf = (MyFloat *) forcegrid;
  pos_recvbuf = pos_sendbuf + 4 * buf_capacity;

  to_slab_fac_folded = to_slab_fac * POWERSPEC_FOLDFAC;

  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;

  istart = 0;

  do
    {
      t0 = my_second();

      for(i = 0; i < NTask; i++)
	nsend_local[i] = 0;

      for(i = istart, nbuf = 0; i < NumPart; i++)
	{
	  if(typelist[P[i].Type] == 0)
	    continue;

	  if(nbuf + 1 >= buf_capacity)
	    break;


	  /* make sure that particles are properly box-wrapped */
	  pp[0] = P[i].Pos[0];
        pp[0] = WRAP_POSITION_UNIFORM_BOX(pp[0]);

	  slab_x = to_slab_fac_folded * pp[0];
	  slab_xx = slab_x + 1;
	  slab_x %= PMGRID;
	  slab_xx %= PMGRID;

	  nsend_local[slab_to_task[slab_x]]++;
	  nbuf++;

	  if(slab_to_task[slab_x] != slab_to_task[slab_xx])
	    {
	      nsend_local[slab_to_task[slab_xx]]++;
	      nbuf++;
	    }
	}

      for(i = 1, nsend_offset[0] = 0; i < NTask; i++)
	nsend_offset[i] = nsend_offset[i - 1] + nsend_local[i - 1];

      for(i = 0; i < NTask; i++)
	nsend_local[i] = 0;

      for(i = istart, nbuf = 0; i < NumPart; i++)
	{
	  if(typelist[P[i].Type] == 0) continue;
        if(P[i].Mass <= 0) continue;

	  if(nbuf + 1 >= buf_capacity)
	    break;

	  /* make sure that particles are properly box-wrapped */
	  pp[0] = P[i].Pos[0];
        pp[0] = WRAP_POSITION_UNIFORM_BOX(pp[0]);

	  slab_x = to_slab_fac_folded * pp[0];
	  slab_xx = slab_x + 1;
	  slab_x %= PMGRID;
	  slab_xx %= PMGRID;

	  for(j = 0; j < 3; j++)
	    pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_x]] + nsend_local[slab_to_task[slab_x]]) + j] =
	      P[i].Pos[j];

	  pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_x]] + nsend_local[slab_to_task[slab_x]]) + 3] =
	    P[i].Mass;

	  nsend_local[slab_to_task[slab_x]]++;
	  nbuf++;

	  if(slab_to_task[slab_x] != slab_to_task[slab_xx])
	    {
	      for(j = 0; j < 3; j++)
		pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_xx]] + nsend_local[slab_to_task[slab_xx]]) +
			    j] = P[i].Pos[j];

	      pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_xx]] + nsend_local[slab_to_task[slab_xx]]) +
			  3] = P[i].Mass;


	      nsend_local[slab_to_task[slab_xx]]++;
	      nbuf++;
	    }
	}

      istart = i;


      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      t1 = my_second();
#ifndef IO_REDUCED_MODE
      if(ThisTask == 0)
	{
	  printf("buffer filled (took %g sec)\n", timediff(t0, t1));
	  fflush(stdout);
	}
#endif
      t0 = my_second();
      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(recvTask != sendTask)
		{
		  MPI_Sendrecv(&pos_sendbuf[4 * nsend_offset[recvTask]],
			       4 * nsend_local[recvTask] * sizeof(MyFloat), MPI_BYTE,
			       recvTask, TAG_PM_FOLD,
			       &pos_recvbuf[0],
			       4 * nsend[recvTask * NTask + ThisTask] * sizeof(MyFloat), MPI_BYTE,
			       recvTask, TAG_PM_FOLD, MPI_COMM_WORLD, &status);

		  pos = &pos_recvbuf[0];
		  count = nsend[recvTask * NTask + ThisTask];
		}
	      else
		{
		  pos = &pos_sendbuf[4 * nsend_offset[ThisTask]];
		  count = nsend_local[ThisTask];
		}

	      for(n = 0; n < count; n++, pos += 4)
		{
		  /* make sure that particles are properly box-wrapped */
		  for(j = 0; j < 3; j++)
		    {
		      pp[j] = pos[j];
                pp[j] = WRAP_POSITION_UNIFORM_BOX(pp[j]);
		    }

		  slab_x = to_slab_fac_folded * pp[0];
		  dx = to_slab_fac_folded * pp[0] - slab_x;
		  slab_xx = slab_x + 1;
		  slab_x %= PMGRID;
		  slab_xx %= PMGRID;

		  slab_y = to_slab_fac_folded * pp[1];
		  dy = to_slab_fac_folded * pp[1] - slab_y;
		  slab_yy = slab_y + 1;
		  slab_y %= PMGRID;
		  slab_yy %= PMGRID;

		  slab_z = to_slab_fac_folded * pp[2];
		  dz = to_slab_fac_folded * pp[2] - slab_z;
		  slab_zz = slab_z + 1;
		  slab_z %= PMGRID;
		  slab_zz %= PMGRID;

		  float mass = pos[3];

		  if(slab_to_task[slab_x] == ThisTask)
		    {
		      slab_x -= first_slab_of_task[ThisTask];

		      rhogrid[(slab_x * PMGRID + slab_y) * PMGRID2 + slab_z] +=
			mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
		      rhogrid[(slab_x * PMGRID + slab_yy) * PMGRID2 + slab_z] +=
			mass * (1.0 - dx) * dy * (1.0 - dz);
		      rhogrid[(slab_x * PMGRID + slab_y) * PMGRID2 + slab_zz] +=
			mass * (1.0 - dx) * (1.0 - dy) * dz;
		      rhogrid[(slab_x * PMGRID + slab_yy) * PMGRID2 + slab_zz] += mass * (1.0 - dx) * dy * dz;
		    }

		  if(slab_to_task[slab_xx] == ThisTask)
		    {
		      slab_xx -= first_slab_of_task[ThisTask];

		      rhogrid[(slab_xx * PMGRID + slab_y) * PMGRID2 + slab_z] +=
			mass * (dx) * (1.0 - dy) * (1.0 - dz);
		      rhogrid[(slab_xx * PMGRID + slab_yy) * PMGRID2 + slab_z] +=
			mass * (dx) * dy * (1.0 - dz);
		      rhogrid[(slab_xx * PMGRID + slab_y) * PMGRID2 + slab_zz] +=
			mass * (dx) * (1.0 - dy) * dz;
		      rhogrid[(slab_xx * PMGRID + slab_yy) * PMGRID2 + slab_zz] += mass * (dx) * dy * dz;
		    }

		}
	    }
	}

      count = NumPart - istart;	/* local remaining particles */
      MPI_Allreduce(&count, &rest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      iter++;

      t1 = my_second();
#ifndef IO_REDUCED_MODE
      if(ThisTask == 0)
	{
	  printf("particles exchanged and binned. (took %g sec) max-rest=%d\n", timediff(t0, t1), rest);
	  fflush(stdout);
	}
#endif
    }
  while(rest > 0);

  tend = my_second();

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("folded density field assembled (took %g seconds, iter=%d)\n", timediff(tstart, tend), iter);
      fflush(stdout);
    }
#endif
  tstart = my_second();

  /* Do the FFT of the self-folded density field */
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  tend = my_second();

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("FFT for folded density done (took %g seconds)\n", timediff(tstart, tend));
      fflush(stdout);
    }
#endif
  myfree(nsend);
  myfree(nsend_offset);
  myfree(nsend_local);

  tend = my_second();
}


#ifdef OUTPUT_LONGRANGE_POTENTIAL
void dump_potential(void)
{
  char buf[1000];
  int nprocgroup, masterTask, groupTask, n, i, j, k;
  double asmth, fac, box, tstart, tend;
  float *potential;
  FILE *fd;

  tstart = my_second();

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("Start dumping potential\n");
      fflush(stdout);
    }
#endif
  sprintf(buf, "%s/snapdir_%03d/potential_%03d.%d", All.OutputDir, All.PowerSpecFlag - 1,
	  All.PowerSpecFlag - 1, ThisTask);

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(!(fd = fopen(buf, "w")))
	    {
	      printf("Error. Can't write in file '%s'\n", buf);
	      endrun(11);
	    }

	  n = PMGRID;
	  fwrite(&n, sizeof(int), 1, fd);

	  n = sizeof(float);
	  fwrite(&n, sizeof(int), 1, fd);

	  fwrite(&slabs_per_task[ThisTask], sizeof(int), 1, fd);
	  fwrite(&first_slab_of_task[ThisTask], sizeof(int), 1, fd);

	  box = All.BoxSize;
	  asmth = All.Asmth[0];

	  fwrite(&box, sizeof(double), 1, fd);
	  fwrite(&asmth, sizeof(double), 1, fd);

	  fac = All.G * All.PartMass / (M_PI * All.BoxSize);

	  potential = (float *) forcegrid;

	  for(i = 0; i < slabs_per_task[ThisTask]; i++)
	    for(j = 0; j < PMGRID; j++)
	      for(k = 0; k < PMGRID; k++)
		*potential++ = fac * rhogrid[(i * PMGRID + j) * PMGRID2 + k];

	  potential = (float *) forcegrid;

	  fwrite(potential, sizeof(float), PMGRID * PMGRID * slabs_per_task[ThisTask], fd);

	  fclose(fd);
	}

      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }


  MPI_Barrier(MPI_COMM_WORLD);

  tend = my_second();

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("finished writing potential (took=%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }
#endif
}
#endif



#ifdef KSPACE_NEUTRINOS
#include <gsl/gsl_rng.h>

static gsl_rng *random_generator_neutrinos;
static unsigned int *seedtable;

void kspace_neutrinos_set_seeds(void)
{
  int i, j;

  random_generator_neutrinos = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator_neutrinos, All.KspaceNeutrinoSeed);
  seedtable = mymalloc("seedtable", PMGRID * PMGRID * sizeof(unsigned int));

  for(i = 0; i < PMGRID / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[i * PMGRID + j] = 0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i + 1; j++)
	seedtable[j * PMGRID + i] = 0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i; j++)
	seedtable[(PMGRID - 1 - i) * PMGRID + j] = 0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i + 1; j++)
	seedtable[(PMGRID - 1 - j) * PMGRID + i] = 0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i; j++)
	seedtable[i * PMGRID + (PMGRID - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i + 1; j++)
	seedtable[j * PMGRID + (PMGRID - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i; j++)
	seedtable[(PMGRID - 1 - i) * PMGRID + (PMGRID - 1 - j)] =
	  0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);

      for(j = 0; j < i + 1; j++)
	seedtable[(PMGRID - 1 - j) * PMGRID + (PMGRID - 1 - i)] =
	  0x7fffffff * gsl_rng_uniform(random_generator_neutrinos);
    }
}



void kspace_neutrinos_init(void)
{
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl;
  int x, y, z, xx, yy, ip;

  init_transfer_functions();

  kspace_neutrinos_set_seeds();	/* set seeds */

  Cdata = (fftw_complex *) mymalloc("Cdata", maxfftsize * sizeof(d_fftw_real));	/* this will hold the neutrine waves */

  /* first, clean the array */

  /* note: we use TRANSPOSED_ORDER in pm_periodic, while in N-GenIC we use NORMAL_ORDER */
  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	  Cdata[ip].re = 0;
	  Cdata[ip].im = 0;
	}


  for(x = 0; x < PMGRID; x++)
    for(y = 0; y < PMGRID; y++)
      {
	gsl_rng_set(random_generator_neutrinos, seedtable[x * PMGRID + y]);

	for(z = 0; z < PMGRID / 2; z++)
	  {
	    phase = gsl_rng_uniform(random_generator_neutrinos) * 2 * M_PI;
	    do
	      ampl = gsl_rng_uniform(random_generator_neutrinos);
	    while(ampl == 0);

	    if(x == PMGRID / 2 || y == PMGRID / 2 || z == PMGRID / 2)
	      continue;
	    if(x == 0 && y == 0 && z == 0)
	      continue;

	    if(x < PMGRID / 2)
	      kvec[0] = x * 2 * M_PI / All.BoxSize;
	    else
	      kvec[0] = -(PMGRID - x) * 2 * M_PI / All.BoxSize;

	    if(y < PMGRID / 2)
	      kvec[1] = y * 2 * M_PI / All.BoxSize;
	    else
	      kvec[1] = -(PMGRID - y) * 2 * M_PI / All.BoxSize;

	    if(z < PMGRID / 2)
	      kvec[2] = z * 2 * M_PI / All.BoxSize;
	    else
	      kvec[2] = -(PMGRID - z) * 2 * M_PI / All.BoxSize;

	    kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
	    kmag = sqrt(kmag2);

	    if(All.SphereMode == 1)
	      {
		if(kmag * All.BoxSize / (2 * M_PI) > All.Nsample / 2)	/* select a sphere in k-space */
		  continue;
	      }
	    else
	      {
		if(fabs(kvec[0]) * All.BoxSize / (2 * M_PI) > All.Nsample / 2)
		  continue;
		if(fabs(kvec[1]) * All.BoxSize / (2 * M_PI) > All.Nsample / 2)
		  continue;
		if(fabs(kvec[2]) * All.BoxSize / (2 * M_PI) > All.Nsample / 2)
		  continue;
	      }

	    p_of_k = 1.0;

	    p_of_k *= -log(ampl);

	    delta = sqrt(p_of_k);

	    if(z > 0)
	      {
		if(y >= slabstart_y && y < (slabstart_y + nslab_y))
		  {
		    ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
		    Cdata[ip].re = delta * cos(phase);
		    Cdata[ip].im = delta * sin(phase);
		  }
		else		/* z=0 plane needs special treatment */
		  {
		    if(x == 0)
		      {
			if(y >= PMGRID / 2)
			  continue;
			else
			  {
			    yy = PMGRID - y;	/* note: y!=0 surely holds at this point */

			    if(y >= slabstart_y && y < (slabstart_y + nslab_y))
			      {
				ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;

				Cdata[ip].re = delta * cos(phase);
				Cdata[ip].im = delta * sin(phase);
			      }

			    if(yy >= slabstart_y && yy < (slabstart_y + nslab_y))
			      {
				ip =
				  PMGRID * (PMGRID / 2 + 1) * (yy - slabstart_y) + (PMGRID / 2 + 1) * x + z;

				Cdata[ip].re = delta * cos(phase);
				Cdata[ip].im = -delta * sin(phase);
			      }
			  }
		      }
		    else	/* here comes x!=0 : conjugate can be on other processor! */
		      {
			if(x >= PMGRID / 2)
			  continue;
			else
			  {
			    xx = PMGRID - x;
			    if(xx == PMGRID)
			      xx = 0;
			    yy = PMGRID - y;
			    if(yy == PMGRID)
			      yy = 0;

			    if(y >= slabstart_y && y < (slabstart_y + nslab_y))
			      {
				ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;

				Cdata[ip].re = delta * cos(phase);
				Cdata[ip].im = delta * sin(phase);
			      }

			    if(yy >= slabstart_y && yy < (slabstart_y + nslab_y))
			      {
				ip =
				  PMGRID * (PMGRID / 2 + 1) * (yy - slabstart_y) + (PMGRID / 2 + 1) * xx + z;

				Cdata[ip].re = delta * cos(phase);
				Cdata[ip].im = -delta * sin(phase);
			      }
			  }
		      }
		  }
	      }
	  }
      }
}

#endif /* KSPACE_NEUTRINOS */



#endif
#endif
