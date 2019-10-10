#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*! \file pm_nonperiodic.c
 *  \brief code for non-periodic FFT to compute long-range PM force
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * slightly by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#include "../allvars.h"
#include "../proto.h"

#ifdef PMGRID
#if !defined (BOX_PERIODIC) || defined (PM_PLACEHIGHRESREGION)

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

#ifndef GRIDBOOST
#define GRIDBOOST 2
#endif
#define  GRID  (GRIDBOOST*PMGRID)
#define  GRID2 (2*(GRID/2 + 1))


#if (GRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#define d_fftw_real fftw_real

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

static int slab_to_task[GRID];
static int *slabs_per_task;
static int *first_slab_of_task;

static int slabstart_x, nslab_x, slabstart_y, nslab_y;

static int fftsize, maxfftsize;

static fftw_real *kernel[2], *rhogrid, *forcegrid, *workspace;
static fftw_complex *fft_of_kernel[2], *fft_of_rhogrid;
static d_fftw_real *d_rhogrid, *d_forcegrid, *d_workspace;


#ifdef DM_SCALARFIELD_SCREENING
static fftw_real *kernel_scalarfield[2];
static fftw_complex *fft_of_kernel_scalarfield[2];
#endif

void pm_nonperiodic_transposeA(fftw_real * field, fftw_real * scratch);
void pm_nonperiodic_transposeB(fftw_real * field, fftw_real * scratch);
int pm_nonperiodic_compare_sortindex(const void *a, const void *b);


static struct part_slab_data
{
  large_array_offset globalindex;
  int partindex;
  int localindex;
} *part;

static int *part_sortindex;


/*! This function determines the particle extension of all particles, and for
 *  those types selected with PM_PLACEHIGHRESREGION if this is used, and then
 *  determines the boundaries of the non-periodic FFT-mesh that can be placed
 *  on this region. Note that a sufficient buffer region at the rim of the
 *  occupied part of the mesh needs to be reserved in order to allow a correct
 *  finite differencing using a 4-point formula. In addition, to allow
 *  non-periodic boundaries, the actual FFT mesh used is twice as large in
 *  each dimension compared with PMGRID.
 */
void pm_init_regionsize(void)
{
  double meshinner[2], xmin[2][3], xmax[2][3];
  int i, j;

  /* find enclosing rectangle */

  for(j = 0; j < 3; j++)
    {
      xmin[0][j] = xmin[1][j] = 1.0e36;
      xmax[0][j] = xmax[1][j] = -1.0e36;
    }

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	if(P[i].Pos[j] > xmax[0][j])
	  xmax[0][j] = P[i].Pos[j];
	if(P[i].Pos[j] < xmin[0][j])
	  xmin[0][j] = P[i].Pos[j];

#ifdef PM_PLACEHIGHRESREGION
	if(((1 << P[i].Type) & (PM_PLACEHIGHRESREGION)))
	  {
	    if(P[i].Pos[j] > xmax[1][j])
	      xmax[1][j] = P[i].Pos[j];
	    if(P[i].Pos[j] < xmin[1][j])
	      xmin[1][j] = P[i].Pos[j];
	  }
#endif
      }

  MPI_Allreduce(xmin, All.Xmintot, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, All.Xmaxtot, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for(j = 0; j < 2; j++)
    {
      All.TotalMeshSize[j] = All.Xmaxtot[j][0] - All.Xmintot[j][0];
      All.TotalMeshSize[j] = DMAX(All.TotalMeshSize[j], All.Xmaxtot[j][1] - All.Xmintot[j][1]);
      All.TotalMeshSize[j] = DMAX(All.TotalMeshSize[j], All.Xmaxtot[j][2] - All.Xmintot[j][2]);
#ifdef PM_ENLARGEREGION
      All.TotalMeshSize[j] *= PM_ENLARGEREGION;
#endif

      /* symmetrize the box onto the center */
      for(i = 0; i < 3; i++)
	{
	  All.Xmintot[j][i] = (All.Xmintot[j][i] + All.Xmaxtot[j][i]) / 2 - All.TotalMeshSize[j] / 2;
	  All.Xmaxtot[j][i] = All.Xmintot[j][i] + All.TotalMeshSize[j];
	}
    }

  /* this will produce enough room for zero-padding and buffer region to
     allow finite differencing of the potential  */

  for(j = 0; j < 2; j++)
    {
      meshinner[j] = All.TotalMeshSize[j];
      All.TotalMeshSize[j] *= 2.001 * (GRID) / ((double) (GRID - 2 - 8));
    }

  /* move lower left corner by two cells to allow finite differencing of the potential by a 4-point function */

  for(j = 0; j < 2; j++)
    for(i = 0; i < 3; i++)
      {
	All.Corner[j][i] = All.Xmintot[j][i] - 2.0005 * All.TotalMeshSize[j] / GRID;
	All.UpperCorner[j][i] = All.Corner[j][i] + (GRID / 2 - 1) * (All.TotalMeshSize[j] / GRID);
      }


#ifndef BOX_PERIODIC
  All.Asmth[0] = ASMTH * All.TotalMeshSize[0] / GRID;
  All.Rcut[0] = RCUT * All.Asmth[0];
#endif

#ifdef PM_PLACEHIGHRESREGION
  All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID;
  All.Rcut[1] = RCUT * All.Asmth[1];
#endif

#ifdef PM_PLACEHIGHRESREGION
  if(2 * All.TotalMeshSize[1] / GRID < All.Rcut[0])
    {
      All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID) / ((double) (GRID - 2));

      for(i = 0; i < 3; i++)
	{
	  All.Corner[1][i] = All.Xmintot[1][i] - 1.0001 * All.Rcut[0];
	  All.UpperCorner[1][i] = All.Corner[1][i] + (GRID / 2 - 1) * (All.TotalMeshSize[1] / GRID);
	}

      if(2 * All.TotalMeshSize[1] / GRID > All.Rcut[0])
	{
	  All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID) / ((double) (GRID - 10));

	  for(i = 0; i < 3; i++)
	    {
	      All.Corner[1][i] = All.Xmintot[1][i] - 1.0001 * (All.Rcut[0] + 2 * All.TotalMeshSize[1] / GRID);
	      All.UpperCorner[1][i] = All.Corner[1][i] + (GRID / 2 - 1) * (All.TotalMeshSize[1] / GRID);
	    }
	}

      All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID;
      All.Rcut[1] = RCUT * All.Asmth[1];
    }
#endif

  if(ThisTask == 0)
    {
#ifndef BOX_PERIODIC
      printf("Allowed region for isolated PM mesh (coarse):\n");
      printf("(%g|%g|%g)  -> (%g|%g|%g)   ext=%g  totmeshsize=%g  meshsize=%g\n\n",
	     All.Xmintot[0][0], All.Xmintot[0][1], All.Xmintot[0][2],
	     All.Xmaxtot[0][0], All.Xmaxtot[0][1], All.Xmaxtot[0][2], meshinner[0], All.TotalMeshSize[0],
	     All.TotalMeshSize[0] / GRID);
#endif
#ifdef PM_PLACEHIGHRESREGION
      printf("Allowed region for isolated PM mesh (high-res):\n");
      printf("(%g|%g|%g)  -> (%g|%g|%g)   ext=%g  totmeshsize=%g  meshsize=%g\n\n",
	     All.Xmintot[1][0], All.Xmintot[1][1], All.Xmintot[1][2],
	     All.Xmaxtot[1][0], All.Xmaxtot[1][1], All.Xmaxtot[1][2],
	     meshinner[1], All.TotalMeshSize[1], All.TotalMeshSize[1] / GRID);
#endif
    }

}

/*! Initialization of the non-periodic PM routines. The plan-files for FFTW
 *  are created. Finally, the routine to set-up the non-periodic Greens
 *  function is called.
 */
void pm_init_nonperiodic(void)
{
  int i, slab_to_task_local[GRID];
  double bytes_tot = 0;
  size_t bytes;

  /* Set up the FFTW plan files. */

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, GRID, GRID, GRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, GRID, GRID, GRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);


  for(i = 0; i < GRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, GRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  slabs_per_task = (int *) mymalloc("slabs_per_task", NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc("first_slab_of_task", NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* now allocate memory to hold the FFT fields */

#if !defined(BOX_PERIODIC)
  if(!(kernel[0] = (fftw_real *) mymalloc("kernel[0]", bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-kernel[0]' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;
  fft_of_kernel[0] = (fftw_complex *) kernel[0];
#ifdef DM_SCALARFIELD_SCREENING
  if(!
     (kernel_scalarfield[0] =
      (fftw_real *) mymalloc("kernel_scalarfield[0]", bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-kernel_scalarfield[0]' (%g MB).\n",
	     bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;
  fft_of_kernel_scalarfield[0] = (fftw_complex *) kernel_scalarfield[0];
#endif
#endif

#if defined(PM_PLACEHIGHRESREGION)
  if(!(kernel[1] = (fftw_real *) mymalloc("kernel[1]", bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-kernel[1]' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;
  fft_of_kernel[1] = (fftw_complex *) kernel[1];
#ifdef DM_SCALARFIELD_SCREENING
  if(!
     (kernel_scalarfield[1] =
      (fftw_real *) mymalloc("kernel_scalarfield[1]", bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-kernel_scalarfield[1]' (%g MB).\n",
	     bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;
  fft_of_kernel_scalarfield[1] = (fftw_complex *) kernel_scalarfield[1];
#endif
#endif

  if(ThisTask == 0)
    printf("\nAllocated %g MByte for FFT kernel(s).\n\n", bytes_tot / (1024.0 * 1024.0));

}


/*! This function allocates the workspace needed for the non-periodic FFT
 *  algorithm. Three fields are used, one for the density/potential fields,
 *  one to hold the force field obtained by finite differencing, and finally
 *  an additional workspace which is used both in the parallel FFT itself, and
 *  as a buffer for the communication algorithm.
 */
void pm_init_nonperiodic_allocate(void)
{
  double bytes_tot = 0;
  size_t bytes;

  if(!(rhogrid = (fftw_real *) mymalloc("rhogrid", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  fft_of_rhogrid = (fftw_complex *) rhogrid;

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


  if(ThisTask == 0)
    printf("Using %g MByte for non-periodic FFT computation. (presently allocated %g MB)\n",
	   bytes_tot / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));

  workspace = forcegrid;

  d_rhogrid = (d_fftw_real *) rhogrid;
  d_forcegrid = (d_fftw_real *) forcegrid;
  d_workspace = (d_fftw_real *) workspace;
}


/*! This function frees the memory allocated for the non-periodic FFT
 *  computation. (With the exception of the Greens function(s), which are kept
 *  statically in memory for the next force computation.)
 */
void pm_init_nonperiodic_free(void)
{
  /* deallocate memory */
  myfree(part_sortindex);
  myfree(part);
  myfree(forcegrid);
  myfree(rhogrid);
}


/*! This function sets-up the Greens function for the non-periodic potential
 *  in real space, and then converts it to Fourier space by means of a FFT.
 */
void pm_setup_nonperiodic_kernel(void)
{
  long i, j, k, x, y, z, ip;
  double xx, yy, zz, r, u, fac;
  double kx, ky, kz, k2, fx, fy, fz, ff;

  if(ThisTask == 0)
    printf("Setting up non-periodic PM kernel (GRID=%d)  presently allocated=%g MB).\n", (int) GRID,
	   AllocatedBytes / (1024.0 * 1024.0));

  /* now set up kernel and its Fourier transform */

  pm_init_nonperiodic_allocate();

#if !defined(BOX_PERIODIC)
  for(i = 0; i < fftsize; i++)	/* clear local density field */
    kernel[0][i] = 0;
#ifdef DM_SCALARFIELD_SCREENING
  for(i = 0; i < fftsize; i++)	/* clear local density field */
    kernel_scalarfield[0][i] = 0;
#endif

  for(i = slabstart_x; i < (slabstart_x + nslab_x); i++)
    for(j = 0; j < GRID; j++)
      for(k = 0; k < GRID; k++)
	{
	  xx = ((double) i) / GRID;
	  yy = ((double) j) / GRID;
	  zz = ((double) k) / GRID;

	  if(xx >= 0.5)
	    xx -= 1.0;
	  if(yy >= 0.5)
	    yy -= 1.0;
	  if(zz >= 0.5)
	    zz -= 1.0;

	  r = sqrt(xx * xx + yy * yy + zz * zz);

	  u = 0.5 * r / (((double) ASMTH) / GRID);

	  fac = 1 - erfc(u);

	  if(r > 0)
	    kernel[0][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] = -fac / r;
	  else
	    kernel[0][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =
	      -1 / (sqrt(M_PI) * (((double) ASMTH) / GRID));
#ifdef DM_SCALARFIELD_SCREENING
	  if(r > 0)
	    kernel_scalarfield[0][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =
	      -fac * All.ScalarBeta * exp(-r / All.ScalarScreeningLength) / r;
	  else
	    kernel_scalarfield[0][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =
	      -1 / (sqrt(M_PI) * (((double) ASMTH) / GRID));
#endif
	}

  /* do the forward transform of the kernel */

  rfftwnd_mpi(fft_forward_plan, 1, kernel[0], workspace, FFTW_TRANSPOSED_ORDER);
#ifdef DM_SCALARFIELD_SCREENING
  rfftwnd_mpi(fft_forward_plan, 1, kernel_scalarfield[0], workspace, FFTW_TRANSPOSED_ORDER);
#endif
#endif


#if defined(PM_PLACEHIGHRESREGION)
  for(i = 0; i < fftsize; i++)	/* clear local density field */
    kernel[1][i] = 0;
#ifdef DM_SCALARFIELD_SCREENING
  for(i = 0; i < fftsize; i++)	/* clear local density field */
    kernel_scalarfield[1][i] = 0;
#endif
  for(i = slabstart_x; i < (slabstart_x + nslab_x); i++)
    for(j = 0; j < GRID; j++)
      for(k = 0; k < GRID; k++)
	{
	  xx = ((double) i) / GRID;
	  yy = ((double) j) / GRID;
	  zz = ((double) k) / GRID;

	  if(xx >= 0.5)
	    xx -= 1.0;
	  if(yy >= 0.5)
	    yy -= 1.0;
	  if(zz >= 0.5)
	    zz -= 1.0;

	  r = sqrt(xx * xx + yy * yy + zz * zz);

	  u = 0.5 * r / (((double) ASMTH) / GRID);

	  fac = erfc(u * All.Asmth[1] / All.Asmth[0]) - erfc(u);

	  if(r > 0)
	    kernel[1][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] = -fac / r;
	  else
	    {
	      fac = 1 - All.Asmth[1] / All.Asmth[0];
	      kernel[1][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =
		-fac / (sqrt(M_PI) * (((double) ASMTH) / GRID));
	    }
#ifdef DM_SCALARFIELD_SCREENING
	  if(r > 0)
	    kernel_scalarfield[1][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =
	      -fac * All.ScalarBeta * exp(-r / All.ScalarScreeningLength) / r;
	  else
	    {
	      fac = 1 - All.Asmth[1] / All.Asmth[0];
	      kernel_scalarfield[1][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =
		-fac / (sqrt(M_PI) * (((double) ASMTH) / GRID));
	    }
#endif
	}

  report_memory_usage(&HighMark_pmnonperiodic, "PM_NONPERIODIC_SETUP");

  /* do the forward transform of the kernel */
  rfftwnd_mpi(fft_forward_plan, 1, kernel[1], workspace, FFTW_TRANSPOSED_ORDER);
#ifdef DM_SCALARFIELD_SCREENING
  rfftwnd_mpi(fft_forward_plan, 1, kernel_scalarfield[1], workspace, FFTW_TRANSPOSED_ORDER);
#endif
#endif

  /* deconvolve the Greens function twice with the CIC kernel */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < GRID; x++)
      for(z = 0; z < GRID / 2 + 1; z++)
	{
	  if(x > GRID / 2)
	    kx = x - GRID;
	  else
	    kx = x;
	  if(y > GRID / 2)
	    ky = y - GRID;
	  else
	    ky = y;
	  if(z > GRID / 2)
	    kz = z - GRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / GRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / GRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / GRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      ff = ff * ff * ff * ff;

	      ip = GRID * (GRID / 2 + 1) * (y - slabstart_y) + (GRID / 2 + 1) * x + z;
#if !defined(BOX_PERIODIC)
	      fft_of_kernel[0][ip].re *= ff;
	      fft_of_kernel[0][ip].im *= ff;
#ifdef DM_SCALARFIELD_SCREENING
	      fft_of_kernel_scalarfield[0][ip].re *= ff;
	      fft_of_kernel_scalarfield[0][ip].im *= ff;
#endif
#endif
#if defined(PM_PLACEHIGHRESREGION)
	      fft_of_kernel[1][ip].re *= ff;
	      fft_of_kernel[1][ip].im *= ff;
#ifdef DM_SCALARFIELD_SCREENING
	      fft_of_kernel_scalarfield[1][ip].re *= ff;
	      fft_of_kernel_scalarfield[1][ip].im *= ff;
#endif
#endif
	    }
	}
  /* end deconvolution */

  pm_init_nonperiodic_free();
}

#ifdef PM_PLACEHIGHRESREGION
int pmforce_is_particle_high_res(int type, MyDouble * Pos)
{
#ifndef SPECIAL_GAS_TREATMENT_IN_HIGHRESREGION
  /* standard treatment */
  return (1 << type) & (PM_PLACEHIGHRESREGION);
#else

  if((1 << type) & (PM_PLACEHIGHRESREGION))
    return 1;

  /* special treatment */
  int j, flag = 1;
  for(j = 0; j < 3; j++)
    if(Pos[j] < All.Xmintot[1][j] || Pos[j] > All.Xmaxtot[1][j])
      flag = 0;

  return flag;
#endif
}
#endif

/*! Calculates the long-range non-periodic forces using the PM method.  The
 *  potential is Gaussian filtered with Asmth, given in mesh-cell units. The
 *  potential is finite differenced using a 4-point finite differencing
 *  formula to obtain the force fields, which are then interpolated to the
 *  particle positions. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The CIC kernel is deconvolved.
 */
int pmforce_nonperiodic(int grnr)
{
  double dx, dy, dz;
  double fac, to_slab_fac;
  double re, im, acc_dim;
  long i, j, slab, level, sendTask, recvTask, flag, flagsum;
  long x, y, z, yl, zl, yr, zr, yll, zll, yrr, zrr, ip, dim;
  long xx, yy, zz, task, pindex;
  long slab_x, slab_y, slab_z;
  long slab_xx, slab_yy, slab_zz;
  long num_on_grid, num_field_points;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;
  MPI_Status status;

#ifdef DM_SCALARFIELD_SCREENING
  int phase;
  double kscreening2;

  kscreening2 = pow(All.BoxSize / All.ScalarScreeningLength / (2 * M_PI), 2);
#endif

#ifdef KSPACE_NEUTRINOS
  terminate("this option is not implemented here");
#endif


  if(ThisTask == 0)
    printf("Starting non-periodic PM calculation (grid=%d)  presently allocated=%g MB).\n", grnr,
	   AllocatedBytes / (1024.0 * 1024.0));

  fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3);	/* to get potential */
  fac *= 1 / (2 * All.TotalMeshSize[grnr] / GRID);	/* for finite differencing */

  to_slab_fac = GRID / All.TotalMeshSize[grnr];


  /* first, check whether all particles lie in the allowed region */

  for(i = 0, flag = 0; i < NumPart; i++)
    {
#ifdef PM_PLACEHIGHRESREGION
      if(grnr == 0 || (grnr == 1 && pmforce_is_particle_high_res(P[i].Type, P[i].Pos)))
#endif
	{
	  for(j = 0; j < 3; j++)
	    {
	      if(P[i].Pos[j] < All.Xmintot[grnr][j] || P[i].Pos[j] > All.Xmaxtot[grnr][j])
		{
		  if(flag == 0)
		    {
		      printf
			("Particle Id=%d on task=%d with coordinates (%g|%g|%g) lies outside PM mesh.\nStopping\n",
			 (int) P[i].ID, ThisTask, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }
		  flag++;
		  break;
		}
	    }
	}
    }

  MPI_Allreduce(&flag, &flagsum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(flagsum > 0)
    {
      if(ThisTask == 0)
	{
	  printf("In total %ld particles were outside allowed range.\n", flagsum);
	  fflush(stdout);
	}
      return 1;			/* error - need to return because particles were outside allowed range */
    }

  pm_init_nonperiodic_allocate();

#ifdef DM_SCALARFIELD_SCREENING
  for(phase = 0; phase < 2; phase++)
    {
#endif


      /* determine the cells each particles accesses, and how many particles lie on the grid patch */
      for(i = 0, num_on_grid = 0; i < NumPart; i++)
	{
#ifdef DM_SCALARFIELD_SCREENING
	  if(phase == 1)
	    if(P[i].Type == 0)	/* don't bin baryonic mass in this phase */
	      continue;
#endif
	  if(P[i].Pos[0] < All.Corner[grnr][0] || P[i].Pos[0] >= All.UpperCorner[grnr][0])
	    continue;
	  if(P[i].Pos[1] < All.Corner[grnr][1] || P[i].Pos[1] >= All.UpperCorner[grnr][1])
	    continue;
	  if(P[i].Pos[2] < All.Corner[grnr][2] || P[i].Pos[2] >= All.UpperCorner[grnr][2])
	    continue;

	  slab_x = (int) (to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]));
	  slab_y = (int) (to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]));
	  slab_z = (int) (to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]));

	  for(xx = 0; xx < 2; xx++)
	    for(yy = 0; yy < 2; yy++)
	      for(zz = 0; zz < 2; zz++)
		{
		  slab_xx = slab_x + xx;
		  slab_yy = slab_y + yy;
		  slab_zz = slab_z + zz;

		  offset = ((large_array_offset) GRID2) * (GRID * slab_xx + slab_yy) + slab_zz;

		  part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
		  part[num_on_grid].globalindex = offset;
		  part_sortindex[num_on_grid] = num_on_grid;
		  num_on_grid++;
		}
	}

      /* note: num_on_grid can be up to 8 times larger than the particle number,  
         but num_field_points will generally be much smaller */

      /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
      mysort_pmnonperiodic(part_sortindex, num_on_grid, sizeof(int), pm_nonperiodic_compare_sortindex);
#else
      qsort(part_sortindex, num_on_grid, sizeof(int), pm_nonperiodic_compare_sortindex);
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

	  slab = part[part_sortindex[i]].globalindex / (GRID * GRID2);
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

	  slab_x = (int) (to_slab_fac * (P[pindex].Pos[0] - All.Corner[grnr][0]));
	  slab_y = (int) (to_slab_fac * (P[pindex].Pos[1] - All.Corner[grnr][1]));
	  slab_z = (int) (to_slab_fac * (P[pindex].Pos[2] - All.Corner[grnr][2]));

	  dx = to_slab_fac * (P[pindex].Pos[0] - All.Corner[grnr][0]) - slab_x;
	  dy = to_slab_fac * (P[pindex].Pos[1] - All.Corner[grnr][1]) - slab_y;
	  dz = to_slab_fac * (P[pindex].Pos[2] - All.Corner[grnr][2]) - slab_z;

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
				   MPI_BYTE, recvTask, TAG_NONPERIOD_A, import_d_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_B, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
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
		    first_slab_of_task[ThisTask] * GRID * ((large_array_offset) GRID2);

		  d_rhogrid[offset] += import_d_data[i];
		}

	      if(level > 0)
		{
		  myfree(import_globalindex);
		  myfree(import_d_data);
		}
	    }
	}

      report_memory_usage(&HighMark_pmnonperiodic, "PM_NONPERIODIC");

      /* Do the FFT of the density field */

      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      /* multiply with the Fourier transform of the Green's function (kernel) */

      for(y = 0; y < nslab_y; y++)
	for(x = 0; x < GRID; x++)
	  for(z = 0; z < GRID / 2 + 1; z++)
	    {
	      ip = GRID * (GRID / 2 + 1) * y + (GRID / 2 + 1) * x + z;

#ifdef DM_SCALARFIELD_SCREENING
	      if(phase == 1)
		{
		  re =
		    fft_of_rhogrid[ip].re * fft_of_kernel_scalarfield[grnr][ip].re -
		    fft_of_rhogrid[ip].im * fft_of_kernel_scalarfield[grnr][ip].im;

		  im =
		    fft_of_rhogrid[ip].re * fft_of_kernel_scalarfield[grnr][ip].im +
		    fft_of_rhogrid[ip].im * fft_of_kernel_scalarfield[grnr][ip].re;
		}
	      else
#endif
		{
		  re =
		    fft_of_rhogrid[ip].re * fft_of_kernel[grnr][ip].re -
		    fft_of_rhogrid[ip].im * fft_of_kernel[grnr][ip].im;

		  im =
		    fft_of_rhogrid[ip].re * fft_of_kernel[grnr][ip].im +
		    fft_of_rhogrid[ip].im * fft_of_kernel[grnr][ip].re;
		}

	      fft_of_rhogrid[ip].re = re;
	      fft_of_rhogrid[ip].im = im;
	    }

      /* get the potential by inverse FFT */

      rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      /* Now rhogrid holds the potential */


#ifdef EVALPOTENTIAL		/* get the potential if desired */

      /* send the potential components to the right processors */

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
				   MPI_BYTE, recvTask, TAG_NONPERIOD_C, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
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
		    first_slab_of_task[ThisTask] * GRID * ((large_array_offset) GRID2);
		  import_data[i] = rhogrid[offset];
		}

	      if(level > 0)
		{
		  MPI_Sendrecv(import_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A,
			       localfield_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		  myfree(import_globalindex);
		  myfree(import_data);
		}
	    }
	}

      /* read out the potential values which all have been assembled in localfield_data */

      double pot;

      for(i = 0, j = 0; i < NumPart; i++)
	{
#ifdef PM_PLACEHIGHRESREGION
	  if(grnr == 1)
	    if(!(pmforce_is_particle_high_res(P[i].Type, P[i].Pos)))
	      continue;
#endif
	  while(j < num_on_grid && (part[j].partindex >> 3) != i)
	    j++;

	  slab_x = (int) (to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]));
	  dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;

	  slab_y = (int) (to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]));
	  dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;

	  slab_z = (int) (to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]));
	  dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

	  pot =
	    +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
	    + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
	    + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
	    + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
	    + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
	    + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
	    + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
	    + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	  P[i].PM_Potential += pot * fac * (2 * All.TotalMeshSize[grnr] / GRID);	/* compensate the finite differencing factor * */
	}
#endif


      /* get the force components by finite differencing the potential for each dimension, 
         and send back the results to the right CPUs */

      for(dim = 2; dim >= 0; dim--)	/* Calculate each component of the force. */
	{			/* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */
	  if(dim == 0)
	    pm_nonperiodic_transposeA(rhogrid, forcegrid);	/* compute the transpose of the potential field */

	  for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
	    if(xx >= 2 && xx < GRID / 2 - 2)
	      for(y = 2; y < GRID / 2 - 2; y++)
		for(z = 2; z < GRID / 2 - 2; z++)
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
			break;
		      case 2:
			zr = z + 1;
			zl = z - 1;
			zrr = z + 2;
			zll = z - 2;
			break;
		      }

		    if(dim == 0)
		      {
			forcegrid[GRID / 2 * (x + y * nslab_x) + z]
			  =
			  fac * ((4.0 / 3) *
				 (rhogrid[GRID / 2 * (x + yl * nslab_x) + zl] -
				  rhogrid[GRID / 2 * (x + yr * nslab_x) + zr]) -
				 (1.0 / 6) * (rhogrid[GRID / 2 * (x + yll * nslab_x) + zll] -
					      rhogrid[GRID / 2 * (x + yrr * nslab_x) + zrr]));
		      }
		    else

		      forcegrid[GRID2 * (GRID * x + y) + z]
			=
			fac * ((4.0 / 3) *
			       (rhogrid[GRID2 * (GRID * x + yl) + zl] -
				rhogrid[GRID2 * (GRID * x + yr) + zr]) -
			       (1.0 / 6) * (rhogrid[GRID2 * (GRID * x + yll) + zll] -
					    rhogrid[GRID2 * (GRID * x + yrr) + zrr]));
		  }

	  if(dim == 0)
	    pm_nonperiodic_transposeB(forcegrid, rhogrid);	/* compute the transpose of the potential field */

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
			(fftw_real *) mymalloc("import_data", localfield_togo[recvTask * NTask + ThisTask] *
					       sizeof(fftw_real));
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
				       recvTask, TAG_NONPERIOD_C, import_globalindex,
				       localfield_togo[recvTask * NTask +
						       sendTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
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
			first_slab_of_task[ThisTask] * GRID * ((large_array_offset) GRID2);
		      import_data[i] = forcegrid[offset];
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(import_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_A,
				   localfield_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

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
#ifdef PM_PLACEHIGHRESREGION
	      if(grnr == 1)
		if(!(pmforce_is_particle_high_res(P[i].Type, P[i].Pos)))
		  continue;
#endif
	      while(j < num_on_grid && (part[j].partindex >> 3) != i)
		j++;

	      slab_x = (int) (to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]));
	      dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;

	      slab_y = (int) (to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]));
	      dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;

	      slab_z = (int) (to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]));
	      dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

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

  pm_init_nonperiodic_free();

  if(ThisTask == 0)
    printf("done PM.\n");

  return 0;
}


void pm_nonperiodic_transposeA(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < GRID / 2; z++)
	  {
	    scratch[GRID / 2 * (first_slab_of_task[task] * nslab_x +
				x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z] =
	      field[GRID2 * (GRID * x + y) + z];
	  }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + GRID / 2 * first_slab_of_task[task] * nslab_x,
		GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + GRID / 2 * first_slab_of_task[task] * nslab_x,
		GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
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
	  MPI_Sendrecv(scratch + GRID / 2 * first_slab_of_task[task] * nslab_x,
		       GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       field + GRID / 2 * first_slab_of_task[task] * nslab_x,
		       GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif

}



void pm_nonperiodic_transposeB(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + GRID / 2 * first_slab_of_task[task] * nslab_x,
		GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + GRID / 2 * first_slab_of_task[task] * nslab_x,
		GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
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
	  MPI_Sendrecv(field + GRID / 2 * first_slab_of_task[task] * nslab_x,
		       GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       scratch + GRID / 2 * first_slab_of_task[task] * nslab_x,
		       GRID / 2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }

#endif

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < GRID / 2; z++)
	  {
	    field[GRID2 * (GRID * x + y) + z] =
	      scratch[GRID / 2 * (first_slab_of_task[task] * nslab_x +
				  x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z];
	  }
}


/* a full transpose */
/*
void pm_nonperiodic_transpose(fftw_real *field, fftw_real *scratch)
{
  MPI_Request *requests;
  int x, y, z, k, task;
  int nrequests = 0;


requests =  (MPI_Request *)  mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < GRID2; z++)
	  {
	    scratch[GRID2 *  (first_slab_of_task[task] * nslab_x   +
			      x *  slabs_per_task[task] + (y - first_slab_of_task[task])) + z] = 
	      field[GRID2 * (GRID * x + y) + z];
	  }
  

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + GRID2 * first_slab_of_task[task] * nslab_x,
		GRID2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real), 
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD,
		&requests[nrequests++]);

      MPI_Irecv(field + GRID2 * first_slab_of_task[task] * nslab_x,
		GRID2 * nslab_x * slabs_per_task[task] * sizeof(fftw_real), 
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD,
		&requests[nrequests++]);
    }


  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);


  memcpy(scratch, field, GRID2*GRID*nslab_x *sizeof(fftw_real));

  for(x = 0; x < nslab_x; x++)
    for(y = 0; y < GRID; y++)
      for(z = 0; z < GRID; z++)
	{
	  field[GRID2 * (GRID * x + y) + z] =
	    scratch[GRID2 * (x + y * nslab_x) + z];
	}
  
  myfree(requests);
}
*/



/*! Calculates the long-range non-periodic potential using the PM method.  The
 *  potential is Gaussian filtered with Asmth, given in mesh-cell units.  We
 *  carry out a CIC charge assignment, and compute the potenial by Fourier
 *  transform methods. The CIC kernel is deconvolved.
 */
int pmpotential_nonperiodic(int grnr)
{
  double dx, dy, dz;
  double fac, to_slab_fac;
  double re, im, pot;
  long i, j, slab, level, sendTask, recvTask, flag, flagsum;
  long x, y, z, ip, xx, yy, zz, task, pindex;
  long slab_x, slab_y, slab_z, slab_xx, slab_yy, slab_zz;
  long num_on_grid, num_field_points;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;
  MPI_Status status;

  if(ThisTask == 0)
    printf("Starting non-periodic PM-potential calculation (grid=%d)  presently allocated=%g MB).\n", grnr,
	   AllocatedBytes / (1024.0 * 1024.0));

  fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3);	/* to get potential */

  to_slab_fac = GRID / All.TotalMeshSize[grnr];


  /* first, check whether all particles lie in the allowed region */

  for(i = 0, flag = 0; i < NumPart; i++)
    {
#ifdef PM_PLACEHIGHRESREGION
      if(grnr == 0 || (grnr == 1 && pmforce_is_particle_high_res(P[i].Type, P[i].Pos)))
#endif
	{
	  for(j = 0; j < 3; j++)
	    {
	      if(P[i].Pos[j] < All.Xmintot[grnr][j] || P[i].Pos[j] > All.Xmaxtot[grnr][j])
		{
		  if(flag == 0)
		    {
		      printf
			("Particle Id=%d on task=%d with coordinates (%g|%g|%g) lies outside PM mesh.\nStopping\n",
			 (int) P[i].ID, ThisTask, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }
		  flag++;
		  break;
		}
	    }
	}
    }

  MPI_Allreduce(&flag, &flagsum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(flagsum > 0)
    {
      if(ThisTask == 0)
	{
	  printf("In total %ld particles were outside allowed range.\n", flagsum);
	  fflush(stdout);
	}
      return 1;			/* error - need to return because particles were outside allowed range */
    }

  pm_init_nonperiodic_allocate();

  /* determine the cells each particles accesses, and how many particles lie on the grid patch */
  for(i = 0, num_on_grid = 0; i < NumPart; i++)
    {
      if(P[i].Pos[0] < All.Corner[grnr][0] || P[i].Pos[0] >= All.UpperCorner[grnr][0])
	continue;
      if(P[i].Pos[1] < All.Corner[grnr][1] || P[i].Pos[1] >= All.UpperCorner[grnr][1])
	continue;
      if(P[i].Pos[2] < All.Corner[grnr][2] || P[i].Pos[2] >= All.UpperCorner[grnr][2])
	continue;

      slab_x = (int) (to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]));
      slab_y = (int) (to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]));
      slab_z = (int) (to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]));

      for(xx = 0; xx < 2; xx++)
	for(yy = 0; yy < 2; yy++)
	  for(zz = 0; zz < 2; zz++)
	    {
	      slab_xx = slab_x + xx;
	      slab_yy = slab_y + yy;
	      slab_zz = slab_z + zz;

	      offset = GRID2 * (GRID * slab_xx + slab_yy) + slab_zz;

	      part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
	      part[num_on_grid].globalindex = offset;
	      part_sortindex[num_on_grid] = num_on_grid;
	      num_on_grid++;
	    }
    }

  /* note: num_on_grid can be up to 8 times larger than the particle number,  
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
  mysort_pmnonperiodic(part_sortindex, num_on_grid, sizeof(int), pm_nonperiodic_compare_sortindex);
#else
  qsort(part_sortindex, num_on_grid, sizeof(int), pm_nonperiodic_compare_sortindex);
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

      slab = part[part_sortindex[i]].globalindex / (GRID * GRID2);
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

      slab_x = (int) (to_slab_fac * (P[pindex].Pos[0] - All.Corner[grnr][0]));
      slab_y = (int) (to_slab_fac * (P[pindex].Pos[1] - All.Corner[grnr][1]));
      slab_z = (int) (to_slab_fac * (P[pindex].Pos[2] - All.Corner[grnr][2]));

      dx = to_slab_fac * (P[pindex].Pos[0] - All.Corner[grnr][0]) - slab_x;
      dy = to_slab_fac * (P[pindex].Pos[1] - All.Corner[grnr][1]) - slab_y;
      dz = to_slab_fac * (P[pindex].Pos[2] - All.Corner[grnr][2]) - slab_z;

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
			       recvTask, TAG_NONPERIOD_A,
			       import_d_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_B, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
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
		import_globalindex[i] - first_slab_of_task[ThisTask] * GRID * ((large_array_offset) GRID2);

	      d_rhogrid[offset] += import_d_data[i];
	    }

	  if(level > 0)
	    {
	      myfree(import_globalindex);
	      myfree(import_d_data);
	    }
	}
    }

  report_memory_usage(&HighMark_pmnonperiodic, "PM_NONPERIODIC_POTENTIAL");

  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with the Fourier transform of the Green's function (kernel) */

  for(y = 0; y < nslab_y; y++)
    for(x = 0; x < GRID; x++)
      for(z = 0; z < GRID / 2 + 1; z++)
	{
	  ip = GRID * (GRID / 2 + 1) * y + (GRID / 2 + 1) * x + z;

	  re =
	    fft_of_rhogrid[ip].re * fft_of_kernel[grnr][ip].re -
	    fft_of_rhogrid[ip].im * fft_of_kernel[grnr][ip].im;

	  im =
	    fft_of_rhogrid[ip].re * fft_of_kernel[grnr][ip].im +
	    fft_of_rhogrid[ip].im * fft_of_kernel[grnr][ip].re;

	  fft_of_rhogrid[ip].re = re;
	  fft_of_rhogrid[ip].im = im;
	}

  /* get the potential by inverse FFT */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* Now rhogrid holds the potential */

  /* send the potential components to the right processors */

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
			       MPI_BYTE, recvTask, TAG_NONPERIOD_C, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
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
		import_globalindex[i] - first_slab_of_task[ThisTask] * GRID * ((large_array_offset) GRID2);
	      import_data[i] = rhogrid[offset];
	    }

	  if(level > 0)
	    {
	      MPI_Sendrecv(import_data,
			   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_NONPERIOD_A,
			   localfield_data + localfield_offset[recvTask],
			   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

	      myfree(import_globalindex);
	      myfree(import_data);
	    }
	}
    }


  /* read out the potential values which all have been assembled in localfield_data */

  for(i = 0, j = 0; i < NumPart; i++)
    {
#ifdef PM_PLACEHIGHRESREGION
      if(grnr == 1)
	if(!(pmforce_is_particle_high_res(P[i].Type, P[i].Pos)))
	  continue;
#endif
      while(j < num_on_grid && (part[j].partindex >> 3) != i)
	j++;

      slab_x = (int) (to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]));
      dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;

      slab_y = (int) (to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]));
      dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;

      slab_z = (int) (to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]));
      dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

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
      P[i].Potential += fac * pot;
#endif
    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm_init_nonperiodic_free();

  if(ThisTask == 0)
    printf("done PM potential.\n");

  return 0;
}


int pm_nonperiodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *) a].globalindex < part[*(int *) b].globalindex)
    return -1;

  if(part[*(int *) a].globalindex > part[*(int *) b].globalindex)
    return +1;

  return 0;
}

static void msort_pmnonperiodic_with_tmp(int *b, size_t n, int *t)
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

  msort_pmnonperiodic_with_tmp(b1, n1, t);
  msort_pmnonperiodic_with_tmp(b2, n2, t);

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

void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
  const size_t size = n * s;

  int *tmp = (int *) mymalloc("int *tmp", size);

  msort_pmnonperiodic_with_tmp((int *) b, n, tmp);

  myfree(tmp);
}


#endif
#endif
