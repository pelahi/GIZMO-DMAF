#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../allvars.h"
#include "../proto.h"


/*! \file potential.c
 *  \brief Computation of the gravitational potential of particles
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * slightly by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#if !defined(EVALPOTENTIAL) && (defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL))

/*! This function computes the gravitational potential for ALL the particles.
 *  First, the (short-range) tree potential is computed, and then, if needed,
 *  the long range PM potential is added.
 */
void compute_potential(void)
{
  int i;

#ifndef SELFGRAVITY_OFF
  int j, k, ret, recvTask;
  int ndone, ndone_flag, dummy;
  int ngrp, place, nexport, nimport;
  double fac;
  MPI_Status status;
  double r2;

  if(All.ComovingIntegrationOn)
    set_softenings();

#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("Start computation of potential for all particles...\n");
      fflush(stdout);
    }
#endif
  CPU_Step[CPU_MISC] += measure_time();


  if(TreeReconstructFlag)
    {
#ifndef IO_REDUCED_MODE
      if(ThisTask == 0) printf("Tree construction.\n");
#endif
      CPU_Step[CPU_MISC] += measure_time();

      rearrange_particle_sequence();

      force_treebuild(NumPart, NULL);

      CPU_Step[CPU_TREEBUILD] += measure_time();

      TreeReconstructFlag = 0;

#ifndef IO_REDUCED_MODE
      if(ThisTask == 0) printf("Tree construction done.\n");
#endif
    }


  /* allocate buffers to arrange communication */
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct gravdata_in) + sizeof(struct potdata_out) +
					     sizemax(sizeof(struct gravdata_in),sizeof(struct potdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

#ifndef FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_AND_POTENTIAL
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_current != All.Ti_Current)
      drift_particle(i, All.Ti_Current);
#endif
  i = 0;			/* beginn with this index */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      for(nexport = 0; i < NumPart; i++)
	{
	  ret = force_treeevaluate_potential(i, 0, &nexport, Send_count);
	  if(ret < 0)
	    break;		/* export buffer has filled up */
	}

      MYSORT_DATAINDEX(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      GravDataGet = (struct gravdata_in *) mymalloc("GravDataGet", nimport * sizeof(struct gravdata_in));
      GravDataIn = (struct gravdata_in *) mymalloc("GravDataIn", nexport * sizeof(struct gravdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

        for(k = 0; k < 3; k++)
        {
            GravDataIn[j].Pos[k] = P[place].Pos[k];
        }

        GravDataIn[j].Type = P[place].Type;
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
        GravDataIn[j].Mass = P[place].Mass;
#endif
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        double h_place = PPP[place].Hsml;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        h_place = PPP[place].AGS_Hsml;
#endif
        if(h_place > All.ForceSoftening[P[place].Type])
        {
            GravDataIn[j].Soft = h_place;
        } else {
            GravDataIn[j].Soft = All.ForceSoftening[P[place].Type];
        }
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(RT_USE_GRAVTREE)
        if((P[place].Type == 0) && (PPP[place].Hsml > All.ForceSoftening[P[place].Type]))
        {
            GravDataIn[j].Soft = PPP[place].Hsml;
        } else {
            GravDataIn[j].Soft = All.ForceSoftening[P[place].Type];
        }
#endif
        GravDataIn[j].OldAcc = P[place].OldAcc;
        
        for(k = 0; k < NODELISTLENGTH; k++)
        {
            GravDataIn[j].NodeList[k] = DataNodeList[DataIndexTable[j].IndexGet].NodeList[k];
        }
	}


      /* exchange particle data */

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&GravDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
			       recvTask, TAG_POTENTIAL_A,
			       &GravDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
			       recvTask, TAG_POTENTIAL_A, MPI_COMM_WORLD, &status);
		}
	    }
	}

      myfree(GravDataIn);
      PotDataResult = (struct potdata_out *) mymalloc("PotDataResult", nimport * sizeof(struct potdata_out));
      PotDataOut = (struct potdata_out *) mymalloc("PotDataOut", nexport * sizeof(struct potdata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	{
	  force_treeevaluate_potential(j, 1, &dummy, &dummy);
	}

      if(i >= NumPart)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&PotDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct potdata_out),
			       MPI_BYTE, recvTask, TAG_POTENTIAL_B,
			       &PotDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct potdata_out),
			       MPI_BYTE, recvTask, TAG_POTENTIAL_B, MPI_COMM_WORLD, &status);
		}
	    }

	}

      /* add the results to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  P[place].Potential += PotDataOut[j].Potential;
	}

      myfree(PotDataOut);
      myfree(PotDataResult);
      myfree(GravDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);

#ifndef ADAPTIVE_GRAVSOFT_FORALL
    /* add correction to exclude self-potential */
    for(i = 0; i < NumPart; i++)
    {
        /* remove self-potential */
        P[i].Potential += P[i].Mass / All.SofteningTable[P[i].Type];
        
#ifdef BOX_PERIODIC
        if(All.ComovingIntegrationOn)
            P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
            pow(All.Omega0 * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G), 1.0 / 3);
#endif
    }
#endif

  /* multiply with the gravitational constant */

  for(i = 0; i < NumPart; i++)
    P[i].Potential *= All.G;


#ifdef PMGRID

#ifdef BOX_PERIODIC
  pmpotential_periodic();
#ifdef PM_PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(1);	/* try again */
    }
  if(i == 1)
    endrun(88686);
#endif
#else
  i = pmpotential_nonperiodic(0);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(0);	/* try again */
    }
  if(i == 1)
    endrun(88687);
#ifdef PM_PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();

      i = pmpotential_nonperiodic(1);
    }
  if(i != 0)
    endrun(88688);
#endif
#endif

#endif



  if(All.ComovingIntegrationOn)
    {
#ifndef BOX_PERIODIC
      fac = -0.5 * All.Omega0 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits;

      for(i = 0; i < NumPart; i++)
	{
	  for(k = 0, r2 = 0; k < 3; k++)
	    r2 += P[i].Pos[k] * P[i].Pos[k];

	  P[i].Potential += fac * r2;
	}
#endif
    }
  else
    {
      fac = -0.5 * All.OmegaLambda * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits;
      if(fac != 0)
	{
	  for(i = 0; i < NumPart; i++)
	    {
	      for(k = 0, r2 = 0; k < 3; k++)
		r2 += P[i].Pos[k] * P[i].Pos[k];

	      P[i].Potential += fac * r2;
	    }
	}
    }


#ifndef IO_REDUCED_MODE
  if(ThisTask == 0)
    {
      printf("potential done.\n");
      fflush(stdout);
    }
#endif

#else
  for(i = 0; i < NumPart; i++)
    P[i].Potential = 0;
#endif

  CPU_Step[CPU_POTENTIAL] += measure_time();
}


#endif
