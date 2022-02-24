#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef USE_FFTW3
#include <fftw3-mpi.h>
#elif USE_PFFT
#else
#ifdef SINGLEPRECISION_FFTW2
#include <srfftw_mpi.h>
#else
#include <rfftw_mpi.h>
#endif
#endif

#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#include <time.h>

void displacement_fields(void)
{
  MPI_Request request;
  MPI_Status status;
  gsl_rng *random_generator;
  int i, j, k, ii, jj, kk, axes;
  int n;
  int sendTask, recvTask;
  double fac, vel_prefac, vel_prefac2;
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl, hubble_a;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis, dis2, maxdisp, meandisp, maxdisp2, maxvel, meanvel, meandisp2, max_disp_glob, max_disp2_glob, max_vel_glob, mean_vel_glob, mean_disp_glob,mean_disp2_glob;
  unsigned int *seedtable;

  time_t time_x;

  time_x = time(NULL);

  unsigned int bytes,dimsize;
  int coord;
  fftw_complex *(cdisp[3]), *(cdisp2[3]) ; /* ZA and 2nd order displacements */
  fftw_complex *(cdigrad[6]);
#ifdef USE_PFFT
#elif USE_FFTW3
  double *(disp[3]), *(disp2[3]) ;
  double *(digrad[6]);
#else
  fftw_real *(disp[3]), *(disp2[3]) ;
  fftw_real *(digrad[6]);
#endif
	
#ifdef CORRECT_CIC
  double fx, fy, fz, ff, smth;
#endif

  FILE *fd;
  fd=fopen("out.txt","w");
  
  if(ThisTask == 0)
    {
      printf("\nComputing displacement fields...\n");
      fflush(stdout);
    }

  hubble_a =
    Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);
  
  vel_prefac = InitTime * hubble_a * F_Omega(InitTime);
  vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime);
  
  vel_prefac /= sqrt(InitTime);	/* converts to Gadget velocity */
  vel_prefac2 /= sqrt(InitTime);

  fac = pow(2 * PI / Box, 1.5);
  
  if(ThisTask == 0)
    {
      fprintf(stdout,"vel_prefac= %g, vel_prefac2= %g,  hubble_a=%g fom=%g fac=%g\n", vel_prefac, vel_prefac2, hubble_a, F_Omega(InitTime),fac);
      fflush(stdout);
    }
  
  maxdisp  = 0;
  meandisp = 0;
  maxdisp2 = 0;
  maxvel   = 0;
  meanvel  = 0;
  meandisp2= 0;

  // Generate random seeds
  
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  
  gsl_rng_set(random_generator, Seed);

  if(!(seedtable = malloc(Nmesh * Nmesh * sizeof(unsigned int))))
    FatalError(4);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
      
      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }

  if(ThisTask==0)
    {
      fprintf(stdout,"Finished seed table...\n");
      fflush(stdout);
    }
  
  for(axes=0,bytes=0,dimsize=0; axes < 3; axes++)
    {
#ifdef USE_FFTW3
      bytes+=sizeof(double)*TotalSizePlusAdditional;
      dimsize += TotalSizePlusAdditional;
      cdisp[axes] = fftw_alloc_complex(dimsize);
      cdisp[axes] = (fftw_complex *) malloc(bytes += sizeof(double) * TotalSizePlusAdditional);
      disp[axes] = (double *) cdisp[axes];
#else
      cdisp[axes] = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      disp[axes] = (fftw_real *) cdisp[axes];
#endif
    }
	
  ASSERT_ALLOC(cdisp[0] && cdisp[1] && cdisp[2]);
	
#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {
      // Initialise the arrays
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    for(axes = 0; axes < 3; axes++)
	      {
#ifdef USE_FFTW3
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][0] = 0;
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][1] = 0;
#else
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
#endif
	      }

      if(ThisTask==0)
	{
	  fprintf(stdout,"Initialised displacement arrays...\n");
	  fflush(stdout);
	}
      
      for(i = 0; i < Nmesh; i++)
	{
	  ii = Nmesh - i;
	  if(ii == Nmesh)
	    ii = 0;
	  if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
	     (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
	    {
	      for(j = 0; j < Nmesh; j++)
		{
		  gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);
		  
		  for(k = 0; k < Nmesh / 2; k++)
		    {
		      phase = gsl_rng_uniform(random_generator) * 2 * PI;
		      do
			ampl = gsl_rng_uniform(random_generator);
		      while(ampl == 0);
		      
		      if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
			continue;
		      if(i == 0 && j == 0 && k == 0)
			continue;
		      
		      if(i < Nmesh / 2)
			kvec[0] = i * 2 * PI / Box;
		      else
			kvec[0] = -(Nmesh - i) * 2 * PI / Box;
		      
		      if(j < Nmesh / 2)
			kvec[1] = j * 2 * PI / Box;
		      else
			kvec[1] = -(Nmesh - j) * 2 * PI / Box;
		      
		      if(k < Nmesh / 2)
			kvec[2] = k * 2 * PI / Box;
		      else
			kvec[2] = -(Nmesh - k) * 2 * PI / Box;
		      
		      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
		      kmag = sqrt(kmag2);
		      
		      if(SphereMode == 1)
			{
			  if(kmag * Box / (2 * PI) > Nsample / 2)	/* select a sphere in k-space */
			    continue;
			}
		      else
			{
			  if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			}
		      
		      p_of_k = PowerSpec(kmag);
		      
		      p_of_k *= -log(ampl);
		      
		      delta = fac * sqrt(p_of_k) / Dplus;	/* scale back to starting redshift */

		      if(k > 0)
			{
			  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
			    for(axes = 0; axes < 3; axes++)
			      {
#ifdef USE_FFTW3
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k][0] =
				  -kvec[axes] / kmag2 * delta * sin(phase);
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k][1] =
				  kvec[axes] / kmag2 * delta * cos(phase);
#else
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				  -kvec[axes] / kmag2 * delta * sin(phase);
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				  kvec[axes] / kmag2 * delta * cos(phase);
#endif
				
			      }
			}
		      else	/* k=0 plane needs special treatment */
			{
			  if(i == 0)
			    {
			      if(j >= Nmesh / 2)
				continue;
			      else
				{
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    {
				      jj = Nmesh - j;	/* note: j!=0 surely holds at this point */
				      
				      for(axes = 0; axes < 3; axes++)
					{
#ifdef USE_FFTW3
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k][0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k][1] =
					    kvec[axes] / kmag2 * delta * cos(phase);
					  
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k][0] =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k][1] =
					    -kvec[axes] / kmag2 * delta * cos(phase);
#else
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase);
					  
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					    -kvec[axes] / kmag2 * delta * cos(phase);
#endif
					}
				    }
				}
			    }
			  else	/* here comes i!=0 : conjugate can be on other processor! */
			    {
			      if(i >= Nmesh / 2)
				continue;
			      else
				{
				  ii = Nmesh - i;
				  if(ii == Nmesh)
				    ii = 0;
				  jj = Nmesh - j;
				  if(jj == Nmesh)
				    jj = 0;
				  
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
#ifdef USE_FFTW3
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k][0] =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k][1] =
					  kvec[axes] / kmag2 * delta * cos(phase);
#else
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					  kvec[axes] / kmag2 * delta * cos(phase);
#endif
				      }
				  
				  if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
#ifdef USE_FFTW3
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k][0] =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k][1] =
					  -kvec[axes] / kmag2 * delta * cos(phase);
#else
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					  -kvec[axes] / kmag2 * delta * cos(phase);
#endif
				      }
				}
			    }
			}
		    }
		}
	    }
	}
      
      /* At this point, cdisp[axes] contains the complex Zeldovich displacement */
      
      if(ThisTask == 0) printf("Done Zeldovich.\n");
      
      /* Compute displacement gradient */
      
      for(i = 0; i < 6; i++)
	{
#ifdef USE_FFTW3
	  bytes = sizeof(double) * TotalSizePlusAdditional;	  
	  cdigrad[i] = fftw_alloc_complex(TotalSizePlusAdditional);
	  digrad[i] = (double *) cdigrad[i];
#else
	  cdigrad[i] = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  digrad[i] = (fftw_real *) cdigrad[i];
#endif
	  ASSERT_ALLOC(cdigrad[i]);
	  
          if(ThisTask==0)printf("\n Allocating memory for vector number ..%i\n", i);
	}
      
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh/2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
	      
	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;
	      
	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;
	      
	      /* Derivatives of ZA displacement  */
	      /* d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i */
#ifdef USE_FFTW3	      
	      cdigrad[0][coord][0] = -cdisp[0][coord][1] * kvec[0]; /* disp0,0 */
	      cdigrad[0][coord][1] = cdisp[0][coord][0] * kvec[0];
	      
	      cdigrad[1][coord][0] = -cdisp[0][coord][1] * kvec[1]; /* disp0,1 */
	      cdigrad[1][coord][1] = cdisp[0][coord][0] * kvec[1];
	      
	      cdigrad[2][coord][0] = -cdisp[0][coord][1] * kvec[2]; /* disp0,2 */
	      cdigrad[2][coord][1] = cdisp[0][coord][0] * kvec[2];
	      
	      cdigrad[3][coord][0] = -cdisp[1][coord][1] * kvec[1]; /* disp1,1 */
	      cdigrad[3][coord][1] = cdisp[1][coord][0] * kvec[1];
	      
	      cdigrad[4][coord][0] = -cdisp[1][coord][1] * kvec[2]; /* disp1,2 */
	      cdigrad[4][coord][1] = cdisp[1][coord][0] * kvec[2];
	      
	      cdigrad[5][coord][0] = -cdisp[2][coord][1] * kvec[2]; /* disp2,2 */
	      cdigrad[5][coord][1] = cdisp[2][coord][0] * kvec[2];
#else
	      cdigrad[0][coord].re = -cdisp[0][coord].im * kvec[0]; /* disp0,0 */
	      cdigrad[0][coord].im = cdisp[0][coord].re * kvec[0];
	      
	      cdigrad[1][coord].re = -cdisp[0][coord].im * kvec[1]; /* disp0,1 */
	      cdigrad[1][coord].im = cdisp[0][coord].re * kvec[1];
	      
	      cdigrad[2][coord].re = -cdisp[0][coord].im * kvec[2]; /* disp0,2 */
	      cdigrad[2][coord].im = cdisp[0][coord].re * kvec[2];
	      
	      cdigrad[3][coord].re = -cdisp[1][coord].im * kvec[1]; /* disp1,1 */
	      cdigrad[3][coord].im = cdisp[1][coord].re * kvec[1];
	      
	      cdigrad[4][coord].re = -cdisp[1][coord].im * kvec[2]; /* disp1,2 */
	      cdigrad[4][coord].im = cdisp[1][coord].re * kvec[2];
	      
	      cdigrad[5][coord].re = -cdisp[2][coord].im * kvec[2]; /* disp2,2 */
	      cdigrad[5][coord].im = cdisp[2][coord].re * kvec[2];	      
#endif
	    }      
      
      if(ThisTask==0) printf("\n Begin Fourier transforming displacement gradient \n");
      //get_time(time_x);
      
      if(ThisTask == 0) printf("Fourier transforming displacement gradient...");

#ifdef USE_FFTW3
      for(i = 0; i < 6; i++)
	{
	  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdigrad[i],digrad[i],
					  MPI_COMM_WORLD, FFTW_ESTIMATE);      
	  fftw_execute(Inverse_plan);
	  if(ThisTask==0) printf("\n Done Fourier transforming displacement gradient number  %i\n", i);
	}
#else
      for(i = 0; i < 6; i++)
	{
	  rfftwnd_mpi(Inverse_plan, 1, digrad[i], Workspace, FFTW_NORMAL_ORDER);
	  if(ThisTask==0) printf("\n Done Fourier transforming displacement gradient number  %i\n", i);
	}
#endif
      
      /* Compute second order source and store it in digrad[3]*/
      
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k < Nmesh; k++)
	    {
	      coord = (i * Nmesh + j) * (2*(Nmesh / 2 + 1)) + k;
	      
	      digrad[3][coord] =
		digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])
		+digrad[3][coord]*digrad[5][coord]
                -digrad[1][coord]*digrad[1][coord]
		-digrad[2][coord]*digrad[2][coord]
		-digrad[4][coord]*digrad[4][coord];
	    }
      
      
      if(ThisTask == 0) printf("Fourier transforming second order source...");
#ifdef USE_FFTW3
      Forward_plan = fftw_mpi_plan_dft_r2c_3d(Nmesh,Nmesh,Nmesh,digrad[3],cdigrad[3],
					      MPI_COMM_WORLD, FFTW_ESTIMATE);      
      fftw_execute(Forward_plan);
#else
      rfftwnd_mpi(Forward_plan, 1, digrad[3], Workspace, FFTW_NORMAL_ORDER);
#endif
      if(ThisTask == 0) printf("Done.\n");
      
      if(ThisTask==0)printf("\n Done the next Fourier transform in ...i\n");
	//get_time(time_x);
      
      /* The memory allocated for cdigrad[0], [1], and [2] will be used for 2nd order displacements */
      /* Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later */
      
      for(axes = 0; axes < 3; axes++) 
	{
	  cdisp2[axes] = cdigrad[axes];
#ifdef USE_FFTW3
	  disp2[axes] = (double *) cdisp2[axes];
#else
	  disp2[axes] = (fftw_real *) cdisp2[axes];
#endif
	}
      
      free(cdigrad[4]); free(cdigrad[5]); 
      
      /* Solve Poisson eq. and calculate 2nd order displacements */
      
      for(i = 0; i < Local_nx; i++){
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
	      
	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;
	      
	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;
	      
	      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
#ifdef CORRECT_CIC
	      /* calculate smooth factor for deconvolution of CIC interpolation */
	      fx = fy = fz = 1;
	      if(kvec[0] != 0)
		{
		  fx = (kvec[0] * Box / 2) / Nmesh;
		  fx = sin(fx) / fx;
		}
	      if(kvec[1] != 0)
		{
		  fy = (kvec[1] * Box / 2) / Nmesh;
		  fy = sin(fy) / fy;
		}
	      if(kvec[2] != 0)
		{
		  fz = (kvec[2] * Box / 2) / Nmesh;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth = ff * ff;
#endif
	      /* cdisp2 = source * k / (sqrt(-1) k^2) */
	      for(axes = 0; axes < 3; axes++)
		{
		  if(kmag2 > 0.0) 
		    {
#ifdef USE_FFTW3
		      cdisp2[axes][coord][0] = cdigrad[3][coord][1] * kvec[axes] / kmag2;
		      cdisp2[axes][coord][1] = -cdigrad[3][coord][0] * kvec[axes] / kmag2;
#else
		      cdisp2[axes][coord].re = cdigrad[3][coord].im * kvec[axes] / kmag2;
		      cdisp2[axes][coord].im = -cdigrad[3][coord].re * kvec[axes] / kmag2;
#endif
		    }
		  else
#ifdef USE_FFTW3
		    cdisp2[axes][coord][0] = cdisp2[axes][coord][1] = 0.0;
#ifdef CORRECT_CIC
		  cdisp[axes][coord][0] *= smth;   cdisp[axes][coord][1] *= smth;
		  cdisp2[axes][coord][0] *= smth;  cdisp2[axes][coord][1] *= smth;
#endif
#else
		  cdisp2[axes][coord].re = cdisp2[axes][coord].im = 0.0;
#ifdef CORRECT_CIC
		  cdisp[axes][coord].re *= smth;   cdisp[axes][coord].im *= smth;
		  cdisp2[axes][coord].re *= smth; cdisp2[axes][coord].im *= smth;
#endif
#endif
		}
	    }

	//if(ThisTask==0)printf("\n Calculating 2nd order displacement i in local_nx %i\n", i);
	//get_time(time_x);

      }
      /* Free cdigrad[3] */
      free(cdigrad[3]);
      
      /* Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements */

      for(axes = 0; axes < 3; axes++)
	{
          if(ThisTask == 0) printf("Fourier transforming displacements, axis %d.\n",axes);

#ifdef USE_FFTW3
	  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp[axes],disp[axes],
						  MPI_COMM_WORLD, FFTW_ESTIMATE);      
	  fftw_execute(Inverse_plan);
	  Inverse_plan = fftw_mpi_plan_dft_c2r_3d(Nmesh,Nmesh,Nmesh,cdisp2[axes],disp2[axes],
						  MPI_COMM_WORLD, FFTW_ESTIMATE);      	  
	  fftw_execute(Inverse_plan);
#else
	  rfftwnd_mpi(Inverse_plan, 1, disp[axes], Workspace, FFTW_NORMAL_ORDER);
	  rfftwnd_mpi(Inverse_plan, 1, disp2[axes], Workspace, FFTW_NORMAL_ORDER);
#endif
	  /* now get the plane on the right side from neighbour on the right, 
	     and send the left plane */
	  
	  recvTask = ThisTask;
	  do
	    {
	      recvTask--;
	      if(recvTask < 0)
		recvTask = NTask - 1;
	    }
	  while(Local_nx_table[recvTask] == 0);
	  
	  sendTask = ThisTask;
	  do
	    {
	      sendTask++;
	      if(sendTask >= NTask)
		sendTask = 0;
	    }
	  while(Local_nx_table[sendTask] == 0);
	  
	  /* use non-blocking send */
	  
	  if(Local_nx > 0)
	    {
	      /* send ZA disp */
#ifdef USE_FFTW3
	      MPI_Isend(&(disp[axes][0]),
			sizeof(double) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(double) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
#else
	      MPI_Isend(&(disp[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
#endif
	      MPI_Wait(&request, &status);
	      
	      
	      /* send 2nd order disp */
#ifdef USE_FFTW3
	      MPI_Isend(&(disp2[axes][0]),
			sizeof(double) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp2[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(double) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
#else
	      MPI_Isend(&(disp2[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp2[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
#endif
	      MPI_Wait(&request, &status);
	    }
	  
          //if(ThisTask==0)printf("\n Fourier transforming both displacements back on axes.... %i\n", axes);
          //get_time(time_x);
	}

#ifdef USE_FFTW3      
      fftw_destroy_plan(Inverse_plan);
      fftw_destroy_plan(Forward_plan);
#else
      fftw_free(Workspace);
      rfftwnd_mpi_destroy_plan(Inverse_plan);
      rfftwnd_mpi_destroy_plan(Forward_plan);      
#endif      
      
      /* read-out displacements */
      long long nmesh3 = ((long long)Nmesh ) * ((long long)Nmesh) *((long long)Nmesh);
      
      for(n = 0; n < NumPart; n++)
	{
#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
	  if(P[n].Type == Type)
#endif
	    {
	      u = P[n].Pos[0] / Box * Nmesh;
	      v = P[n].Pos[1] / Box * Nmesh;
	      w = P[n].Pos[2] / Box * Nmesh;
	      
	      i = (int) u;
	      j = (int) v;
	      k = (int) w;
	      
	      if(i == (Local_x_start + Local_nx))
		i = (Local_x_start + Local_nx) - 1;
	      if(i < Local_x_start)
		i = Local_x_start;
	      if(j == Nmesh)
		j = Nmesh - 1;
	      if(k == Nmesh)
		k = Nmesh - 1;
	      
	      u -= i;
	      v -= j;
	      w -= k;
	      
	      i -= Local_x_start;
	      ii = i + 1;
	      jj = j + 1;
	      kk = k + 1;
	      
	      if(jj >= Nmesh)
		jj -= Nmesh;
	      if(kk >= Nmesh)
		kk -= Nmesh;
	      
	      f1 = (1 - u) * (1 - v) * (1 - w);
	      f2 = (1 - u) * (1 - v) * (w);
	      f3 = (1 - u) * (v) * (1 - w);
	      f4 = (1 - u) * (v) * (w);
	      f5 = (u) * (1 - v) * (1 - w);
	      f6 = (u) * (1 - v) * (w); 
	      f7 = (u) * (v) * (1 - w);
	      f8 = (u) * (v) * (w);
	     
	      for(axes = 0; axes < 3; axes++)
		{
		  dis = disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

		  dis2 = disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;
		  dis2 /= (double) nmesh3;
#ifdef ONLY_ZA
		  P[n].Pos[axes] += dis;
		  P[n].Vel[axes] = dis * vel_prefac;
#else
		  P[n].Pos[axes] += dis - 3./7. * dis2;
		  P[n].Vel[axes] = dis * vel_prefac - 3./7. * dis2 * vel_prefac2;
#endif
		  
		  P[n].Pos[axes] = periodic_wrap(P[n].Pos[axes]);
		  
		  if(dis - 3./7. * dis2 > maxdisp)
		    maxdisp = dis;
		  if(dis2 > maxdisp2)
		    maxdisp2 = dis2;
                  meandisp2+=fabs(dis2);
                  meandisp+=fabs( dis - 3./7. * dis2);
		}
	    }
	}
      for(n = 0; n < NumPart; n++)
	{
          double vi = sqrt(P[n].Vel[0]*P[n].Vel[0]+P[n].Vel[1]*P[n].Vel[1]+P[n].Vel[2]*P[n].Vel[2]);
          meanvel  += vi;
          if(vi > maxvel) maxvel = vi;
        }
      
      if(ThisTask==0)printf("\n Reading out displacements\n");
      //get_time(time_x);
      
    }

  for(axes = 0; axes < 3; axes++) free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) free(cdisp2[axes]);
  
  gsl_rng_free(random_generator);
  
  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxdisp2, &max_disp2_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxvel, &max_vel_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&meanvel, &mean_vel_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&meandisp, &mean_disp_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&meandisp2, &mean_disp2_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if(ThisTask == 0)
    {
      printf("\nMaximum displacement: %g gadget units, mean displacement= %g\n",
	     max_disp_glob, mean_disp_glob / ((double)TotNumPart));
      printf("\nMaximum 2nd order displacement: %g gadget units, mean absolute value = %g\n",
	     max_disp2_glob, mean_disp2_glob / ((double)TotNumPart));
      printf("\nMaximum speed: %g gadget units, and mean speed = %g\n",
	     max_vel_glob, mean_vel_glob / ((double)TotNumPart));
    }
}


void initialize_ffts(void)
{
  int i, additional;
  int local_ny_after_transpose, local_y_start_after_transpose;
  int *slab_to_task_local;
  size_t bytes;
  
#ifdef USE_FFTW3
  const ptrdiff_t L=Nmesh, M=Nmesh, N=Nmesh;
  ptrdiff_t local_n0, local_n0_start, total_size;
#else
  int total_size;
#endif
  
  additional = (Nmesh) * (2 * (Nmesh / 2 + 1));    /* additional plane on the right side */  
  
#ifdef USE_FFTW3
  fftw_mpi_init();
  
  total_size = fftw_mpi_local_size_3d(L,M,N/2+1,MPI_COMM_WORLD,&local_n0,&local_n0_start);
  Local_nx=local_n0;
  Local_x_start=local_n0_start;
  
  TotalSizePlusAdditional = total_size + additional;
  
  bytes=sizeof(double)*total_size;  
#else
  Forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  
  Inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  
  rfftwnd_mpi_local_sizes(Forward_plan, &Local_nx, &Local_x_start,
			  &local_ny_after_transpose, &local_y_start_after_transpose, &total_size);
  
  Workspace = (fftw_real *) malloc(bytes = sizeof(fftw_real) * total_size);    
  
  ASSERT_ALLOC(Workspace);
  
  TotalSizePlusAdditional = total_size + additional;      
#endif
  
  // Set up slabs to contain particles...
  
  Local_nx_table = malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);
  
  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	fprintf(stdout,"Task=%d Local_nx=%d\n", i, Local_nx_table[i]);
      fflush(stdout);
    }
  
  Slab_to_task = malloc(sizeof(int) * Nmesh);
  slab_to_task_local = malloc(sizeof(int) * Nmesh);
  
  for(i = 0; i < Nmesh; i++)
    slab_to_task_local[i] = 0;
  
  for(i = 0; i < Local_nx; i++)
    slab_to_task_local[Local_x_start + i] = ThisTask;
  
  MPI_Allreduce(slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  free(slab_to_task_local);
}

