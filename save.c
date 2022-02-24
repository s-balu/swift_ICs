#include <math.h>
#include <stdlib.h>
#include "hdf5.h"
#include "mpi.h"

#include "allvars.h"
#include "proto.h"

void write_header_attributes_in_hdf5(hid_t handle);

#ifdef SWIFT_ICS
void write_units_attributes_in_hdf5(hid_t handle);
#endif

static double meanspacing, shift_gas, shift_dm, shift;
static int blockmaxlen, maxidlen, maxlongidlen;
static size_t bytes;

void write_properties_block(FILE *fd, int mode);
void write_properties_hdf5(hid_t handle[6], char *buf, int mode);

void write_particle_data(void)
{
  int nprocgroup, groupTask, masterTask;
  
  if(ThisTask == 0)
    printf("\nwriting initial conditions... \n");
  
  if((NTask < NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      FatalError(24131);
    }
  
  nprocgroup = NTask / NumFilesWrittenInParallel;
  
  if((NTask % NumFilesWrittenInParallel))
    nprocgroup++;
  
  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  MPI_Barrier(MPI_COMM_WORLD);  
#ifndef SINGLE_HDF5
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask)) // this processor's turn
#endif
	save_local_data();
      
      MPI_Barrier(MPI_COMM_WORLD); // wait inside the group
#ifndef SINGLE_HDF5
    }
#endif
//        MPI_Barrier(MPI_COMM_WORLD); // wait inside the group

  if(ThisTask == 0)
    printf("done with writing initial conditions.\n");
}

void save_local_data(void)
{
#define BUFFER 250

  float *block;
  int *blockid;
  long long *blocklongid;
  int4byte dummy;
  FILE *fd;
  char buf[300],filename[300];
  char hdf5_grpname[6][20]={"Coordinates","Velocities","ParticleIDs",
			    "Masses","InternalEnergy","SmoothingLength"};
  int i, k, pc;

  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
#ifdef SWIFT_ICS
  hid_t hdf5_unitsgrp = 0;
#endif
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank, pcsum = 0;
  int mpi_status;
  MPI_Info FILE_INFO_TEMPLATE;
#ifdef SINGLE_HDF5
  hid_t hdf5_file_ap_list;
  hid_t hdf5_data_cp_list;
#endif
  int type;
#ifdef NO64BITID
  int id_offset;
#else
  long long id_offset;
#endif

  /*
  if(NumPart == 0)
    return;
  */

#ifdef SINGLE_HDF5  
  MPI_Allreduce(&NumPart,&MaxNumPart,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if(NumPart<MaxNumPart)
    {
      fprintf(stdout,"MaxNumPart %d Task %d\n",MaxNumPart,ThisTask);
      fflush(stdout);
    }
#endif  
  
  // Define variables to be set in the header...
  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.n_all_high[i] = 0;
      header.mass[i] = 0;
    }
  
  header.npartTotal[1] = 1; //Put this before doing the check as npartTotal is obtained from the glass file (which we're not using) PA. Use a GlassTileFac = N_mesh instead as npartTotal is not a long int.
  
  size_t n_part_total[N_GADGET_TYPE];
  size_t totnumpart_check;
  for(i = 0,totnumpart_check=0; i < 3; i++){
    n_part_total[i] = (size_t)(header.npartTotal[i]) * (size_t)GlassTileFac * (size_t)GlassTileFac * (size_t)GlassTileFac;
    totnumpart_check+=n_part_total[i];
  }
  
  if(totnumpart_check!=TotNumPart){
    printf("totnumpart_check!=TotNumPart (ie. %zu!=%lld) on rank %d\n",totnumpart_check,TotNumPart,ThisTask);
  }
  
  for(i = 0; i < 3; i++){
    header.n_all_high[i]=(uint4byte)((n_part_total[i])>>32);
    header.n_all_high[i]=0;
#ifdef LONGIDS    
    header.npartTotal[i]=(n_part_total[i]);//-(((size_t)(header.n_all_high[i]))<<32));
#else
    header.npartTotal[i]=(uint4byte)(n_part_total[i]);//-(((size_t)(header.n_all_high[i]))<<32));
#endif
  }
  
  header.time = InitTime;
  header.redshift = 1.0 / InitTime - 1;
  
  header.flag_sfr = header.flag_feedback = header.flag_cooling =
    header.flag_stellarage = header.flag_metals = 0;
  
#ifndef SINGLE_HDF5
  header.num_files = NTaskWithN;
#else
  header.num_files = 1;
#endif
  
  header.BoxSize = Box;
  header.Omega0 = Omega;
  header.OmegaLambda = OmegaLambda;
  header.HubbleParam = HubbleParam;
  
  header.npart[1] = NumPart;
  header.npartTotal[1] = TotNumPart;
  
#ifndef PRODUCEGAS
  header.mass[1] = (Omega) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPart;
#else
  header.npart[0] = NumPart;
  header.npartTotal[0] = TotNumPart;
  header.mass[0] = (OmegaBaryon) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPart;
  header.mass[1] = (Omega - OmegaBaryon) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPart;
#endif
  
#ifdef PRODUCEGAS
  meanspacing = Box / pow(TotNumPart, 1.0 / 3);
  shift_gas = -0.5 * (Omega - OmegaBaryon) / (Omega) * meanspacing;
  shift_dm = +0.5 * OmegaBaryon / (Omega) * meanspacing;
#endif
  
  // Now start the process of writing data to file...
  bytes = BUFFER * 1024 * 1024;
  
  blockmaxlen = bytes / (3 * sizeof(float));
  maxidlen = bytes / (sizeof(int));
  maxlongidlen = bytes / (sizeof(long long));
  
  if(NTaskWithN > 1)
    if(OutputFormat<3) {
      sprintf(filename, "%s/%s.%d", OutputDir, FileBase, ThisTask);
    } else {
#ifndef SINGLE_HDF5
      sprintf(filename, "%s/%s.%d.hdf5", OutputDir, FileBase, ThisTask);
#else
      sprintf(filename, "%s/%s.hdf5", OutputDir, FileBase);
#endif
    }
  else
    if(OutputFormat<3) {
      sprintf(filename, "%s/%s", OutputDir, FileBase);
    } else {
      sprintf(filename, "%s/%s.hdf5", OutputDir, FileBase);
    }
  
  if(OutputFormat<3) {
    if(!(fd = fopen(filename, "w")))
      {
	printf("Error. Can't write in file '%s'\n", buf);
	FatalError(10);
      }
  } else {
    
#ifdef SINGLE_HDF5
    if(ThisTask==0)
#endif      
      if((hdf5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT))<0)
	{
	  fprintf(stderr,"Error while opening %s on %d\n",buf,ThisTask);
	  fflush(stderr);
	  exit(1245);
	}    
  }
  
  // Write header information
  if(OutputFormat<3) {
    dummy = sizeof(header);
    my_fwrite(&dummy, sizeof(dummy), 1, fd);
    my_fwrite(&header, sizeof(header), 1, fd);
    my_fwrite(&dummy, sizeof(dummy), 1, fd);
  } else {
    
#ifdef SINGLE_HDF5
    if(ThisTask==0)
      {
#endif	
	if((hdf5_headergrp = H5Gcreate1(hdf5_file, "/Header", 0))<0)
	  {
	    fprintf(stderr,"Error writing file header on %d\n",ThisTask);
	    fflush(stderr);
	    exit(1234);
	  }
	write_header_attributes_in_hdf5(hdf5_headergrp);
	H5Gclose(hdf5_headergrp);
	
	fprintf(stdout,"Written header information (Task %d)\n",ThisTask);
	fflush(stdout);
#ifdef SINGLE_HDF5
      }
#endif
    
#ifdef SWIFT_ICS
    if(ThisTask==0)
      {
	if((hdf5_unitsgrp = H5Gcreate1(hdf5_file, "/Units", 0))<0)
	  {
	    fprintf(stderr,"Error writing file units on %d\n",ThisTask);
	    fflush(stderr);
	    exit(1234);
	  }
	write_units_attributes_in_hdf5(hdf5_unitsgrp);
	H5Gclose(hdf5_unitsgrp);
	fprintf(stdout,"Written units information (Task %d)\n",ThisTask);
	fflush(stdout);
      }
#endif
  }
  
#ifdef SINGLE_HDF5
  if(ThisTask==0)
    H5Fclose(hdf5_file);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
#ifdef SINGLE_HDF5
  hdf5_file_ap_list = H5Pcreate(H5P_FILE_ACCESS);
  mpi_status = MPI_Info_create(&FILE_INFO_TEMPLATE);
  
  if((hdf5_status=H5Pset_fapl_mpio(hdf5_file_ap_list, MPI_COMM_WORLD, FILE_INFO_TEMPLATE))<0)
    {
      fprintf(stdout,"Error setting parallel I/O... [Task %d]\n",ThisTask);
      fflush(stdout);
      exit(1234);
    }
  
  //if((hdf5_status=H5Pset_alignment(hdf5_file_ap_list,1024,4096))<0)
  //  {
  //    fprintf(stdout,"Error setting HDF5 alignment... [Task %d]\n",ThisTask);
  //    fflush(stdout);
  //    exit(1234);
  //  }
  
  hdf5_file = H5Fopen(filename, H5F_ACC_RDWR, hdf5_file_ap_list);
  H5Pclose(hdf5_file_ap_list);
#endif
  
  if(OutputFormat>=3) {
    for(i = 0; i < 6; i++)
#ifdef SINGLE_HDF5
      if(header.npartTotal[i] > 0)
#else
	if(header.npart[i] > 0)
#endif	
	  {
	    sprintf(buf, "/PartType%d", i);
	    hdf5_grp[i] = H5Gcreate1(hdf5_file, buf, 0);
	  }
  }
  
  // Now write the various properties blocks
  // Write out positions  
  fprintf(stdout,"Writing position on Task %d\n",ThisTask);
  fflush(stdout);
  if(OutputFormat<3)
    write_properties_block(fd,1);
  else
    write_properties_hdf5(hdf5_grp,hdf5_grpname[0],1);
  
  // Write out velocities  
  fprintf(stdout,"Writing velocities on Task %d\n",ThisTask);
  fflush(stdout);
  if(OutputFormat<3)
    write_properties_block(fd,2);
  else
    write_properties_hdf5(hdf5_grp,hdf5_grpname[1],2);

  // Write out IDs  
  fprintf(stdout,"Writing particle IDs on Task %d\n",ThisTask);
  fflush(stdout);
  if(OutputFormat<3)
    write_properties_block(fd,3);
  else
    write_properties_hdf5(hdf5_grp,hdf5_grpname[2],3);

#if defined(WRITE_MASSES) || defined(SWIFT_ICS)
  // Write out masses
  fprintf(stdout,"Writing masses on Task %d\n",ThisTask);
  fflush(stdout);
  if(OutputFormat<3)
    write_properties_block(fd,4);
  else
    write_properties_hdf5(hdf5_grp,hdf5_grpname[3],4);
#endif

#ifdef PRODUCEGAS
  //write zero temperatures if needed
  
  fprintf(stdout,"Writing internal energies on Task %d\n",ThisTask);
  fflush(stdout);
  if(OutputFormat<3)
    write_properties_block(fd,5);
  else
    write_properties_hdf5(hdf5_grp,hdf5_grpname[4],5);    
#endif

#if defined(SWIFT_ICS) && defined(PRODUCEGAS)
  fprintf(stdout,"Writing smoothing lengths on Task %d\n",ThisTask);
  fflush(stdout);    
  write_properties_hdf5(hdf5_grp,hdf5_grpname[5],6);
#endif

  if(OutputFormat>=3)
    {
      for(type = 5; type >= 0; type--)
#ifdef SINGLE_HDF5
	if(header.npartTotal[type] > 0)
#else
        if(header.npart[type] > 0)
#endif	  
	  hdf5_status=H5Gclose(hdf5_grp[type]);
      hdf5_status=H5Fclose(hdf5_file);
    } else {
    fclose(fd);
  }
}

void write_header_attributes_in_hdf5(hid_t handle)
{
    hsize_t adim[1] = { 6 };
    hid_t hdf5_dataspace, hdf5_attribute;
    herr_t hdf5_status;

#ifdef SWIFT_ICS
    double BoxSize_SWIFT[3] = {header.BoxSize,header.BoxSize,header.BoxSize};
    hsize_t bdim[1] = { 3 };
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    hdf5_status=H5Sset_extent_simple(hdf5_dataspace, 1, bdim, NULL);
#else
    hdf5_dataspace = H5Screate(H5S_SCALAR);
#endif
    hdf5_attribute = H5Acreate1(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
#ifndef SWIFT_ICS
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
#else
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, BoxSize_SWIFT);
#endif
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

#ifdef SWIFT_ICS
    int dimension = {3};
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Dimension", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &dimension);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
#endif
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_Entropy_ICs", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_entropyICs);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_Cooling", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_cooling);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_Sfr", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_sfr);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_Feedback", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_feedback);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_Metals", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_metals);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_StellarAge", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_stellarage);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);


    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Flag_DoublePrecision", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_doubleprecision);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

    //#ifndef SWIFT_ICS
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    hdf5_status=H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate1(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    hdf5_status=H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
#ifndef SINGLE_HDF5
#ifndef LONGIDS
    hdf5_attribute = H5Acreate1(handle, "NumPart_ThisFile", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
#else 
    hdf5_attribute = H5Acreate1(handle, "NumPart_ThisFile", H5T_NATIVE_ULLONG, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_ULLONG, header.npart);
#endif
#else
#ifndef LONGIDS
    hdf5_attribute = H5Acreate1(handle, "NumPart_ThisFile", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);    
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
#else
    hdf5_attribute = H5Acreate1(handle, "NumPart_ThisFile", H5T_NATIVE_ULLONG, hdf5_dataspace, H5P_DEFAULT);    
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_ULLONG, header.npartTotal);
#endif
#endif
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    //#endif    
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    hdf5_status=H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
#ifndef LONGIDS
    hdf5_attribute = H5Acreate1(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
#else
    hdf5_attribute = H5Acreate1(handle, "NumPart_Total", H5T_NATIVE_ULLONG, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_ULLONG, header.npartTotal);
#endif
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SIMPLE);
    hdf5_status=H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
    hdf5_attribute = H5Acreate1(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.n_all_high);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    //#ifndef SWIFT_ICS
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    //#endif

}

#ifdef SWIFT_ICS
void write_units_attributes_in_hdf5(hid_t handle)
{
    hsize_t adim[1] = { 6 };
    hid_t hdf5_dataspace, hdf5_attribute;
    herr_t hdf5_status;
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Unit length in cgs (U_L)", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &UnitLength_in_cm);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Unit mass in cgs (U_M)", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &UnitMass_in_g);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Unit time in cgs (U_t)", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &UnitTime_in_s);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    double dummy=1.0;
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Unit current in cgs (U_I)", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);
    
    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate1(handle, "Unit temperature in cgs (U_T)", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    hdf5_status=H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy);
    hdf5_status=H5Aclose(hdf5_attribute);
    hdf5_status=H5Sclose(hdf5_dataspace);        
}
#endif

void write_properties_hdf5(hid_t handle[6], char *buf, int mode)
{
  float *block;
#ifdef NO64BITID
  int *blockid;
#else
  long long *blockid;
#endif
  long long i, k, pc;
  
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
#ifdef SWIFT_ICS
  hid_t hdf5_unitsgrp = 0;
#endif
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  unsigned long long dims[2], count[2], start[2];
  int rank;
  long long pcsum = 0;
  int mpi_status;
  MPI_Info FILE_INFO_TEMPLATE;
#ifdef SINGLE_HDF5
  hid_t hdf5_file_ap_list;
  hid_t hdf5_data_cp_list;
#endif
  int type;
#ifdef NO64BITID
  int id_offset;
#else
  long long id_offset;
#endif

  if(!(block = malloc(bytes)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double)bytes);
      FatalError(24);
    }

  switch(mode)
    {
    case 1:
    case 2:
    case 4:
    case 5:
    case 6:      
      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
      break;
    case 3:
#ifdef NO64BITID
      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
#else
      hdf5_datatype = H5Tcopy(H5T_NATIVE_LLONG);
#endif
      break;
    }

  switch(mode)
    {      
    case 1:
    case 2:
      for(type=0;type<6;type++)
	{
#ifndef SINGLE_HDF5
	  if(header.npart[type]==0) continue;
	  dims[0] = header.npart[type];
#else
	  if(header.npartTotal[type]==0) continue;
	  dims[0] = header.npartTotal[type];
#endif
	  dims[1] = 3;
	  rank = 2;
	  
	  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
	  hdf5_dataset = H5Dcreate1(handle[type], buf, hdf5_datatype, hdf5_dataspace_in_file,H5P_DEFAULT);
	  hdf5_status=H5Sclose(hdf5_dataspace_in_file);	  
	  
	  pcsum=0;
	  
#ifdef SINGLE_HDF5
	  pcsum=(unsigned long long)ThisTask*header.npart[type];
	  if(header.npart[type]<MaxNumPart)
	    if(header.npart[type]>0)	    
	      pcsum=(unsigned long long)ThisTask*MaxNumPart;
#endif
	  
	  if(mode==1) {
	    shift=shift_dm;
	    if(type==0) shift=shift_gas;
	  }
	  
	  for(i = 0, pc = 0; i < header.npart[type]; i++)
	    {
	      if(mode==1) {
		for(k = 0; k < 3; k++)
		  {
		    //block[3 * pc + k] = P[i].Pos[k];
		    block[3 * pc + k] = periodic_wrap(P[i].Pos[k] + shift);
		  }
	      } else {
		for(k = 0; k < 3; k++)
		  {
		    block[3 * pc + k] = P[i].Vel[k];
		  }
	      }

	      pc++;
	      
	      if(pc == blockmaxlen)
		{
		  start[0] = pcsum;
		  start[1] = 0;
		  
		  count[0] = pc;
		  count[1] = 3;

		  hdf5_dataspace_memory = H5Screate_simple(rank, count, NULL);

		  hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
		  
		  hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						    H5S_SELECT_SET,
						    start, NULL, count, NULL);
		  
#ifdef SINGLE_HDF5
		  hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPIO_COLLECTIVE
		  hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);		  
#endif
		  if(header.npart[type]>0)
		    {
#endif		  

		      hdf5_status =
		  	H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
#ifndef SINGLE_HDF5		  
				 hdf5_dataspace_in_file, H5P_DEFAULT, block);
#else		  
		                 hdf5_dataspace_in_file, hdf5_file_ap_list, block);
#endif	       
	      
#ifdef SINGLE_HDF5
		    }
	          hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif
	          hdf5_status = H5Sclose(hdf5_dataspace_in_file);	      
		  hdf5_status = H5Sclose(hdf5_dataspace_memory);

		  pcsum+=pc;
		  pc=0;
		}
	    }
	  
	  if(pc>0)
	    {
	      start[0] = pcsum;
	      start[1] = 0;
	      
	      count[0] = pc;
	      count[1] = 3;

	      hdf5_dataspace_memory = H5Screate_simple(rank, count, NULL);

	      hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
	      
	      hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						H5S_SELECT_SET,
						start, NULL, count, NULL);
	      
#ifdef SINGLE_HDF5
	      hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPIO_COLLECTIVE
	      hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);		  
#endif
		  if(header.npart[type]>0)
		    {
#endif	      

		  hdf5_status =
		    H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
#ifndef SINGLE_HDF5		  
			     hdf5_dataspace_in_file, H5P_DEFAULT, block);
#else		  
		             hdf5_dataspace_in_file, hdf5_file_ap_list, block);
#endif

#ifdef SINGLE_HDF5
	            }
	      hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif
	      hdf5_status = H5Sclose(hdf5_dataspace_in_file);	      
	      hdf5_status = H5Sclose(hdf5_dataspace_memory);

	    }

      hdf5_status = H5Dclose(hdf5_dataset);
      }
      break;
    case 3:
#ifdef NO64BITID
      blockid = (int *)block;
#else
      blockid = (long long *)block;
#endif      
      for(type=0;type<6;type++)
	{
#ifndef SINGLE_HDF5
	  if(header.npart[type]==0) continue;
	  dims[0] = header.npart[type];
#else
	  if(header.npartTotal[type]==0) continue;
	  dims[0] = header.npartTotal[type];
#endif
	  dims[1] = 1;
	  rank = 1;

	  hdf5_dataspace_in_file = H5Screate_simple(rank, dims , NULL);	  
	  hdf5_dataset = H5Dcreate(handle[type], buf, hdf5_datatype, hdf5_dataspace_in_file,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  hdf5_status=H5Sclose(hdf5_dataspace_in_file);
	  
	  pcsum=0;

#ifdef SINGLE_HDF5
	  pcsum=(unsigned long long)ThisTask*header.npart[type];
	  if(header.npart[type]<MaxNumPart)
	    if(header.npart[type]>0)	    
	      pcsum=(unsigned long long)ThisTask*MaxNumPart;
#endif	  
	  id_offset=0;
	  
#ifdef PRODUCEGAS
	  if(type==1) id_offset+=TotNumPart;
#endif
	  for(i = 0, pc = 0; i < header.npart[type]; i++)
	    {
	      blockid[pc] = P[i].ID+id_offset;
	      pc++;
	      
#ifdef NO64BITID	      
	      if(pc == maxidlen)
#else		
	      if(pc == maxlongidlen)
#endif		
		{
		  start[0] = pcsum;
		  start[1] = 0;
		  
		  count[0] = pc;
		  count[1] = 1;
		  
		  hdf5_dataspace_memory = H5Screate_simple(rank, count, NULL);
		  
		  hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
		  
		  hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						    H5S_SELECT_SET,
						    start, NULL, count, NULL);
		  
#ifdef SINGLE_HDF5
		  hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPIO_COLLECTIVE
		  hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);
#endif
		  if(header.npart[type]>0)
		    {
#endif		  

		  hdf5_status=
		    H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
#ifndef SINGLE_HDF5		  
			     hdf5_dataspace_in_file, H5P_DEFAULT, blockid);
#else		  
		             hdf5_dataspace_in_file, hdf5_file_ap_list, blockid);
#endif	       
	      
#ifdef SINGLE_HDF5
		    } 
	          hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif
		  hdf5_status = H5Sclose(hdf5_dataspace_in_file);
		  hdf5_status = H5Sclose(hdf5_dataspace_memory);
		  
		  pcsum += pc;
		  pc=0;
	        }
              }

            if(pc>0)
	      {
		start[0] = pcsum;
		start[1] = 0;
		
		count[0] = pc;
		count[1] = 1;
		
		hdf5_dataspace_memory = H5Screate_simple(rank, count, NULL);
		
		hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
		
		hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						  H5S_SELECT_SET,
						  start, NULL, count, NULL);
		
#ifdef SINGLE_HDF5
		hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPIO_COLLECTIVE
		hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);
#endif
		if(header.npart[type]>0)
		  {		
#endif		

		hdf5_status=
		  H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
#ifndef SINGLE_HDF5		  
			   hdf5_dataspace_in_file, H5P_DEFAULT, blockid);
#else		  
		           hdf5_dataspace_in_file, hdf5_file_ap_list, blockid);
#endif	       
		  
#ifdef SINGLE_HDF5
	           }
                hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif
                hdf5_status = H5Sclose(hdf5_dataspace_in_file);
	        hdf5_status = H5Sclose(hdf5_dataspace_memory);
               }
            hdf5_status = H5Dclose(hdf5_dataset);
            }
            break;
    case 4:
    case 5:
    case 6:
      for(type=0;type<6;type++)
        {
	  if((mode>4)&&(type>0)) break;
#ifndef SINGLE_HDF5
	  if(header.npart[type]==0) continue;
	  dims[0] = header.npart[type];
#else
	  if(header.npartTotal[type]==0) continue;
	  dims[0] = header.npartTotal[type];
#endif
	  dims[1] = 1;
	  rank = 1;

	  hdf5_dataspace_in_file = H5Screate_simple(rank, dims , NULL);	  
	  hdf5_dataset = H5Dcreate(handle[type], buf, hdf5_datatype, hdf5_dataspace_in_file,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	  hdf5_status=H5Sclose(hdf5_dataspace_in_file);
	  
	  pcsum=0;
	  
#ifdef SINGLE_HDF5
	  pcsum=(unsigned long long)ThisTask*header.npart[type];
	  if(header.npart[type]<MaxNumPart)
	    if(header.npart[type]>0)	    
	      pcsum=(unsigned long long)ThisTask*MaxNumPart;
#endif
	  
	  for(i = 0, pc = 0; i < header.npart[type]; i++)
            {
	      switch(mode)
		{
		case 4:
		  block[pc] = header.mass[type];
		  break;
		case 5:
		  block[pc] = 0.0;
		  break;
		case 6:
		  block[pc] = meanspacing;
		  break;
		}
	      
	      pc++;
              
	      if(pc == blockmaxlen)
		{
		  start[0] = pcsum;
		  start[1] = 0;
                  
		  count[0] = pc;
		  count[1] = 1;

		  hdf5_dataspace_memory = H5Screate_simple(rank, count, NULL);

		  hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
		  
		  hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						    H5S_SELECT_SET,
						    start, NULL, count, NULL);
#ifdef SINGLE_HDF5
		  hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPIO_COLLECTIVE
		  hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);
#endif
		  if(header.npart[type]>0)
		    {		  
#endif		  
		  hdf5_status=
		    H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
#ifndef SINGLE_HDF5		  
			     hdf5_dataspace_in_file, H5P_DEFAULT, block);
#else		  
			     hdf5_dataspace_in_file, hdf5_file_ap_list, block);
#endif	       

#ifdef SINGLE_HDF5
		    }
	          hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif		  
		  hdf5_status = H5Sclose(hdf5_dataspace_in_file);
		  hdf5_status = H5Sclose(hdf5_dataspace_memory);

		  pcsum += pc;
		  pc=0;
		}
	    }

	  if(pc>0)
            {
	      start[0] = pcsum;
	      start[1] = 0;
              
	      count[0] = pc;
	      count[1] = 1;
	      
	      hdf5_dataspace_memory = H5Screate_simple(rank, count, NULL);
	      
	      hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);	      
	      
	      hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						H5S_SELECT_SET,
						start, NULL, count, NULL);
#ifdef SINGLE_HDF5
	      hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
#ifdef MPIO_COLLECTIVE
	      hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);		  
#endif
	      if(header.npart[type]>0)
		{	      
#endif	      
	      
	      hdf5_status=
		H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
#ifndef SINGLE_HDF5		  
			 hdf5_dataspace_in_file, H5P_DEFAULT, block);
#else		  
	                 hdf5_dataspace_in_file, hdf5_file_ap_list, block);
#endif	       

#ifdef SINGLE_HDF5
	         }
              hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif
	      hdf5_status = H5Sclose(hdf5_dataspace_memory);
	      hdf5_status = H5Sclose(hdf5_dataspace_in_file);
	    }
      hdf5_status = H5Dclose(hdf5_dataset);
      }
      break;
    }

   hdf5_status = H5Tclose(hdf5_datatype);

   free(block);
}


// Write properties blocks for GADGET SnapFormat=1 output

void write_properties_block(FILE *fd, int mode)
{
  float *block;
#ifdef NO64BITID
  int *blockid;
#else
  long long *blockid;
#endif
  int4byte dummy;
  int i, k, pc;

  if(!(block = malloc(bytes)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double)bytes);
      FatalError(24);
    }

  switch(mode)
    {
    case 1:
      dummy = sizeof(float) * 3 * NumPart;
#ifdef  PRODUCEGAS
      dummy *= 2;
#endif      
      break;
    case 2:
      dummy = sizeof(float) * 3 * NumPart;
#ifdef  PRODUCEGAS
      dummy *= 2;
#endif      
      break;
    case 3:
#ifdef NO64BITID
      dummy = sizeof(int) * NumPart;
#else
      dummy = sizeof(long long) * NumPart;
#endif
#ifdef  PRODUCEGAS
      dummy *= 2;
#endif      
      break;
    case 4:
      dummy = sizeof(float) * NumPart;
#ifdef  PRODUCEGAS
      dummy *= 2;
#endif
      break;
    case 5:
      dummy = sizeof(float) * NumPart;
      break;
    }
  
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  switch(mode)
    {
    case 1:
      for(i = 0, pc = 0; i < NumPart; i++)
	{
	  for(k = 0; k < 3; k++)
	    {
	      block[3 * pc + k] = P[i].Pos[k];
#ifdef  PRODUCEGAS
	      block[3 * pc + k] = periodic_wrap(P[i].Pos[k] + shift_gas);
#endif
	    }
	  
	  pc++;
	  
	  if(pc == blockmaxlen)
	    {
	      my_fwrite(block, sizeof(float), 3 * pc, fd);
	      pc = 0;
	    }
	}
      
      if(pc > 0)
	my_fwrite(block, sizeof(float), 3 * pc, fd);

#ifdef  PRODUCEGAS
      for(i = 0, pc = 0; i < NumPart; i++)
	{
	  for(k = 0; k < 3; k++)
	    {
	      block[3 * pc + k] = periodic_wrap(P[i].Pos[k] + shift_dm);
	    }
	  
	  pc++;
	  
	  if(pc == blockmaxlen)
	    {
	      my_fwrite(block, sizeof(float), 3 * pc, fd);
	      pc = 0;
	    }
	}
      if(pc > 0)
	my_fwrite(block, sizeof(float), 3 * pc, fd);
#endif
      break;
    case 2:     // Velocities
      for(i = 0, pc = 0; i < NumPart; i++)
	{
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = P[i].Vel[k];
#ifdef MULTICOMPONENTGLASSFILE
	  if(WDM_On == 1 && WDM_Vtherm_On == 1 && P[i].Type == 1)
	    add_WDM_thermal_speeds(&block[3 * pc]);
#else
#ifndef PRODUCEGAS
	  if(WDM_On == 1 && WDM_Vtherm_On == 1)
	    add_WDM_thermal_speeds(&block[3 * pc]);
#endif
#endif
	  pc++;
	  
	  if(pc == blockmaxlen)
	    {
	      my_fwrite(block, sizeof(float), 3 * pc, fd);
	      pc = 0;
	    }
	}
      if(pc > 0)
	my_fwrite(block, sizeof(float), 3 * pc, fd);
#ifdef PRODUCEGAS
      for(i = 0, pc = 0; i < NumPart; i++)
	{
	  for(k = 0; k < 3; k++)
	    block[3 * pc + k] = P[i].Vel[k];
	  
	  if(WDM_On == 1 && WDM_Vtherm_On == 1)
	    add_WDM_thermal_speeds(&block[3 * pc]);
	  
	  pc++;
	  
	  if(pc == blockmaxlen)
	    {
	      my_fwrite(block, sizeof(float), 3 * pc, fd);
	      pc = 0;
	    }
	}
      if(pc > 0)
	my_fwrite(block, sizeof(float), 3 * pc, fd);
#endif
      break;
    case 3:     // IDs
#ifdef NO64BITID
      blockid = (int *)block;
#else
      blockid = (long long *)block;
#endif      
      for(i = 0, pc = 0; i < NumPart; i++)
	{
	  blockid[pc] = P[i].ID;

	  pc++;
#ifdef NO64BITID	  
	  if(pc == maxidlen)
#else
	  if(pc == maxlongidlen)
#endif
	    {
#ifdef NO64BITID
	      my_fwrite(blockid, sizeof(int), pc, fd);
#else
	      my_fwrite(blockid, sizeof(long long), pc, fd);
#endif
	      pc = 0;
	    }
	}
      
      if(pc > 0)
	{
#ifdef NO64BITID
	  my_fwrite(blockid, sizeof(int), pc, fd);
#else
	  my_fwrite(blockid, sizeof(long long), pc, fd);
#endif
	}
      
#ifdef PRODUCEGAS
      for(i = 0, pc = 0; i < NumPart; i++)
	{
	  blockid[pc] = P[i].ID + TotNumPart;

	  pc++;
#ifdef NO64BITID
	  if(pc == maxidlen)	  
#else
	  if(pc == maxlongidlen)
#endif
	    {
#ifdef NO64BITID
	      my_fwrite(blockid, sizeof(int), pc, fd);
#else
	      my_fwrite(blockid, sizeof(long long), pc, fd);
#endif
	      pc = 0;
	    }
	}
      if(pc > 0)
	{
#ifdef NO64BITID
	  my_fwrite(blockid, sizeof(int), pc, fd);
#else
	  my_fwrite(blockid, sizeof(long long), pc, fd);
#endif
	}
#endif
      break;
    case 4:     // Masses
      for(i = 0, pc = 0; i < NumPart; i++)
        {
#ifdef PRODUCEGAS
	  block[pc] = header.mass[0];
#else
	  block[pc] = header.mass[1];
#endif
	  pc++;
          
	  if(pc == blockmaxlen)
            {
	      my_fwrite(block, sizeof(float), pc, fd);
	      pc = 0;
            }
        }
      
      if(pc > 0)
	my_fwrite(block, sizeof(float), pc, fd);

#ifdef PRODUCEGAS
      for(i = 0, pc = 0; i < NumPart; i++)
        {
	  block[pc] = header.mass[1];
	  
	  pc++;
          
	  if(pc == blockmaxlen)
            {
	      my_fwrite(block, sizeof(float), pc, fd);
	      pc = 0;
            }
        }

      if(pc > 0)
	my_fwrite(block, sizeof(float), pc, fd);
#endif

      break;
    case 5:    // Internal energy
      for(i = 0, pc = 0; i < NumPart; i++)
	{
          block[pc] = 0;
          pc++;
	  
          if(pc == blockmaxlen)
	    {
              my_fwrite(block, sizeof(float), pc, fd);
              pc = 0;
	    }
	}
      if(pc > 0)
	my_fwrite(block, sizeof(float), pc, fd);
      break;      
    }

  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  free(block);
}

/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      FatalError(777);
    }
  return nwritten;
}

/* This catches I/O errors occuring for fread(). In this case we better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      FatalError(778);
    }
  return nread;
}
/*
void write_properties_hdf5(hid_t handle[6], char *buf, int mode)
{
  float *block;
#ifdef NO64BITID
  int *blockid;
#else
  long long *blockid;
#endif
  int i, k, pc;
  
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
#ifdef SWIFT_ICS
  hid_t hdf5_unitsgrp = 0;
#endif
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank, pcsum = 0;
  int mpi_status;
  MPI_Info FILE_INFO_TEMPLATE;
#ifdef SINGLE_HDF5
  hid_t hdf5_file_ap_list;
  hid_t hdf5_data_cp_list;
#endif
  int type;
#ifdef NO64BITID
  int id_offset;
#else
  long long id_offset;
#endif

  if(!(block = malloc(bytes)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double)bytes);
      FatalError(24);
    }
  
#ifdef SINGLE_HDF5
  hdf5_data_cp_list = H5Pcreate(H5P_DATASET_CREATE);
#endif
  
  switch(mode)
    {
    case 1:
    case 2:
    case 4:
    case 5:
    case 6:      
      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
      break;
    case 3:
#ifdef NO64BITID
      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
#else
      hdf5_datatype = H5Tcopy(H5T_NATIVE_LLONG);
#endif
      break;
    }

  switch(mode)
    {
    case 1:
    case 2:
    case 4:
    case 5:
    case 6:
      for(type=0;type<6;type++)
        {
	  if(mode>4 && type==0) continue;
#ifndef SINGLE_HDF5
	  if(header.npart[type]==0) continue;
	  dims[0] = header.npart[type];
#else
	  if(header.npartTotal[type]==0) continue;
	  dims[0] = header.npartTotal[type];
#endif
	  dims[1] = 1;
	  rank = 1;
	  
#ifdef SINGLE_HDF5
	  hdf5_status = H5Pset_chunk(hdf5_data_cp_list, rank, dims);
#endif
	  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
#ifndef SINGLE_HDF5
	  hdf5_dataset = H5Dcreate1(handle[type], buf, hdf5_datatype, hdf5_dataspace_in_file,H5P_DEFAULT);
#else
	  hdf5_dataset = H5Dcreate(handle[type], buf, hdf5_datatype, hdf5_dataspace_in_file,H5P_DEFAULT,hdf5_data_cp_list,H5P_DEFAULT);
#endif

#ifdef SINGLE_HDF5
	  hdf5_status = H5Pclose(hdf5_data_cp_list);
#endif
	  hdf5_status = H5Sclose(hdf5_dataspace_in_file);
	  
#ifdef SINGLE_HDF5
	  hdf5_file_ap_list = H5Pcreate(H5P_DATASET_XFER);
	  hdf5_status = H5Pset_dxpl_mpio(hdf5_file_ap_list, H5FD_MPIO_COLLECTIVE);
#endif
	  pcsum=0;
#ifdef SINGLE_HDF5
	  pcsum+=ThisTask*NumPart;
#endif
	  for(i = 0, pc = 0; i < NumPart; i++)
            {
	      switch(mode)
		{
		case 4:
		  block[pc] = header.mass[type];
		  break;
		case 5:
		  block[pc] = 0.0;
		  break;
		case 6:
		  block[pc] = meanspacing;
		  break;
		}
	      
	      pc++;
              
	      if(pc == blockmaxlen)
		{
		  start[0] = pcsum;
		  start[1] = 0;
                  
		  count[0] = pc;
		  count[1] = 1;
		  pcsum += pc;

		  hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
		  
		  hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,start, NULL, count, NULL);
                  
		  dims[0] = pc;
		  dims[1] = 1;
		  hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);
		  
		  hdf5_status=
#ifndef SINGLE_HDF5
		    H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
			     hdf5_dataspace_in_file, H5P_DEFAULT, block);
#else
		  H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
			   hdf5_dataspace_in_file, hdf5_file_ap_list, block);
#endif
		  hdf5_status = H5Sclose(hdf5_dataspace_memory);
		  hdf5_status = H5Sclose(hdf5_dataspace_in_file);
		  
		  pc=0;
		}
            }

	  if(pc>0)
            {
	      start[0] = pcsum;
	      start[1] = 0;
              
	      count[0] = pc;
	      count[1] = 1;
	      pcsum+=pc;

	      hdf5_dataspace_in_file = H5Dget_space(hdf5_dataset);
	      
	      hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
						H5S_SELECT_SET,
						start, NULL, count, NULL);
	      
	      dims[0] = pc;
	      dims[1] = 1;
	      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);
              
	      hdf5_status=
#ifndef SINGLE_HDF5
                H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
                         hdf5_dataspace_in_file, H5P_DEFAULT, block);
#else
	      H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
		       hdf5_dataspace_in_file, hdf5_file_ap_list, block);
#endif
	      
	      hdf5_status = H5Sclose(hdf5_dataspace_memory);
	      hdf5_status = H5Sclose(hdf5_dataspace_in_file);
            }

	  
	  hdf5_status = H5Dclose(hdf5_dataset);

#ifdef SINGLE_HDF5
	  hdf5_status = H5Pclose(hdf5_file_ap_list);
#endif
        }
      break;
      
    }
  
  hdf5_status = H5Tclose(hdf5_datatype);
    
  free(block);
}
*/
