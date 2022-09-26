/* DNSdataReader.c */
/* Wolfgang Tichy and Michal Pirog 8, 2022 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include <sys/stat.h>

/* include cactus stuff */
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

/* include my things */
#include "DNSdataReader.h"
#include "utilities.h"


#define NDATAMAX 23
#define STRLEN 16384
#define MAXPIECES 10

/* stuff for pwp */
int    DNS_EoS_n_pieces;
double DNS_EoS_rho0[MAXPIECES];
double DNS_EoS_kappa[MAXPIECES];
double DNS_EoS_n[MAXPIECES];
double DNS_EoS_k[MAXPIECES];
double DNS_EoS_q[MAXPIECES];
double DNS_EoS_P[MAXPIECES];

/* global vars that are set to sgrid pars */
double sgrid_x_CM; // sgrid's CM of system"
double Omega; // angular vel."
double ecc; // eccentricity"
double rdot; // radial vel."
double m01; // rest mass of star1"
double m02; // rest mass of star2"
double xmax1; // pos. of max density in star1"
double xmax2; // pos. of max density in star2"

/* global counters */
int level_l; /* counter that simulates level->l of bam */

/* position filepointer after the string str */
int DNS1_position_fileptr_after_str(FILE *in, char *str)
{
  char line[STRLEN];
  
  while(fgets(line, STRLEN-1, in)!=NULL)
  {
    if(strstr(line, str)!=NULL) return 1; //break;
  }

//  printf("DNS1_position_fileptr_after_str: done_01 \n");
  return EOF;
}


/* Compute DNSdataReader data */
void DNS1_dataReader(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  /* WT: I think cactus already has double *x, *y, *z; and other vars */
  /* for some reason there is a vel array and not vx, vy, vz vars: */
  CCTK_REAL *matter_vx = vel;
  CCTK_REAL *matter_vy = vel + npoints;
  CCTK_REAL *matter_vz = vel + npoints*2;
  /* DR: backup the original value of lapse and shift here */
  CCTK_REAL *alp_def, *betax_def, *betay_def, *betaz_def;
  int MPIrank, MPIsize;
  int pr=0;
  int i, j, n;
  int ndata = NDATAMAX;
  double xb, yb, zb;
  FILE *fp1, *fp2;
  char gridfile[STRLEN], call_interpolator[STRLEN], IDfile[STRLEN];
  char sgridoutdir[STRLEN], sgridoutdir_previous[STRLEN];
  char sgridcheckpoint_indir[STRLEN], IDfile_new[STRLEN];
  char sgridparfile[STRLEN];
  char *sgridargs;
  char *stringptr;
  double rdotor = rdot/(xmax1-xmax2); 
  double xmax;
  int ret;
  double s180 = (1 - 2*rotation180);
  int use_interpolator = use_Interpolator;
  char command[10000];
  int status = 0;

  /* which variables to set */
  int set_lapse = CCTK_EQUALS(initial_lapse, "DNSdata");
  int set_dtlapse = CCTK_EQUALS(initial_dtlapse, "DNSdata");
  int set_shift = CCTK_EQUALS(initial_shift, "DNSdata");
  int set_dtshift = CCTK_EQUALS(initial_dtshift, "DNSdata");
  int set_hydro = CCTK_EQUALS(initial_hydro, "DNSdata");

  /* backup lapse and shift if needed */
  {
    if(!set_lapse)
    {
      alp_def = malloc(npoints*sizeof(*alp_def));
      memcpy(alp_def, alp, npoints*sizeof(*alp_def));
      printf("DNS1_dataReader: Lapse has been backed up. \n");
    }

    if(!set_shift)
    {
      betax_def = malloc(npoints*sizeof(*betax_def));
      betay_def = malloc(npoints*sizeof(*betay_def));
      betaz_def = malloc(npoints*sizeof(*betaz_def));
      memcpy(betax_def, betax, npoints*sizeof(*betax_def));
      memcpy(betay_def, betay, npoints*sizeof(*betay_def));
      memcpy(betaz_def, betaz, npoints*sizeof(*betaz_def));
      printf("DNS1_dataReader: Shift has been backed up. \n");
    }
  }

  /* get refinement level number and MPI rank and size */
  level_l = ilogb(cctk_levfac[0]);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank); /* find MPI rank of this proc */
  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

  printf("DNS1_dataReader: <CCTK_INFO> was called for output \n");
  CCTK_INFO ("Setting up DNS initial data");

  /* say where we are */
  printf("DNS1_dataReader: level_l=%d  rank=%d \n", level_l, MPIrank);

  /* initialize file names */
  sprintf(gridfile, "%s/grid_level_%d_proc_%d.dat", out_dir, level_l, MPIrank);
  sprintf(sgridoutdir, "%s/sgrid_level_%d_proc_%d", out_dir, level_l, MPIrank);
  sprintf(sgridoutdir_previous, "%s/sgrid_level_%d_proc_%d_previous", out_dir, level_l, MPIrank);
  snprintf(sgridcheckpoint_indir, STRLEN-1, "%s", sgrid_datadir);
  stringptr = strrchr(sgrid_datadir, '/'); /* find last / */
  if(stringptr==NULL) /* no / found in DNSdataReader_sgrid_datadir */
    snprintf(sgridparfile, STRLEN-1, "%s.par", sgrid_datadir);
  else
    snprintf(sgridparfile, STRLEN-1, "%s%s", stringptr+1, ".par");
  
  /* dump the grid points covered by this processor on this level in gridfile */
  fp1 = fopen(gridfile, "wb");
  if(fp1==NULL) DNS2_errorexits("could not open %s", gridfile);
  fprintf(fp1, "%s", "# this pointsfile contains the (x+sgrid_x_CM, y, z) of the bam grid points\n");
  fprintf(fp1, "%s", "$BEGIN_data:\n");
  for(i=0; i<npoints; i++)
  {
    double xs;

    xb = x[i] * s180; /* multiply by -1 if 180 degree rotation */
    yb = y[i] * s180;
    zb = z[i];
    xs = xb + sgrid_x_CM; /* shift */
    fwrite(&xs, sizeof(double), 1, fp1);
    fwrite(&yb, sizeof(double), 1, fp1);
    fwrite(&zb, sizeof(double), 1, fp1);
  }
  fclose(fp1);
  fflush(stdout);

  /* decide if we call the interpolator, or if we just read data */
  if(!use_interpolator)
  {
    /* If you want to use previously generated ID use these two lines
       Set IDfiles_dir to the dir where you put the "ID*.dat"
       files generated in a previous run */

    sprintf(IDfile, "%s/ID_level_%d_proc_%d.dat", IDfiles_dir, level_l, MPIrank);
    printf("DNS1_dataReader: Looking for file %s\n", IDfile);

    /* call interpolator if file cannot be opened */
    fp2 = fopen(IDfile, "rb");
    if(fp2==NULL)
    {
      printf("DNS1_dataReader: Cannot open %s\n", IDfile);
      /* check if other IDfile exists */
      sprintf(IDfile, "%s/ID_level_%d_proc_%d.dat", out_dir, level_l, MPIrank);
      fp2 = fopen(IDfile, "rb");
      if(fp2==NULL)
      {
        printf("DNS1_dataReader: Cannot open %s\n", IDfile);
        printf("DNS1_dataReader: Will call interpolator.\n");
        use_interpolator = 1;
      }
      else
      {
        printf("DNS1_dataReader: %s exists, no need to call interpolator.\n", IDfile);
        use_interpolator = 0;
        //copy_IDfile_from_previous = 1; /* copy it to current outdir later */
        fclose(fp2);
      }
    }
    else
    {
      printf("DNS1_dataReader: %s exists, no need to call interpolator.\n", IDfile);
      use_interpolator = 0;
      fclose(fp2);
    }
  }
  else 
    use_interpolator = 1;

  /* check if we call the interpolator, or if we just read data */
  if(use_interpolator)
  {
    sprintf(IDfile, "%s/ID_level_%d_proc_%d.dat", out_dir, level_l, MPIrank);
    sprintf(IDfile_new, "%s%s", IDfile, "_new");

    /* info */
    printf("DNS1_dataReader: sgridoutdir = %s\n", sgridoutdir);
    printf("DNS1_dataReader: sgridcheckpoint_indir = %s\n", sgridcheckpoint_indir);

    DNS2_remove_dir(sgridoutdir);
    mkdir("mkdir", S_IRWXU | S_IRWXG);
    printf("DNS1_dataReader:datareader: New directory created \n");

    fflush(stdout);

    /* call the interpolator that will generate the ID */
    sprintf(call_interpolator, "%s %s/%s "
            "--modify-par:BNSdata_Interpolate_pointsfile=%s "
            "--modify-par:BNSdata_Interpolate_output=%s "
            "--modify-par:outdir=%s "
            "--modify-par:checkpoint_indir=%s",
            sgrid_exe, sgrid_datadir, sgridparfile,
            gridfile, IDfile_new, sgridoutdir, sgridcheckpoint_indir);
    if(Interpolate_verbose)
      strcat(call_interpolator, " --modify-par:BNSdata_Interpolate_verbose=yes");
    if(Interpolate_max_xyz_diff>0.0)
    {
      char str[STRLEN];
      sprintf(str, " --modify-par:BNSdata_Interpolate_max_xyz_diff=%g",
              Interpolate_max_xyz_diff);
      strcat(call_interpolator, str);
    }
    if(!Interpolate_make_finer_grid2)
      strcat(call_interpolator, " --modify-par:BNSdata_Interpolate_make_finer_grid2_forXYZguess=no");
    if(!keep_sgrid_output)
      strcat(call_interpolator, " > /dev/null");
    fflush(stdout);
    ret = DNS2_run(call_interpolator);
    printf("DNS1_dataReader: DNS2_run returned: %d\n", ret);
    fflush(stdout);

    /* move IDfile_new to IDfile */
    rename(IDfile_new, IDfile);

    /* clean up */

    DNS2_remove_dir(sgridoutdir_previous);
    if(!keep_sgrid_output)
    DNS2_remove_dir(sgridoutdir);
  }
  else
  {
    printf("DNS1_dataReader: Skipping interpolator, reading data from %s\n", IDfile);
  }

  /* info about filenames */
  printf("DNS1_dataReader: gridfile = %s\n", gridfile);
  printf("DNS1_dataReader: IDfile = %s\n", IDfile);
  fflush(stdout);

  /* read ADM variables from files generated by the interpolator */
  fp2 = fopen(IDfile, "rb");
  if(fp2==NULL) DNS2_errorexits("could not open %s", IDfile);
  ndata = 20;
  j=DNS1_position_fileptr_after_str(fp2, "$BEGIN_data:\n");
  if(j==EOF) DNS2_errorexits("could not find $BEGIN_data: in %s", IDfile);
  /* go over all points */
  for(i=0; i<npoints; i++)
  {
    double detg;
    double dataline[NDATAMAX];
    double q, VRx,VRy,VRz;

    /* read ndata doubles, i.e. one line of:
       alpha DNSdata_Bx DNSdata_By DNSdata_Bz gxx gxy gxz gyy gyz gzz
       Kxx Kxy Kxz Kyy Kyz Kzz q VRx VRy VRz */
    n=fread(dataline, sizeof(double), ndata, fp2);
    j=0;
    alp[i]=dataline[j];  j++;
    betax[i]=dataline[j];  j++;
    betay[i]=dataline[j];  j++;
    betaz[i]=dataline[j];  j++;
    gxx[i]=dataline[j];  j++;
    gxy[i]=dataline[j];  j++;
    gxz[i]=dataline[j];  j++;
    gyy[i]=dataline[j];  j++;
    gyz[i]=dataline[j];  j++;
    gzz[i]=dataline[j];  j++;
    kxx[i]=dataline[j];  j++;
    kxy[i]=dataline[j];  j++;
    kxz[i]=dataline[j];  j++;
    kyy[i]=dataline[j];  j++;
    kyz[i]=dataline[j];  j++;
    kzz[i]=dataline[j];  j++;
    q = dataline[j];  j++;
    VRx = dataline[j];  j++;
    VRy = dataline[j];  j++;
    VRz = dataline[j];  j++;
    if(j!=n) DNS2_errorexits("error reading dataline from %s", IDfile);

    /* transform some tensor components, if we have a 180 degree rotation */
    betax[i] *= s180;
    betay[i] *= s180;
    gxz[i] *= s180;
    gyz[i] *= s180;
    kxz[i] *= s180;
    kyz[i] *= s180;
    VRx *= s180;
    VRy *= s180;

    /* convert q, VR to rho, press, eps, matter_vx */
    if(set_hydro)
    {
      if(q>0.0)
      {
        double vIx,vIy,vIz;
        double xix,xiy,xiz;
        double DNS_n, DNS_kappa, DNS_k, hmk;

        DNS1_select_polytrope_n_kappa_k_of_hm1(q, &DNS_n, &DNS_kappa, &DNS_k);
        hmk = q + (1.0 - DNS_k);  /* hmk = hm1 + 1 - k , hm1 = q*/
        rho[i] = pow(hmk/(DNS_kappa*(DNS_n+1.0)), DNS_n);
        press[i]   = rho[i] * hmk/((DNS_n)+1.0);
        eps[i]= q - press[i]/rho[i];
        // printf("q %e r %e q %e eps %e n %e kappa %e k %e\n", q, rho[i],press[i],eps[i],DNS_n, DNS_kappa, DNS_k);

        xb = x[i] * s180; /* in case there is 180 deg. rot. */
        xmax = (xb>0)?xmax1:xmax2;
        /* construct KV xi from Omega, ecc, rdot, xmax1-xmax2 */
        xix = -Omega*y[i] + x[i]*rdotor; /* CM is at (0,0,0) in bam */
        xiy =  Omega*(x[i] - ecc*xmax) + y[i]*rdotor;
        xiz = z[i]*rdotor;
        /* vI^i = VR^i + xi^i */
        vIx = VRx + xix;
        vIy = VRy + xiy;
        vIz = VRz + xiz;
        /* Note: vI^i = u^i/u^0 in DNSdata,
                 while matter_v^i = u^i/(alpha u^0) + beta^i / alpha
           ==> matter_v^i = (vI^i + beta^i)/alpha                     */
        matter_vx[i] = (vIx + betax[i])/alp[i];
        matter_vy[i] = (vIy + betay[i])/alp[i];
        matter_vz[i] = (vIz + betaz[i])/alp[i];
      }
      else
      {
        rho[i] = press[i] = eps[i] = 0.0;
        matter_vx[i] = matter_vy[i] = matter_vz[i] = 0.0;
      }
    }

    /* print g_ij, K_ij, beta^i, alpha */
    if(pr)
    {
      printf("DNS1_dataReader: (x,y,z)=(%g,%g,%g)\n", x[i],y[i],z[i]);

      printf("DNS1_dataReader: alpha=%g beta=%g %g %g\n",
	alp[i], betax[i], betay[i], betaz[i]);  

      printf("DNS1_dataReader: g=%g %g %g %g %g %g \n",
	gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i]);
    
      printf("DNS1_dataReader: k=%g %g %g %g %g %g\n",
	kxx[i], kxy[i], kxz[i], kyy[i], kyz[i], kzz[i]);

      printf("DNS1_dataReader: rho=%g p=%g epsl=%g v=%g %g %g\n",
        rho[i], press[i], eps[i],
        matter_vx[i], matter_vy[i], matter_vz[i]);

      printf("DNS1_dataReader: \n");
    }
    /* check if data makes sense */
    detg=(2.*gxy[i]*gxz[i]*gyz[i] + gxx[i]*gyy[i]*gzz[i] -
             gzz[i]*gxy[i]*gxy[i] - gyy[i]*gxz[i]*gxz[i] -
             gxx[i]*gyz[i]*gyz[i]);

    if(detg<=0)
    {
      printf("DNS1_dataReader: det(g_ij)=%g at ccc=i=%d:  x=%g y=%g z=%g\n",
             detg, i, x[i], y[i], z[i]);
       //errorexit("found a point with det(g_ij)<=0.");
    }
  } /* end of all points loop */

  /* set time derivs of lapse and shift if desired */
  {
    if(set_dtlapse)
    {
      printf("DNS1_dataReader: <CCTK_INFO> was called for output \n");
      CCTK_INFO("Setting time derivatives of lapse");
      DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_PASS_CTOC, alp, dtalp, Omega);
    }

    if(set_dtshift)
    {
      printf("DNS1_dataReader: <CCTK_INFO> was called for output \n");
      CCTK_INFO("Setting time derivatives of shift");
      DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_PASS_CTOC, betax, dtbetax, Omega);
      DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_PASS_CTOC, betay, dtbetay, Omega);
      DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_PASS_CTOC, betaz, dtbetaz, Omega);
    }
  }

  /* restore original lapse and shift if desired */
  {
    if(!set_lapse) {
      memcpy(alp, alp_def, npoints*sizeof(*alp));
      free(alp_def);
    }

    if(!set_shift) {
      memcpy(betax, betax_def, npoints*sizeof(*betax));
      memcpy(betay, betay_def, npoints*sizeof(*betay));
      memcpy(betaz, betaz_def, npoints*sizeof(*betaz));
      free(betax_def);
      free(betay_def);
      free(betaz_def);
    }
  }

  /* increase level_l counter */
  level_l++;
    
//  printf("DNS1_dataReader: done_02 \n");
}

/* get DNSdata pars */
void DNS1_dataPars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS; 
  DECLARE_CCTK_PARAMETERS;
  FILE *fp1;
  char str[STRLEN];
  char strn[STRLEN], strrho0[STRLEN];
  char datadir[STRLEN];
  double ret;
  int i, j, count, start; 

  CCTK_INFO ("DNS1_dataPars: Reading pars for DNS initial data");

  /* get DNSdataReader_sgrid_datadir and remove any trailing / */
  snprintf(datadir, STRLEN-1, "%s", sgrid_datadir);
  j = strlen(datadir);
  if(datadir[j-1]=='/')
  { 
    datadir[j-1]=0;
    ret=CCTK_ParameterSet(sgrid_datadir, "DNSdata", datadir);
    if(ret!=0) printf("DNS1_dataPars: CCTK_ParameterSet returned non-zero exit code!\n");
  }
  strcat(datadir, "/BNSdata_properties.txt");

  /* open file */
  fp1 = fopen(datadir, "r");
  if(fp1==NULL) DNS2_errorexits("could not open %s", datadir);

  /* move fp1 to place where time = 0 is */
  j=DNS1_position_fileptr_after_str(fp1, "NS data properties (time = 0):\n");
  if(j==EOF) DNS2_errorexits("could not find (time = 0) in %s", datadir);

  /* initialize pwp, set everything to default */  
  for(i=0; i<MAXPIECES; i++) 
  {
    DNS_EoS_rho0[i]    = 0;
    DNS_EoS_kappa[i]   = 0;
    DNS_EoS_n[i]       = 0;
    DNS_EoS_k[0]       = 1.;
  }  

  /* get pars from file */
  fgotonext(fp1, "n_list"); 
  fscanline(fp1, strn);

  fgotonext(fp1, "rho0_list");
  fscanline(fp1, strrho0);

  fgetparameter(fp1, "kappa", str);
  DNS_EoS_kappa[0] = atof(str);

  /* process string with n_list */
  start = 2;   count = 0;
  while(sscanf(strn+start, "%s", str)==1)
  {
    start += strlen(str);
    if(strn[start]==' ') start++;
    DNS_EoS_n[count] = atof(str);
    count=count+1;
  }
  DNS_EoS_n_pieces = count;

  if(DNS_EoS_n_pieces == 1)
  {
    printf("DNS1_dataPars: Only 1 entry in n_list  ->  single polytrope\n");
    printf("DNS1_dataPars: rho0_list %s\n", strrho0);
  }
  else
  {
    printf("DNS1_dataPars: pwp with %d pieces\n", DNS_EoS_n_pieces);
    /* process string with rho0_list */
    start = 1;
    /* use count=1 here, because DNS_EoS_rho0[1] is now the first entry
       in rho0_list of DNSdata_properties.txt */
    count = 1; 
    while(sscanf(strrho0+start, "%s", str)==1)
    {
      start += strlen(str);
      if(strrho0[start]==' ') start++;
      DNS_EoS_rho0[count] = atof(str);
      count=count+1;
    }
  }
  /* compute the constant k, k = a+1 according to arXiv:0812.2163*/
  for(i=0;i<DNS_EoS_n_pieces;i++)
  {
    if(i > 0)
    {  
      DNS_EoS_kappa[i] = DNS_EoS_kappa[i-1]*pow(DNS_EoS_rho0[i],(DNS_EoS_n[i]-DNS_EoS_n[i-1])/
                     (DNS_EoS_n[i-1]*DNS_EoS_n[i]));

      DNS_EoS_k[i] = ((DNS_EoS_k[i-1]*DNS_EoS_rho0[i]+DNS_EoS_n[i-1]
               * DNS_EoS_kappa[i-1]*pow(DNS_EoS_rho0[i],1.+1./DNS_EoS_n[i-1]))/DNS_EoS_rho0[i])
               - DNS_EoS_kappa[i]*DNS_EoS_n[i]*pow(DNS_EoS_rho0[i],1./DNS_EoS_n[i]);
    }
    DNS_EoS_q[i]  = DNS_EoS_k[i] + (DNS_EoS_n[i]+1)*DNS_EoS_kappa[i]*pow(DNS_EoS_rho0[i],1./DNS_EoS_n[i])-1.;
    DNS_EoS_P[i]  = DNS_EoS_kappa[i]*pow(DNS_EoS_rho0[i],1.+1./(DNS_EoS_n[i]));
  
    /* print rho, kappa, n, k in log-file*/ 
  }
  /* copy last column and set large end-value*/  
  DNS_EoS_rho0[DNS_EoS_n_pieces]  = DNS_EoS_rho0[(DNS_EoS_n_pieces-1)] + 1.e10;
  DNS_EoS_q[DNS_EoS_n_pieces]     = DNS_EoS_q[(DNS_EoS_n_pieces-1)]    + 1.e10;
  DNS_EoS_P[DNS_EoS_n_pieces]     = DNS_EoS_P[(DNS_EoS_n_pieces-1)]    + 1.e10;
  DNS_EoS_kappa[DNS_EoS_n_pieces] = DNS_EoS_kappa[(DNS_EoS_n_pieces-1)];
  DNS_EoS_n[DNS_EoS_n_pieces]     = DNS_EoS_n[(DNS_EoS_n_pieces-1)];
  DNS_EoS_k[DNS_EoS_n_pieces]     = DNS_EoS_k[(DNS_EoS_n_pieces-1)];

  /* get rest of sgrid pars */
  fgetparameter(fp1, "x_CM", str);
  sgrid_x_CM = atof(str);
  fgetparameter(fp1, "Omega", str);
  Omega = atof(str);
  fgetparameter(fp1, "ecc", str);
  ecc = atof(str);
  fgetparameter(fp1, "rdot", str);
  rdot = atof(str);
  fgetparameter(fp1, "m01", str);
  m01 = atof(str);
  fgetparameter(fp1, "m02", str);
  m02 = atof(str);
  /* shift xmax1/2 such that CM is at 0 */
  fgetparameter(fp1, "xmax1", str);
  xmax1 = atof(str)-sgrid_x_CM;
  fgetparameter(fp1, "xmax2", str);
  xmax2 = atof(str)-sgrid_x_CM;

  /* close file */
  fclose(fp1);

  /* info */
  printf("DNS1_dataPars: assuming q:=h-1\n");
  printf("DNS1_dataPars: Use pwp with %d pieces \n", DNS_EoS_n_pieces);
  printf("DNS1_dataPars: rho0            q           P           kappa          n           k \n");
  for(i=0;i<=DNS_EoS_n_pieces;i++)
  {
    printf("DNS1_dataPars: %e %e %e %e %e %e \n",
    DNS_EoS_rho0[i], DNS_EoS_q[i], DNS_EoS_P[i], DNS_EoS_kappa[i], DNS_EoS_n[i], DNS_EoS_k[i]);
  }
  printf("DNS1_dataPars: Omega = %g\n", Omega);
  printf("DNS1_dataPars: ecc = %g\n", ecc);
  printf("DNS1_dataPars: rdot = %g\n", rdot);
  printf("DNS1_dataPars: m01 = %g\n", m01);
  printf("DNS1_dataPars: m02 = %g\n", m02);
  printf("DNS1_dataPars: xmax1 - sgrid_x_CM = %g\n", xmax1);
  printf("DNS1_dataPars: xmax2 - sgrid_x_CM = %g\n", xmax2);
  printf("DNS1_dataPars: sgrid_x_CM = %g\n", sgrid_x_CM);

  /* init level_l counter */
  level_l = 0;
  
//  printf("DNS1_dataPars: done_03 \n");
}

/* select pwp pars */
void DNS1_select_polytrope_n_kappa_k_of_hm1(double hm1,
                                           double *n, double *kappa, double *k)
{
  int m;
  /* find m such that q=hm1<DNS_EoS_q[m+1] */
  for(m=0; m< DNS_EoS_n_pieces-1; m++)
    if(hm1 < DNS_EoS_q[m+1]) break;
  /* Note, if EoS_n_pieces=1 we always get m=0 */
  *n     = DNS_EoS_n[m];
  *kappa = DNS_EoS_kappa[m];
  *k     = DNS_EoS_k[m];
//  printf("DNS1_select_polytrope_n_kappa_k_of_hm1: done_04 \n");
}


/* Calculate the time deriv of a grid var in the inertial frame from it's
   spatial derivs assuming a helical Killing vector.
   The way it's done below should work for scalars. */
void DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_ARGUMENTS,
                                              CCTK_REAL *var, CCTK_REAL *dtvar,
                                              CCTK_REAL Omega)
{
  DECLARE_CCTK_ARGUMENTS;
  int i;
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  CCTK_REAL *dxvar;
  CCTK_REAL *dyvar;

  /* alloc space for x and y deriv of var */
  dxvar = calloc(npoints, sizeof(CCTK_REAL));
  dyvar = calloc(npoints, sizeof(CCTK_REAL));

  /* get x and y deriv of var */
  Diff_gv(cctkGH, 0, var, dxvar, -1);
  Diff_gv(cctkGH, 1, var, dyvar, -1);

  for(i=0; i<npoints; ++i)
  {
    CCTK_REAL ephix = -y[i];
    CCTK_REAL ephiy = +x[i];
    CCTK_REAL dphi_var = ephix * dxvar[i] + ephiy * dyvar[i];
    dtvar[i] = -Omega * dphi_var;
  }
  /* free derivs */
  free(dyvar);
  free(dxvar);
    
//  printf("DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV: done_05 \n");
}
