/*@Start***********************************************************/
/* GESBG C source File
*
*    File Name:  basic.e
*    Developer:  karla miller
*
*  Id: ssfpepijf1.e,v 1.1 2010-08-19 15:12:14 jfnielse Exp 
*/
 
/*@Synopsis 
*
* basic EPIC sequence for teaching
*
*/     
/*@End*********************************************************/
  
  
/*
 * Loads all gradient waveforms from external files.
 * jfnielse@umich.edu 
 */


@inline epic.h

@global


#include "em_psd_ermes.in"
#include "math.h"
#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "filter.h"

#ifndef HW_IO
#include "stdio.h"
#else
#include "stdioLib.h"
#endif

#include <string.h>
#include "support_func.h"
#include "epic_error.h"
#include "epicfuns.h"

#include "grad_rf_pdssfp.globals.h"
#include "ssfpepi.h"

#include "jfn_globaldefs.h"  
#include "jfn_files.h"        /* functions for reading external readout waveform files */
#include "jfn_rf_files.h"        /* functions for reading external RF waveform files */

@inline Prescan.e PSglobal
int debugstate = 1;





















@ipgexport
/*********************************************************************
 *                        IPGEXPORT SECTION                          *
 *                                                                   *
 * Standard C variables of _any_ type common for both the Host and   *
 * IPG PSD processes. Declare here all the complex type, e.g.,       *
 * structures, arrays, files, etc.                                   *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
@inline Prescan.e PSipgexport
RF_PULSE_INFO rfpulseInfo[RF_FREE];
EXTERN_FILENAME directory;

/* variables and structures for external waveforms */
gradfilehdr roinfo, crusherinfo, gxpdinfo, gxinfo, gyblipinfo, gyinfo;        /* readout (x and y) */
char g_fname[80];
char rffilename[100];    /* external rf file name */
char readoutfilename[100];    /* external readout file name (prefix) */


























@cv
/*********************************************************************
 *                             CV SECTION                            *
 *                                                                   *
 * Standard C variables of _limited_ types common for both the Host  *
 * and IPG PSD processes. Declare here all the simple types, e.g,    *
 * int, float, and C structures containing the min and max values,   *
 * and ID description, etc.                                          *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the IPG sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/

@inline loadrheader.e rheadercv

/* gyromagnetic ratio in Hz/G */
/* #define GAMMA 4258 */

/* Constants for RF pulses designed in Matlab */
#define RAMPTIME    272         /* gradient ramp time (us) from 0 to max */

#define MAXWAIT     1024        /* max width of wait pulses; determines maximum amount of EPI echo-time shifting */
#define MINWAIT     4           /* minimum width of space between pulses */
#define MINPW       8           /* minimum pulse width */


/*********************************************/
/* CVs: overall psd control                  */

int   mode      = 0   with {0,,,INVIS, "Scan mode: GRE BOLD (1), SSFP (2,3), RF-SSFP (4,5), pdSSFP (6), manual (0)", };
int   isSSFP    = 1   with {0,1,,INVIS, "Run as SSFP sequence (1) or GRE (0)? ",      };
int   rfindex   = 0   with {0,,,INVIS, "Index into /usr/g/bin/jfnielse/ssfpepi/rffiles.txt",            };
int   readoutindex = 0   with {0,,,INVIS, "Index into /usr/g/bin/jfnielse/ssfpepi/readoutfiles.txt",    };
int   rf_cycle  = 180 with {0,360,,INVIS, "RF phase cycling, degrees",                };
int   scandur   = 120 with {0,360,,INVIS, "Scan duration (sec)",                };
int   greflip   = 40 with {0,90,,INVIS, "GRE BOLD flip angle (degrees)",                };
int   ssfpflip  = 60 with {0,90,,INVIS, "SSFP flip angle (degrees)",                };

int   myopzres  = 16  with {4,16,,INVIS, "# z partition encodes",                     };
float scaleYfov = 1.0 with {0,1.0,,INVIS, "Amplitude of gy phase-encodes (0 to 1.0)", };
float scaleZfov = 1.0 with {0,1.0,,INVIS, "Amplitude of gz phase-encodes (0 to 1.0)", };
float scalegssx = 0   with {-1.0,1.0,,INVIS, "Amplitude of gx linear B0 correction blip in dual RF pulse", };
float scalegssy = 0   with {-1.0,1.0,,INVIS, "Amplitude of gy linear B0 correction blip in dual RF pulse", };

float multgxpd  = 0.0 with {0,1.0,,INVIS, "Amplitude of partial dephaser (0.0-1.0)",   };

/* excitation pulse timing */
int start_rf1    = 256 with {0,,,, "Start location of RF excitation.",};
int sspwait_rf1  = 256 with {0,,,, "Time for SSP pulses after RF excitation.",};

int nviews;
int npasses;
/* remember to set opfphases manually on console */


/*********************************************/
/* CVs: sequence timing & related parameters */

/* sequence timing */
int ssi_time = 50;         /* ssi time is 4*ssi_time us */ 
int end_pass_time = 2ms;   /* time for endpass sequence */
int start_ro;              /* start time for readout gradients */

/* disdaqs */
int ps_ndisdaq = 50;       /* prescan: # dummy cycles before acquire data */
int ndisdaq = 20;          /* scan: # dummy cycles before acquire data */


/*********************************/
/* CVs: filter-related variables */

/* settings for 2DFT readout echo */
int filter_echo1;
int rcvr_echo1;
int prefill_echo1;

/* filter slots */
int num_filter_slots;
int scanslot = 0 with {0,7,4,VIS, "Scan filter slot number",};
int prescanslot = 0 with {0,7,4,VIS, "MPS2/APS2 filter slot number",};


/**********************/
/* CVs: RF excitation */

float rf_spoil_seed = 117 with {0,,0,VIS,"Phase increment (degrees) for RF spoiling."};
int   rf_spoil_seed_cnt = 0;
float rf_phase = 0.0;
float gssamp = 0;              /* slice-select amplitude on plateau during RF transmission (G/cm) */

/* duration of RF pulses */
int res_rf1;     /* number of points in excitation. */
int rf1_dur;     /* Total duration of excitation. */

/* slice select RF */
float a_rf1;   int ia_rf1;    int pw_rf1;
int off_rf1 = 0;
float alpha_rf1 = 0.46;  float thk_rf1;
float gscale_rf1 = 1.0;  float flip_rf1;

/* slice select gradient */
float a_gzrf1; int ia_gzrf1;
int pw_gzrf1a; int pw_gzrf1d; int pw_gzrf1;

/* dephaser for slice select gradient */
float a_grfdep; int ia_grfdep;
int pw_grfdepa; int pw_grfdep; int pw_grfdepd;

/* rephaser for slice select gradient */
float a_grfrep; int ia_grfrep;
int pw_grfrepa; int pw_grfrep; int pw_grfrepd;

/* B1 scaling */
float xmtaddScan;



/*********************/
/* CVs: 2DFT readout */

/* readout dephaser */
float a_gxdep; int ia_gxdep;
int pw_gxdepa; int pw_gxdep; int pw_gxdepd;

/* readout gradient */
float a_gxw; int ia_gxw;
int pw_gxwa; int pw_gxw; int pw_gxwd;

/* phase encode dephaser */
float a_gydep; float a_gydepa; float a_gydepb;
int ia_gydep;  int ia_gydepwa; int ia_gydepwb;
int pw_gydepa; int pw_gydepd;  int pw_gydep;
int pw_gydep_tot;
float a_gydep2; float a_gydepa2; float a_gydepb2;
int ia_gydep2;  int ia_gydepwa2; int ia_gydepwb2;

/* phase encode loop params */
int endview_iamp;      /* last instruction phase amp */
float endview_scale;   /* ratio of last instruction amp to maximum value */
int gyFlag = 1;


/****************************/
/* CVs: miscellaneous stuff */

/* timing (NECESSARY??) */
int pos_start = 0 with {0,,,INVIS, "Start time for sequence. ",};
int tlead = 25us with {0,,25us,INVIS, "Init deadtime",};

int obl_debug = 0 with {0,1,0,INVIS,
    "On(=1) to print messages for obloptimize",};
int obl_method = 0 with {0,1,0,INVIS,
    "On(=1) to optimize the targets based on actual rotation matrices",};

float yfov_aspect = 1.0 with {0,,,INVIS, "acquired Y FOV aspect ratio to X",};

@inline Prescan.e PScvs

















































@host
/************************************************************************
 *                                                                      *
 *                                                                      *
 *				HOST SECTION 				*
 *                                                                      *
 * Write here the code unique to the Host PSD process. The following    *
 * functions must be declared here: cvinit(), cveval(), cvcheck(),      *
 * and predownload().                                                   *
 *                                                                      *
 *                                                                      *
 ************************************************************************/

#include "sar_pm.h"
#include "grad_rf_pdssfp.h"
#include "sokPortable.h"

FILTER_INFO *echo1_filt;
	/* These will point to a structure defining parameters of 
	the filter used for the 1st echo and 2nd through N echos */

FILTER_INFO echo1_rtfilt, echo1_rtfilt_frac;
	/* for V6 use real time  filters, so allocate space for 
	them instead of trying to point to an infinite number of 
	structures in filter.h. */

@inline Prescan.e PShostVars       /* added with new filter calcs */

/** Load PSD Header **/
abstract("Balanced SSFP sequence for fMRI with 3D stack-of-EPI readout");
psdname("ssfpepi2");

int num_conc_grad = 3;          /* always three for basic */

static char supfailfmt[] = "Support routine %s failed";

#define MAX_ENTRY_POINTS 9

/* peak B1 amplitudes */
float maxB1[MAX_ENTRY_POINTS], maxB1Seq;
int entry;



/************************************************************************/
/*       		  My_orderslice    	   	  		*/
/*  Does what orderslice would do for multi-phase sequences, but	*/
/*  does not assume a separate pass for each slice location. The	*/
/*  order is non-sequential, and we just assume that scan() will	*/
/*  put the data in the right place for each slice and phase.		*/
/*  Most of this code is copied from orderslice.c  			*/
/************************************************************************/

void my_orderslice(void) {

  int i, j, slice;

/*
  for (i = 0; i < opslquant*opfphases; i++) {
*/
  for (i = 0; i < opslquant; i++) {
        slice = i % opslquant;
	
        rsp_info[slice].rsptloc = scan_info[slice].optloc;
        rsp_info[slice].rsprloc = scan_info[slice].oprloc;
        rsp_info[slice].rspphasoff = scan_info[slice].opphasoff;
        rsptrigger[slice] = TRIG_INTERN;
	
        data_acq_order[i].slloc = i/opfphases;
        data_acq_order[i].slpass = 0;
        data_acq_order[i].sltime = i;

        for (j = 0; j <= 8; j++) {
            rsprot[i][j] = (INT)scan_info[slice].oprot[j]; 
	}

  }

  return;
}



/************************************************************************/
/*       			CVINIT    				*/
/* Invoked once (& only once) when the PSD host process	is started up.	*/
/* Code which is independent of any OPIO button operation is put here.	*/
/************************************************************************/
STATUS cvinit(void)
{

  fprintf(stderr, "\nwelcome to cvinit!\n");


  /*************************************************/
  /* CVINIT: initialize system                     */
  /*************************************************/


#ifdef ERMES_DEBUG
  use_ermes = 0;
#else
  use_ermes = 1;
#endif
  
  /* Initialize internal configuration variables */	
  EpicConf();

  /* Initialize physical and logical gradient structures, and
     physical limits on the system. The physical gradient 
     structure defines the physical quantities of the system
     itself - ie gradient strengths etc.  The logical gradient
     structure has z=slice select, x=readout, y=phase encode. */
  inittargets(&loggrd, &phygrd);

  /* Initialize receive filter */
  initfilter();

  /* Calculate limits/targets based on config control variables	*/
  if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant),
		  exist(opplane), exist(opcoax), obl_method, obl_debug,
		  &opnewgeo, cfsrmode)==FAILURE)
    return FAILURE;
  
  TARDIS_FREQ_OFFSET=RCV_FREQ_CERD;

  /* set some cv limits (MOVE TO PREDOWNLOAD??) */
  cvmax(rhfrsize,15000);
  cvmax(opnecho,32);
  cvmin(opxres,64);
  cvmin(opyres,16);


  /*************************************************/
  /* CVINIT: set prescription options              */
  /*************************************************/

  /* turn off some options (TE2 & BW) */
  pite2nub=0;
  pircb2nub = 0;
  pircbnub = 0;

  /* flip angle */
  pifanub = 5;         
  pifaval2 = 30;
  pifaval3 = 45;
  pifaval4 = 60;
  pifaval5 = 90;

  /* x res */
  pixresnub = use4-1;  
  pixresval2 = 128;
  pixresval3 = 192;
  pixresval4 = 256;

  /* y res */
  piyresnub  = use4-1; 
  piyresval2 = 128;
  piyresval3 = 192;
  piyresval4 = 256;

  /* TR */
  pitrnub = 5;         
  pitrval2 = 50ms;
  pitrval3 = 100ms;
  pitrval4 = 200ms;
  pitrval5 = 400ms;

  /* TE */
  pite1nub=5;
  pite1val2 = 25ms;
  pite1val3 = 50ms;
  pite1val4 = 100ms;
  pite1val5 = 200ms;

  /* FOV */
  pifovnub = 4;        
  pifovval2 = 160;
  pifovval3 = 200;
  pifovval4 = 240;

  /* slice thickness */
  pistnub  = 4;        
  pistval2 = 3.0;
  pistval3 = 5.0;
  pistval4 = 10.0;

  cvmod( opuser3, 0, 10, 3, 
     "Z-gradient crusher: number of phase cycles", 
     EM_PSD_SUPPORT_FAILURE, "Not possible.");
  cvmod( opuser16, 0, 1000000, 0, 
     "extra TR - us", 
     EM_PSD_SUPPORT_FAILURE, "Not possible.");

  if (setsysparms() == FAILURE) {
    epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE, 1,
               STRING_ARG, "setsysparms");
    return FAILURE;
  }


@inline Prescan.e PScvinit




  /*************************************************/
  /* CVINIT: initialize RF & echo params           */
  /*         MOVE TO PREDOWNLOAD??                 */
  /*************************************************/

  /* initialize echo pulse */
  filter_echo1 = scanslot;
  rcvr_echo1 = 0;
  prefill_echo1 = DEFPREFILLS;

#include "cvinit.in"	/* Runs the code generated by macros in preproc.*/

  return SUCCESS;
}




/************************************************************************/
/*       			CVEVAL    				                            */
/* Called w/ every OPIO button push which has a corresponding CV. 	    */
/* CVEVAL should only contain code which impacts the advisory panel -   */
/* put other code in cvinit or predownload				                */
/************************************************************************/
STATUS cveval(void)
{
  int temp;
  rfstruct rfinfo;   /* contains rf header info, and the waveforms themselves */
  char tmpfname[100];

@inline Prescan.e PScveval
  fprintf(stderr, "\nwelcome to cveval!\n");

  /* copied from Chuck's flprofile.e */
  #ifdef SIM
    strcpy(directory, "");
	 ndisdaq = 0;
  #else
    strcpy(directory, "/usr/g/bin/jfnielse/ssfpepi/");
  #endif 


  /********************************************************************/
  /* CVEVAL: set scan parameters for scan mode chosen by operator     */
  /********************************************************************/

  temp = _opflip.fixedflag; 
  _opflip.fixedflag = 0; 
  switch (mode)
  {
    case 1:                  /* GRE BOLD fMRI */
      isSSFP = 0;
      multgxpd = 0;
	  readoutindex = 0;
      rfindex = 0;
      opflip = greflip;
      break;
    case 2:                  /* SSFP fMRI */
      isSSFP = 1;
      multgxpd = 0;
	  readoutindex = 1;
      rfindex  = 0;
      rf_cycle = 180;
      opflip = ssfpflip;
      break;
    case 3:                  /* SSFP fMRI */
      isSSFP = 1;
      multgxpd = 0;
	  readoutindex = 1;
      rfindex  = 0;
      rf_cycle = 0;
      opflip = ssfpflip;
      break;
    case 4:                  /* RF-SSFP fMRI */
      isSSFP = 1;
      multgxpd = 0;
	  readoutindex = 1;
      rfindex  = 1;
      rf_cycle = 180;
      opflip = ssfpflip;
      break;
    case 5:                  /* RF-SSFP fMRI */
      isSSFP = 1;
      multgxpd = 0;
	  readoutindex = 1;
      rfindex  = 1;
      rf_cycle = 0;
      opflip = ssfpflip;
      break;
    case 6:                  /* pdSSFP fMRI */
      isSSFP = 1;
      multgxpd = 1;
	  readoutindex = 1;
      rfindex  = 0;
      rf_cycle = 180;
      opflip = ssfpflip;
      break;
  }
  _opflip.fixedflag = temp; 

  /********************************************************************/
  /* CVEVAL: read rf file header, and fill in rfpulse struct          */
  /********************************************************************/

  fprintf(stderr,"reading rf file header\n");
  jfn_rf_getfilename(tmpfname, rfindex, directory);   /* filenames are listed in "rffiles.txt" */
  sprintf(rffilename, "%s%s", directory, tmpfname);
  jfn_rf_readheader(rffilename, &rfinfo);

  pw_rf1  = 4*rfinfo.res;
  rf1_dur = 4*rfinfo.res;
  res_rf1 = rfinfo.res;

  rfpulse[RF1_SLOT].abswidth      = rfinfo.paramsfloat[1];
  rfpulse[RF1_SLOT].effwidth      = rfinfo.paramsfloat[2];
  rfpulse[RF1_SLOT].area          = rfinfo.paramsfloat[3];
  rfpulse[RF1_SLOT].dtycyc        = rfinfo.paramsfloat[4];
  rfpulse[RF1_SLOT].maxpw         = rfinfo.paramsfloat[5];
  rfpulse[RF1_SLOT].num           = rfinfo.paramsfloat[6];
  rfpulse[RF1_SLOT].max_b1        = rfinfo.paramsfloat[7];
  rfpulse[RF1_SLOT].max_int_b1_sq = rfinfo.paramsfloat[8];
  rfpulse[RF1_SLOT].max_rms_b1    = rfinfo.paramsfloat[9];
  rfpulse[RF1_SLOT].nom_fa        = rfinfo.paramsfloat[10];
  rfpulse[RF1_SLOT].nom_pw        = rfinfo.paramsfloat[11];
  rfpulse[RF1_SLOT].nom_bw        = rfinfo.paramsfloat[12];
  gssamp                          = rfinfo.paramsfloat[13];


  /* initialize RF pulse */
  flip_rf1 = opflip;
  a_rf1 = flip_rf1/180;
  pw_gzrf1 = rf1_dur;
  pw_rf1 = rf1_dur;
  thk_rf1 = opslthick;


  /********************************************************************/
  /* CVEVAL: read readout filename prefix                             */
  /********************************************************************/
  jfn_getfilename(tmpfname, readoutindex, directory); /* readout filename prefixes are listed in "readoutfiles.txt" */
  sprintf(readoutfilename, "%s%s", directory, tmpfname);
  fprintf(stderr,"\n\treadoutfilename = %s\n", readoutfilename);


  return SUCCESS;
}




/************************************************************************/
/*       			CVCHECK    				*/
/* Executed on each 'next page' to ensure prescription can proceed 	*/
/* to the next page. Redundant with cveval() at this point.             */
/************************************************************************/
STATUS cvcheck(void)
{
  fprintf(stderr,"\nwelcome to cvcheck!\n");
  return SUCCESS;
}




/************************************************************************/
/*             		    PREDOWNLOAD           		        */
/* Executed prior to a download - all operations not needed for the 	*/
/* advisory panel results.  Execute the	pulsegen macro expansions for	*/
/* the predownload section here.  All internal amps, slice ordering,  	*/
/* prescan slice calc., and SAT placement calculations are performed 	*/
/* in this section.  Time anchor settings for pulsegen are done in this */
/* section too.  				 			*/
/************************************************************************/
STATUS predownload(void)
{
  char filename[100];
  int i, temp;
  int rtime, ftime;
  int mintr, minte;
  float arearf1;
  float areagxw;
  float areagxdep;
  float areagydep;
  float areamult;
  float target;
  float ratiotemp;


  /* earliest possible time for rf pulse to play out.
     determined by addrfbits() pulse, which plays out 264 us before rf pulse,
     and by seqcoresync and seqcorescp ssp packets, which together take 9 us. 
     See EPIC manual and/or MREL twiki site for more info  -- jfn */
  /* note: use fastAddrfbits instead. Allows earlier rf start times -- jfn */
  int rf_ssp_wait = 120;   
                         

  fprintf(stderr,"\nwelcome to predownload!\n");


  /**************************************************************/
  /* PREDOWNLOAD:  read external waveform file headers          */
  /**************************************************************/

  sprintf(g_fname, "%s_gxpd.jfnwav", readoutfilename);
  if (jfn_readextgradhdr(g_fname, &gxpdinfo) == JFN_FAILURE)
    return(FAILURE);
  sprintf(g_fname, "%s_gx.jfnwav", readoutfilename);
  if (jfn_readextgradhdr(g_fname, &roinfo) == JFN_FAILURE)
    return(FAILURE);
  if (jfn_readextgradhdr(g_fname, &gxinfo) == JFN_FAILURE)
    return(FAILURE);
  sprintf(g_fname, "%s_gyblip.jfnwav", readoutfilename);
  if (jfn_readextgradhdr(g_fname, &gyblipinfo) == JFN_FAILURE)
    return(FAILURE);
  sprintf(g_fname, "%s_gy.jfnwav", readoutfilename);
  if (jfn_readextgradhdr(g_fname, &gyinfo) == JFN_FAILURE)
    return(FAILURE);
  sprintf(g_fname, "%scrusher.jfnwav", directory);
  if (jfn_readextgradhdr(g_fname, &crusherinfo) == JFN_FAILURE)
    return(FAILURE);

  nviews = (int) ROUND((float)gxinfo.npix/(float)gxinfo.nechoes);

  fprintf(stderr,"roinfo.res= %d\n", roinfo.res);
  fprintf(stderr,"roinfo.res_k    = %d\n", roinfo.res_k);
  fprintf(stderr,"roinfo.res_kpre = %d\n", roinfo.res_kpre);
  fprintf(stderr,"roinfo.npix = %d\n", roinfo.npix);
  fprintf(stderr,"roinfo.nechoes = %d\n", roinfo.nechoes);
  fprintf(stderr,"roinfo.nechosep = %d\n", roinfo.nechosep);
  fprintf(stderr,"nviews = %d\n", nviews);
  fprintf(stderr,"opyres       = %d\n", opyres);



  /*************************************************/
  /* PREDOWNLOAD: initialize params                */
  /*************************************************/

  /* delay associated with gradients (time from when ask for gradient &
    when it actually turns on). */
  /* psd_grd_wait = get_grad_dly(); */

  /* number of slices */
  pislquant = opslquant;
  slquant1 = 1;

  /* force single nex */
  _opnex.fixedflag = 0;
  /* opnex = 1; nex = 1;  exnex = 1; */
  nex = opnex;  exnex = opnex; 
  acq_type = TYPGRAD;

  /* echoes */
  temp = _opnecho.fixedflag; 
  _opnecho.fixedflag = 0;
  opnecho = 1 ;
  _opnecho.fixedflag = temp;

  /* total number of acquisitions (Pfiles) */
  rhnpasses = 1;
  npasses = 1;
  /* acqs = rhnpasses; */

  /* aspect ratio of FOV along y relative to prescribed FOV */
  yfov_aspect = nop*exist(opphasefov);


  /***************************************************************/
  /* PREDOWNLOAD: frame acquisition parameters                   */
  /***************************************************************/

  /* frame sizes */
  /* rhnframes = opyres*fn*yfov_aspect; */   /* total number of frames */
  rhnframes = nviews;                   /* number of frames in each EPI disc */
  /* rhfrsize = opxres;  */             /* length of each frame */
  rhfrsize = gxinfo.res_k;              /* length of each frame */
	
  rhexecctrl = 9;

  /* receive bandwidth fixed at scanner sampling rate */
  temp = _oprbw.fixedflag;
  _oprbw.fixedflag = 0;
  ratiotemp = (float)rint(125./oprbw);
  oprbw = 125./ratiotemp; 
  oprbw = 125.; 
  _oprbw.fixedflag = temp; 

  if (calcfilter( &echo1_rtfilt,    /* I:   all the filter parameters */
                  exist(oprbw),     /* I/O: desired and final allowable bw */
                  roinfo.res_k,          /* I:   output pts generated by filter */
                  OVERWRITE_NONE )    
	  	  /* oprbw will be updated in the next calcfilter call */
          == FAILURE) {
    epic_error(use_ermes,"%s failed",
         EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"calcfilter:echo1:full");
    return FAILURE;
  }

  echo1_filt = &echo1_rtfilt;

  /* Divide by 0 protection */
  if ((echo1_filt->tdaq == 0) || (echo1_filt->decimation == 0)) {
      epic_error(use_ermes,"echo1 tdaq or decimation = 0",EM_PSD_BAD_FILTER,0);
      return FAILURE;
  }

#include "sokPortable.h"

#include "predownload.in"  /* include 'canned' predownload code */



  /********************************************************/
  /* PREDOWNLOAD: calculate slice select amplitudes       */
  /********************************************************/

  /* RF pulse */
  ia_rf1 = a_rf1 * MAX_PG_IAMP;
  if (res_rf1 == 0) res_rf1 = pw_rf1/(RF_UPDATE_TIME);

  /* Z GRADIENT */
  if (ampslice(&a_gzrf1, (long)(rfpulse[RF1_SLOT].nom_bw),
	       (float)thk_rf1,(float)gscale_rf1,1) == FAILURE) 
     return FAILURE;
  fprintf(stderr, "a_gzrf1 = %f \n", a_gzrf1);

  ia_gzrf1 = a_gzrf1 * MAX_PG_IAMP/loggrd.tz;
  fprintf(stderr, "loggrd.tz = %f \n", loggrd.tz);
  fprintf(stderr, "ia_gzrf1 = %d \n", ia_gzrf1);

  if (optramp(&pw_gzrf1a,a_gzrf1,loggrd.tz,loggrd.zrt,1) == FAILURE)
     return FAILURE;
  if (optramp(&pw_gzrf1d,a_gzrf1,loggrd.tz,loggrd.zft,1)== FAILURE) 
      return FAILURE;

  arearf1 = a_gzrf1*((float)pw_gzrf1 + (float)pw_gzrf1a);

  /* slice-select dephaser/crusher */
  gettarget(&target, ZGRAD, &loggrd);
  getramptime(&rtime, &ftime, ZGRAD, &loggrd);
  rtime = RAMPTIME;
  ftime = RAMPTIME;

  if (isSSFP)  {
    areamult = -0.5;    /* balance Gz 0th moment */
  }
  else {   /* design crusher */
    areamult = (JFN_MAX(1,opuser3)*2*M_PI*1.0e4 / (JFN_GAMMA * 2*M_PI *
               opslthick)) / arearf1;
  }

  /* dephaser prior to RF pulse, or crusher (at end of TR) */
  if ( amppwgrad( areamult*arearf1, target/sqrt(3), 0.0, 0.0, rtime,
	          MIN_PLATEAU_TIME, &a_grfdep,
	          &pw_grfdepa, &pw_grfdep, &pw_grfdepd) == FAILURE )
     return FAILURE;
  ia_grfdep = (a_grfdep / target) * MAX_PG_IAMP;
   
  /* design rephaser following RF pulse */
  if ( amppwgrad( -0.5*arearf1, target/sqrt(3), 0.0, 0.0, rtime,
	          MIN_PLATEAU_TIME, &a_grfrep,
	          &pw_grfrepa, &pw_grfrep, &pw_grfrepd) == FAILURE )
     return FAILURE;
  ia_grfrep = (a_grfrep / target) * MAX_PG_IAMP;



  
  /*****************************************************************/
  /* PREDOWNLOAD: scan & prescan filters                           */
  /*****************************************************************/

  num_filter_slots = 0;
  psd_board_type = PSDCERD;

  /* Set up scan filter */
  setfilter(echo1_filt, SCAN);
  filter_echo1 = echo1_filt->fslot;

/* allocate standard prescan filters */
@inline Prescan.e PSfilter 




  /**************************************************************/
  /* PREDOWNLOAD:  pulse timing calculations                    */
  /**************************************************************/

  fprintf(stderr,"\tstart_rf1 = %d\n", start_rf1);
  fprintf(stderr,"\trf1_dur = %d\n", rf1_dur);
  fprintf(stderr,"\troinfo.res_kpre = %d\n", roinfo.res_kpre);


  /************************************************************/
  /* PREDOWNLOAD:   set TE and TR                             */
  /************************************************************/

  minte = RUP_GRD( rf1_dur/2 + sspwait_rf1 + MINWAIT + MINWAIT + 4*gyblipinfo.res + 4*gxinfo.res/2);
  temp = _opte.fixedflag; 
  _opte.fixedflag = 0; 
  switch (mode) {
    case 1:                    /* GRE BOLD fMRI */
      opte = 30000;            /* 30 msec TE for BOLD imaging */
      break;
    case 2:                     /* SSFP fMRI */
      opte = minte;
      break;
    case 3:                     /* SSFP fMRI */
      opte = minte;
      break;
    case 4:                     /* RF-SSFP fMRI */
      opte = minte;
      break;
    case 5:                     /* RF-SSFP fMRI */
      opte = minte;
      break;
    case 6:                     /* pdSSFP fMRI */
      opte = minte;
      break;
    default:
      opte = MAX(opte,minte);
  }
  _opte.fixedflag = temp;

  start_ro = RUP_GRD(start_rf1 + psd_rf_wait + MINPW + rf1_dur/2 + opte - 4*gxinfo.res/2 - 4*gyblipinfo.res);

  mintr = RUP_GRD(start_ro + 4*(gyblipinfo.res + gxinfo.res + gyblipinfo.res) + MAXWAIT + 6*MINWAIT + 4*ssi_time + opuser16);
  if(!isSSFP) 
    mintr = RUP_GRD(mintr + 4*crusherinfo.res);

  temp = _optr.fixedflag; 
  _optr.fixedflag = 0; 
  switch (mode) {
    case 1:                     /* GRE BOLD fMRI */
      optr = mintr;
      break;
    case 2:                     /* SSFP fMRI */
      optr = mintr;
      break;
    case 3:                     /* SSFP fMRI */
      optr = mintr;
      break;
    case 4:                     /* RF-SSFP fMRI */
      optr = mintr;
      break;
    case 5:                     /* RF-SSFP fMRI */
      optr = mintr;
      break;
    case 6:                     /* pdSSFP fMRI */
      optr = mintr;
      break;
    default:
      optr = MAX(optr,mintr);
  }
  _optr.fixedflag = temp;

  /* number of phases (scans) */
  temp = _opfphases.fixedflag; 
  _opfphases.fixedflag = 0;
  opfphases = ROUND(scandur*1000000.0 / (1.0*optr*nviews*myopzres));
  _opfphases.fixedflag = temp;

  /* set scan clock */
  pitscan = optr*nviews*myopzres*opfphases + optr*ndisdaq;

  /* print timing parameters */
  fprintf(stderr, "\nTiming values:\n");
  fprintf(stderr, "\tpsd_rf_wait   = %d\n", psd_rf_wait);
  fprintf(stderr, "\tpsd_grd_wait  = %d\n", psd_grd_wait);
  fprintf(stderr, "\tpitscan       = %.1f sec\n", pitscan);
  fprintf(stderr, "\n\n");



  /****************************************************************/
  /* PREDOWNLOAD: entry point table                               */
  /****************************************************************/

  if(entrytabinit(entry_point_table, (int)ENTRY_POINT_MAX)  == FAILURE) 
    { epic_error(0,"failed to initialize entry point table in predownload.",0,0);
      return FAILURE; }
  
  /* set up entry points for scan & prescan */
  entry_point_table[L_APS2] =
    entry_point_table[L_MPS2] = entry_point_table[L_SCAN];
  entry_point_table[L_SCAN].epfilter = (unsigned char)scanslot;
  entry_point_table[L_APS2].epfilter = (unsigned char)prescanslot;
  entry_point_table[L_MPS2].epfilter = (unsigned char)prescanslot;
  entry_point_table[L_SCAN].epprexres = rhfrsize;
  entry_point_table[L_APS2].epprexres = 256;
  entry_point_table[L_MPS2].epprexres = 256;
  strcpy(entry_point_table[L_SCAN].epname,"scan");
  strcpy(entry_point_table[L_APS2].epname,"aps2");
  strcpy(entry_point_table[L_MPS2].epname,"mps2");



	/******************************************************************/
	/* PREDOWNLOAD: RF gain                                           */
	/******************************************************************/

	/* First, find the peak B1 for the whole sequence. */
	maxB1Seq = 0.0;
	for (entry=0; entry < MAX_ENTRY_POINTS; entry++) {
		if (peakB1(&maxB1[entry], entry, RF_FREE, rfpulse) == FAILURE) {
			epic_error(0, "peakB1() failed in predownload.",0,0);
			return FAILURE;
		}
		if (maxB1[entry] > maxB1Seq)
			maxB1Seq = maxB1[entry];
	}

	/* Set xmtadd according to maximum B1 and rescale for powermon,
	*  adding additional (audio) scaling if xmtadd is too big.
	*  Add in cfcoilatten, too. 
	*/
	xmtaddScan = -200*log10(maxB1[L_SCAN]/maxB1Seq) + cfcoilatten;
	if (xmtaddScan > cfdbmax) {
		extraScale = (float) pow(10.0, (cfdbmax - xmtaddScan)/200.0);
		xmtaddScan = cfdbmax; 
	} 
	else extraScale = 1.0;

	if (setScale(L_SCAN, RF_FREE, rfpulse, maxB1[L_SCAN], 
	       extraScale) == FAILURE) {
		epic_error(0, "setScale() failed in predownload.",0,0);
		return FAILURE; 
	}
            
	/* ia_rf1 = max_pg_iamp*(*rfpulse[RF1_SLOT].amp); */
	entry_point_table[L_SCAN].epxmtadd = (short) rint((double)xmtaddScan);

	if(orderslice(TYPNCAT, (int)acqs, (int)1, TRIG_INTERN) == FAILURE)
		epic_error(0,"orderslice call failed in predownload.",0,0);

	my_orderslice();


@inline loadrheader.e rheaderinit   /* Recon variables */



  /********************************************************************/
  /* PREDOWNLOAD: recon header variables                              */
  /********************************************************************/

  /* RH variables related to looping */
  /* rhnecho = opnecho; */
  rhnecho = myopzres;                   /* z partition encodes are stored in loaddab echo slot */
  rhnslices = opfphases;                /* temporal phases are stored in loaddab slice slot */
/*
  rhpasframe = rhnframes * rhnecho * rhnslices;
  rhscnframe = rhpasframe * rhnpasses;
*/
  rhbline = 0;
  /* rhrecon = 617; */    /* use /usr/g/bin/recon617 recon script - save Pfiles to
  my directory, instead of /usr/g/mrraw */

  /* RH variables containing bitflags */
  /* rhrawsize = (int)(2*rhptsize*rhfrsize*(rhnframes+rhbline+1)*rhnecho*rhnslices*opfphases); */
  numrecv = rhdab0e-rhdab0s+1;    /* from splx4rf */
  rhrawsize = (int)(2*rhptsize*rhfrsize*(rhnframes+rhbline+1)*rhnecho*rhnslices*numrecv);
  rhdacqctrl = 0;   /* don't flip echoes */
  rhrcctrl = 1;     /* make an image header */
  rhtype = 1;       /* data is not chopped */
  autolock = 1;

  /* RH user variables - values to be written to Pfile */
  rhuser1 = (float)mode;
  rhuser2 = (float)opflip;
  rhuser3 = (float)myopzres;
  rhuser4 = (float)isSSFP;
  rhuser5 = (float)rf_cycle;
  rhuser6 = (float)rfindex;
  rhuser7 = (float)readoutindex;
  rhuser8 = (float)opte;
  rhuser9 = (float)optr;
  rhuser10 = (float)opfphases;
  rhuser11 = (float)scaleYfov;
  rhuser12 = (float)scaleZfov;
  rhuser13 = (float)scalegssx;
  rhuser14 = (float)scalegssy;
  rhuser15 = (float)multgxpd;
  rhuser16 = (float)opxres;
  rhuser17 = (float)opyres;

  scalerotmats(rsprot, &loggrd, &phygrd, (int)(opslquant), obl_debug);

@inline Prescan.e PSpredownload

  fprintf(stderr, "\n\tother: \n");
  fprintf(stderr, "\t\trhnslices = %d\n", rhnslices);
  fprintf(stderr, "\t\topnecho   = %d\n", opnecho);
  fprintf(stderr, "\t\trhnecho   = %d\n", rhnecho);
  fprintf(stderr, "\t\trhnframes = %d\n", rhnframes);

  return SUCCESS;

} /* End-Of-Predownload */

@inline Prescan.e PShost


/* end of HOST section */






































/************************************************************************
 *                                                                      *
 *                   IPG (PULSEGEN) SECTION 				*
 *                                                                      *
 * Write here the functional code that loads hardware sequencer         *
 * memory with data that will allow it to play out the sequence.        * 
 * These functions call pulse generation macros previously defined      *
 * with @pulsedef, and must return SUCCESS or FAILURE.                  *
 *                                                                      *
 ************************************************************************/
@pg
/* #include "rth_keeuk.h" */

/* acquisition pulse */
WF_PULSE echo1   = INITPULSE; 

/* ssp pulses */
WF_PULSE seqcore;  SEQUENCE_ENTRIES  off_seqcore;
WF_PULSE pass;     SEQUENCE_ENTRIES  off_pass;
WF_PULSE endpass = INITPULSE;
#ifndef IPG
  int idx_seqcore;
  int idx_pass;
#endif /* !IPG */

/* slice-select / crusher pulses */
WF_PULSE rf1      = INITPULSE;
WF_PULSE gss,gssx,gssy ;      
WF_PULSE crusher;   /* at end of TR   -- to minimize TR */

/* readout and partial dephasing pulses */
WF_PULSE gx ;  
WF_PULSE gxpd1 ;  
WF_PULSE gxpd2 ;  
WF_PULSE gy ;  
WF_PULSE gyblip1 ;  
WF_PULSE gyblip2 ;  
WF_PULSE gzblip1 ;  
WF_PULSE gzblip2 ;  

WF_PULSE proto = INITPULSE;


/* see also rthawk/rth_pgen.c */
WF_PULSE jfn_makepulse(WF_PROCESSOR wfp, char *pname, char *fname, gradfilehdr ginfo, int tbeg, float mult)
{
	WF_PULSE *ep;
	short *wavespace;
	float target;

#ifdef SIM
	fprintf(stderr, "::ext_gewav: pname = %s, fname = %s, tbeg = %d, mult = %.2f ... ", pname, fname, tbeg, mult);
	fflush(stderr);
#endif
  
	ep = (WF_PULSE *) AllocNode(sizeof(WF_PULSE));
	memcpy((char*) ep, (char*) &proto, sizeof(WF_PULSE));
	pulsename(ep, pname);
	createreserve(ep, wfp, ginfo.res);
  
	wavespace = (short*)AllocNode(ginfo.res*sizeof(short));
	jfn_readextgrad(fname, wavespace, ginfo, mult);
	movewaveimm(wavespace, ep, 0, ginfo.res, TOHARDWARE);
	createinstr(ep, tbeg, 4*ginfo.res, (short)(mult*MAX_PG_IAMP));

	FreeNode(wavespace);

#ifdef SIM
	fprintf(stderr, "OK\n");
#endif

	return(*ep);
}

WF_PULSE jfn_rf_makepulse(WF_PROCESSOR wfp, char *pname, rfstruct *rfinfo, int tbeg)
{
   WF_PULSE *ep;
   WF_PULSE proto = INITPULSE;
   int res = rfinfo->res;
   short* wave;

   ep = (WF_PULSE *) AllocNode(sizeof(WF_PULSE));
   memcpy((char*) ep, (char*) &proto, sizeof(WF_PULSE));
   pulsename(ep, pname);
   createreserve(ep, wfp, res);

   switch (wfp) {
      case TYPXGRAD:
         wave = rfinfo->gx[0];
         break;
      case TYPYGRAD:
         wave = rfinfo->gy[0];
         break;
      case TYPZGRAD:
         wave = rfinfo->gz[0];
         break;
      case TYPRHO1:
         wave = rfinfo->rho[0][0];
         break;
      case TYPTHETA:
         wave = rfinfo->theta[0][0];
         break;
   }

   movewaveimm(wave, ep, 0, res, TOHARDWARE);
   if (wfp == TYPRHO1)
      addrfbits(ep, 0, tbeg, 4*rfinfo->res);
   createinstr(ep, tbeg, 4*res, MAX_PG_IAMP);

   return(*ep);
}


/************************************************************************/
/*             		    PULSEGEN               		        */
/************************************************************************/
void pulsegen(void)
{
	char fname[80];
	int i, j;
	rfstruct rfinfo;   /* contains rf header info, and the waveforms themselves */

	float dx;  /* in-plane shift */
	float dy;

	fprintf(stderr, "\nwelcome to pulsegen!\n");
	psd_board_type = PSDCERD;
	sspinit(psd_board_type);


	/*******************************************************************/
	/* PULSEGEN: RF pulse and slice-select gradients                   */ 
	/*******************************************************************/

	fprintf(stderr, "\tpulsegen: generating 3d rf waveform (psd_rf_wait = %d)\n", psd_rf_wait);

	/* read rf file header, allocate memory, and load waveforms */
	jfn_rf_readheader(rffilename, &rfinfo);
	jfn_rf_allocatemem(&rfinfo);
	jfn_rf_readwaveforms(&rfinfo, 0);

	/* create excitation pulses */
	rf1  = jfn_rf_makepulse(TYPRHO1,  "rf1",   &rfinfo, RUP_GRD(start_rf1+psd_rf_wait));
	gss  = jfn_rf_makepulse(TYPZGRAD, "gss",   &rfinfo, RUP_GRD(start_rf1));
	gssx = jfn_rf_makepulse(TYPXGRAD, "gssx", &rfinfo, RUP_GRD(start_rf1));
	gssy = jfn_rf_makepulse(TYPYGRAD, "gssy", &rfinfo, RUP_GRD(start_rf1));
	setiamp(RUP_GRD(scalegssx*MAX_PG_IAMP),  &gssx,  0);
	setiamp(RUP_GRD(scalegssy*MAX_PG_IAMP),  &gssy,  0);

	/* OK to free local waveform memory, since we're done creating waveforms in hardware */
	jfn_rf_freemem(&rfinfo);

	/* gz crusher at end of TR (if GRE) */
	if (!isSSFP) {
		sprintf(fname, "%scrusher.jfnwav", directory);
		crusher = jfn_makepulse(TYPZGRAD, "crusher", fname, crusherinfo, RUP_GRD(optr-4*ssi_time-4*crusherinfo.res-MINWAIT), 1.0);
	}


	/********************************************************************/
	/* PULSEGEN: gy readout, phase-encode blips, and wait pulses        */
	/********************************************************************/

	fprintf(stderr, "\tpulsegen: generating y gradients and wait pulses \n");

	WAIT(TYPYGRAD, ywait1, RUP_GRD(pend(&rf1,"rf1",0)+sspwait_rf1), MINPW);

	sprintf(fname, "%s_gyblip.jfnwav", readoutfilename);
	gyblip1 = jfn_makepulse(TYPYGRAD, "gyblip1", fname, gyblipinfo, start_ro, 1.0);

	sprintf(fname, "%s_gy.jfnwav", readoutfilename);
	gy = jfn_makepulse(TYPYGRAD, "gy", fname, gyinfo, RUP_GRD(pend(&gyblip1,"gyblip1",0)+MINWAIT), 1.0);

	sprintf(fname, "%s_gyblip.jfnwav", readoutfilename);
	gyblip2 = jfn_makepulse(TYPYGRAD, "gyblip2", fname, gyblipinfo, RUP_GRD(pend(&gy,"gy",0)+MINWAIT), -1.0);

	WAIT(TYPYGRAD, ywait2, RUP_GRD(pend(&gyblip2,"gyblip2",0)+MINWAIT), MAXWAIT);



	/********************************************************************/
	/* PULSEGEN: gx readout, partial dephasing, and wait pulses         */
	/********************************************************************/

	fprintf(stderr, "\tpulsegen: generating x gradients and wait pulses \n");

	WAIT(TYPXGRAD, xwait1, RUP_GRD(pbeg(&ywait1,"ywait1",0)), MINPW);

	sprintf(fname, "%s_gx.jfnwav", readoutfilename);
	gx = jfn_makepulse(TYPXGRAD, "gx", fname, gxinfo, RUP_GRD(pbeg(&gy,"gy",0)), 1.0);

	sprintf(fname, "%s_gxpd.jfnwav", readoutfilename);
	gxpd1 = jfn_makepulse(TYPXGRAD, "gxpd1", fname, gxpdinfo, RUP_GRD(pbeg(&gx,"gx",0)-4*gxpdinfo.res-MINPW), multgxpd);
	gxpd2 = jfn_makepulse(TYPXGRAD, "gxpd2", fname, gxpdinfo, RUP_GRD(pend(&gx,"gx",0)+MINWAIT), multgxpd);

	WAIT(TYPXGRAD, xwait2, RUP_GRD(pbeg(&ywait2,"ywait2",0)), MAXWAIT);

	fprintf(stderr, "\tpulsegen: done with xwait2 \n");




	/********************************************************************/
	/* PULSEGEN: z partition encoding waveforms and wait pulses         */
	/********************************************************************/

	fprintf(stderr, "\tpulsegen: generating z partition encoding gradients and wait pulses \n");

	WAIT(TYPZGRAD, zwait1, RUP_GRD(pbeg(&ywait1,"ywait1",0)), MINPW);

	sprintf(fname, "%s_gyblip.jfnwav", readoutfilename);
	gzblip1 = jfn_makepulse(TYPZGRAD, "gzblip1", fname, gyblipinfo, RUP_GRD(pbeg(&gyblip1,"gyblip1",0)), 1.0);
	gzblip2 = jfn_makepulse(TYPZGRAD, "gzblip2", fname, gyblipinfo, RUP_GRD(pbeg(&gyblip2,"gyblip2",0)), 1.0);

	WAIT(TYPZGRAD, zwait2, RUP_GRD(pbeg(&ywait2,"ywait2",0)), MAXWAIT);


	/*******************************************************************/
	/* PULSEGEN: data acquisition pulse                                */
	/*******************************************************************/

	/* delay by psd_grad_wait to synch with gradients */
	pulsename(&echo1,"echo1");
	acqq(&echo1, (long) RUP_GRD(pbeg(&gx,"gx",0)+4*gxinfo.res_kpre+psd_grd_wait),(long)(DEFAULTPOS),
			(long)(DEFAULTPOS),(long)filter_echo1,
			(TYPDAB_PACKETS)DABNORM);

	/* ssp wait pulse 1, before daq */
	WAIT(SSP,  sspwait1, RUP_GRD(pbeg(&ywait1,"ywait1",0)), MINPW); 

	/* ssp wait pulse 2, after daq */
	WAIT(SSP, sspwait2, RUP_GRD(pbeg(&ywait2,"ywait2",0)), MAXWAIT);


	fprintf(stderr, "\tpulsegen: done daq pulse \n");

		

	/********************************************************************/
	/* PULSEGEN: SSP Packets                                            */
	/********************************************************************/

	fprintf(stderr, "\tpulsegen: generating ssp packets\n");

	/* createseq says that we're done defining a sequence */
	fprintf(stderr, "\t\tpulsegen: generating ssp packets  -- seqcore\n");
	pulsename(&seqcore,"seqcore");
	createseq(&seqcore,optr-4*ssi_time,off_seqcore);


@inline Prescan.e PSpulsegen


	/* endpass is a SSP pulse that tells the signa system
	that some data is available to write to file and/or
	begin reconstruction. 49ms delay is just to avoid 
	timing problems. */

	fprintf(stderr, "\t\tpulsegen: generating ssp packets  -- endpass\n");
	pulsename(&endpass,"endpass");
	createpass(&endpass,(long)(end_pass_time - 1ms)); 

	/* This is to define the sequence with the SSP pulse above.
	Since the pulse is very short, it will be done less than
	1ms after the 49ms delay above.  */

	fprintf(stderr, "\t\tpulsegen: generating ssp packets  -- pass\n");
	pulsename(&pass,"pass");
	createseq(&pass,end_pass_time,off_pass); 


	/* munches on the sequences to create the instructions the
	hardware actually runs */
	buildinstr();


	fprintf(stderr, "pulsegen: done!\n");

	return;

} /* End of pulsegen */



/* For Prescan: Pulse Generation functions */
@inline Prescan.e PSipg

/* End of @pg */







































/************************************************************************
 *             		 REAL-TIME SEQ. PROCESS SECTION 	        *
 *		(The signals for the actual IPG hardware)		*
 * Write here the functional code for the real time processing (IPG     *
 * side). You may declare standard C variables, but of limited types    *
 * short, int, long, float, double, and 1D arrays of those types.       *
 ************************************************************************/
@rsp

/* For IPG Simulator: will generate the entry point list in the IPG tool */
CHAR *entry_name_list[ENTRY_POINT_MAX] = 
	{"scan", 
	 "aps2",
	 "mps2",
@inline Prescan.e PSeplist

int *rf1_freq;
int *receive_freq1;








@rspvar

short viewtable[513];
short zviewtable[128];

int rcphase;

extern PSD_EXIT_ARG psdexitarg;

int slice, dabop, excitation;

int rspent, rspdda, rspbas, rspvus, rspgy1, rspasl;
int rspesl, rspchp, rspnex, rspslq, rspsct;

/* gating parameters */
int phase; 							/* passed to echo slot in loaddab */
short rsp_hrate;      			/* ECG heart rate in beats per minute */
float RRtime;   					/* RR interval time (ms) */
float RRpct = 0.9;				/* percentage of RR time used for acquisition */
int ncardiacphases; 				/* to be calculated on-the-fly */
float timeperphase;    			/* time to acquire one cardiac phase (ms) */
float rate;                   /* used to calculate ncardiacphases */

@inline Prescan.e PSrspvar 










@rsp
#include "pgen_tmpl.h"
#include "epic_loadcvs.h"


STATUS myphasetable(short *phaseTable, int phaseRes) {
	int amax = 32766;
	int amin = -32766;
	int i, ymin, dy;

	/* calculate dy */
	dy = (int) ROUND((amax-amin)/(phaseRes-1.0));
	if (dy * phaseRes/2 - dy/2 > amax) dy--;
	dy &= 0xfffe;  

	/* fill phaseTable array */
	ymin = -dy * phaseRes/2 + dy/2 ;
	ymin &= 0xfffe;
	for(i = 0; i < phaseRes; i++) 
		phaseTable[i] = (short) (ymin + i*dy);

	/* scale y FOV */
/*
	for(i = 0; i < phaseRes; i++) 
		phaseTable[i] = (short) ((float)phaseTable[i] * scaleYfov);
*/

	return SUCCESS;
}


STATUS psdinit(void) 
{
    /* Initialize everything to a known state */
    view = slice = excitation = 0;
    setrfconfig((short) 5);	/* only activate rho1 */
    setssitime(ssi_time);	/* set ssi counter to cv ssi_time */
    rspqueueinit(200);	        /* initialize to 200 entries */
    scopeon(&seqcore);	        /* Activate scope for core */
    syncon(&seqcore);	        /* Activate sync for core */
    syncoff(&pass);		/* Deactivate sync during pass */

    settriggerarray((short)slquant1, rsptrigger);
    setrotatearray((short)slquant1,rsprot[0]);

    /* set filter based on whether prescan or scan called */
    if ((rspent == L_APS2) || (rspent==L_MPS2)) { /* prescan */
      setrfltrs((int)prescanslot,&echo1);
    } else { /* scan */
      setrfltrs((int)scanslot,&echo1);
    }

    boffset(off_seqcore);	/* start the hardware in the 'core' sequence */

    return SUCCESS;
}

@inline Prescan.e PScore


/* For Prescan: MPS2 Function */
STATUS mps2(void) {
    /* Initialize RSP parameters */
    rspent = L_MPS2;	
    rspdda = 2;
    rspbas = 0;
    rspvus = 30000;
    rspgy1 = 0;
    rspnex = 2;
    rspesl = -1;
    rspasl = 0;
    rspslq = slquant1;
    rspsct = 0;
    
    if (psdinit() == FAILURE) return rspexit();
    
    prescanCore();
    return rspexit();
} /* End of MPS2 */


/* For Prescan: APS2 Function */
STATUS aps2(void) 
{
    /* Initialize RSP parameters */
    rspent = L_APS2;	
    rspdda = 2;
    rspbas = 0;
    rspvus = 1024;
    rspgy1 = 0;
    rspnex = 2;
    rspesl = -1;
    rspasl = -1;
    rspslq = slquant1;
    rspsct = 0;
    
    if (psdinit() == FAILURE) return rspexit();
    
    prescanCore();
    return rspexit();
} /* End of APS2 */



void setfrequency_off_x(float freq, WF_PULSE* prf)
{
	float freq_offset = 1000.0 * JFN_GAMMA * a_gxw * rsp_info[0].rsprloc / (10 * TARDIS_FREQ_RES);
	float freq_final = freq + freq_offset;
  	setfrequency(freq_final, prf, 0);
}




/************************************************************************/
/*                          SCAN FUNCTION                               */
/************************************************************************/
STATUS scan(void)
{
	int nexcnt;
	int disdaqcnt;
	int view, zview, phase, pass;

	fprintf(stderr, "\nwelcome to scan!\n");

	/* RPCstart_server(); */	/* start keeukRPC */
	
	rspent = L_SCAN;
	if (psdinit() == FAILURE) return rspexit();
	setrotatearray(opslquant,rsprot[0]);
	settriggerarray(opslquant,rsptrigger);

	settrigger(TRIG_INTERN,0);

	/* set transmit frequency */
	rf1_freq = (int *)AllocNode(opslquant*sizeof(int));
	setupslices(rf1_freq, rsp_info, opslquant, gssamp, (float)1,
		opfov,TYPTRANSMIT);
	setfrequency(rf1_freq[0], &rf1, 0);

	/* set receive frequency */
	receive_freq1 = (int *)AllocNode(opslquant*sizeof(int));
	(*receive_freq1) = (int)(TARDIS_FREQ_OFFSET/TARDIS_FREQ_RES);
	/* setfrequency_off_x(receive_freq1[0], &echo1); */
	setfrequency(receive_freq1[0], &echo1, 0);

	/* set RF amplitude */
    setiamp(ia_rf1, &rf1, 0); 

	myphasetable(viewtable, gxinfo.npix);
	myphasetable(zviewtable, myopzres);


	/********************************/
	/* SCAN: acquire data           */

	fprintf(stderr,"scan: starting acquisition loop\n");
    
	slice = 0;
	setfrequency(rf1_freq[slice], &rf1, 0);

	rf_phase = 0;    
	rf_spoil_seed_cnt = 0;

	/* establish steady-state */
	for (disdaqcnt = 0; disdaqcnt < ndisdaq; disdaqcnt++) {
		doread_single(0, 0, 0, 0);
	}

	/* Acquire multiple sequential Pfiles */

	for (pass = 0; pass < npasses; pass++) {

		boffset(off_seqcore);

		for (phase = 0; phase < opfphases; phase++) {

			for (zview = 0; zview < myopzres; zview++) {
				if (!isSSFP || roinfo.nechoes==1) {
					/* no need to group leaves for GRE or 2DFT imaging */
					for (view = 0; view < rhnframes; view++) 
						doread_single(view, zview, phase, pass);
				} else {
					/* group views to minimize eddy-current artifacts */
					for (view = 0; view < rhnframes; view+=2) 
						doread_single(view, zview, phase, pass);
					for (view = 1; view < rhnframes; view+=2) 
						doread_single(view, zview, phase, pass);
				}
			}
		}
	}   /* do another pass */

	/* switch to "off_pass" sequence, which writes data out to Pfile */
	boffset(off_pass);
	/*
	if (pass < (rhnpasses-1)) 
		setwamp(SSPD+DABPASS, &endpass, 2);
	else 
	*/
	setwamp(SSPD+DABPASS+DABSCAN, &endpass, 2);

	startseq(0, (short)MAY_PAUSE);

	fprintf(stderr,"scan: done!\n");

	FreeNode(rf1_freq);
	FreeNode(receive_freq1);

	return rspexit();

}  /* end-of-Scan */



/* acquire single frame of data */
void doread_single(int view, int zview, int phase, int pass) {
	int dabecho  = zview;    /* z partition-encode number */
	int dabslice = phase;    /* temporal frame number */
	int dabview  = view+1;   /* y phase-encode number */
	int echoshift = (int) ROUND(4.0*(float)gxinfo.nechosep*(float)view/(float)rhnframes);
	int gxiamp = (int) (pow(-1.0,(float)view)*MAX_PG_IAMP);
	
	setiamp(- (short) ((float)viewtable[view]*scaleYfov),  &gyblip1,  0);
	setiamp(  (short) ((float)viewtable[view]*scaleYfov),  &gyblip2,  0);
	setiamp(- (short) ((float)zviewtable[zview]*scaleZfov),  &gzblip1,  0);
	setiamp(  (short) ((float)zviewtable[zview]*scaleZfov),  &gzblip2,  0);

	/* If doing segmented EPI, flip every other view to place ghosts at FOV/2,
	and do echo time shifting */
	if (roinfo.nechoes!=1) {
		setiamp((short)gxiamp,  &gx,  0);

	  	setperiod(RUP_GRD(MINPW  +echoshift),&xwait1,0);
  		setperiod(RUP_GRD(MAXWAIT-echoshift),&xwait2,0);
  		setperiod(RUP_GRD(MINPW  +echoshift),&ywait1,0);
  		setperiod(RUP_GRD(MAXWAIT-echoshift),&ywait2,0);
  		setperiod(RUP_GRD(MINPW  +echoshift),&zwait1,0);
  		setperiod(RUP_GRD(MAXWAIT-echoshift),&zwait2,0);
  		setperiod(RUP_GRD(MINPW  +echoshift),&sspwait1,0);
  		setperiod(RUP_GRD(MAXWAIT-echoshift),&sspwait2,0);
	}

	/* set transmit and receive phase */
	if (isSSFP) {
		rf_phase += (float)(rf_cycle/180.0*M_PI);
	}
	else {
		rf_phase += (rf_spoil_seed * M_PI / 180.0) * rf_spoil_seed_cnt;
		rf_spoil_seed_cnt++;
	}

	setrfphase(rf_phase);

	loaddab(&echo1,dabslice,dabecho,DABSTORE,dabview,(TYPDAB_PACKETS)DABON,PSD_LOAD_DAB_ALL);

	/* play sequence */
	startseq(0, (short)MAY_PAUSE);
	settrigger(TRIG_INTERN,0);
}


void setrfphase(float rfphase)
{
	rfphase = atan2 (sin(rfphase), cos(rfphase));     /* wrap phase to (-pi,pi) range */
	setphase(rfphase, &rf1, 0);
	setphase(rfphase, &echo1, 0);
}






/************************************************************************/
/* 			PRE-SCAN FUNCTION				*/
/************************************************************************/
STATUS prescanCore(void) 
{
	int k, view;

	if (psdinit() == FAILURE) return rspexit();
    
	setrotatearray(1, rsprot[0]);
	settriggerarray(1, rsptrigger);
    
	/* set the transmit frequency */
	rf1_freq = (int *)AllocNode(sizeof(int));
	(*rf1_freq) = GAM * gssamp * rsp_info[rspesl].rsptloc / (10*TARDIS_FREQ_RES);
	setfrequency(rf1_freq[0], &rf1, 0);

	/* set the receive frequency */
	receive_freq1 = (int *)AllocNode(sizeof(int));
	(*receive_freq1) = (int)(TARDIS_FREQ_OFFSET/TARDIS_FREQ_RES);
	setfrequency(receive_freq1[0], &echo1, 0);

    
    loaddab(&echo1, 0, 0, DABSTORE, (int) 0, DABON,PSD_LOAD_DAB_ALL);
	for (k = 0; k < ps_ndisdaq; k++)
		startseq(0, (short)MAY_PAUSE);
/*
		doread(0, 0, 0);
*/
    
/*
	for (view = 1; view < rspvus; view++) {
*/
	for (view = 1; view < opyres; view++) {
    	loaddab(&echo1, 0, 0, DABSTORE, (int) view, DABON,PSD_LOAD_DAB_ALL);
		startseq(0, (short)MAY_PAUSE);
	}
/*
		doread(view, 0, 0); 
*/
    
	return rspexit();
}


/******************************************** 
 * dummylinks                                
 *                                           
 * This routine just pulls in routines from  
 * the archive files by making a dummy call.  
 ********************************************/
void dummylinks()                             
{                                             
    epic_loadcvs("thefile");           
}                                            
