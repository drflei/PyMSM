// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//     MODULE:           PyUNILIB_Interface.h
//
//     Version:          0.B
//     Date:             27/08/2019
//     Author:           Pete Truscott
//     Organisation:     Kallisto Consultancy Ltd, UK
//     Customer:         ESA/ESTEC, NOORDWIJK
//     Contract:         4000117974/16/NL/LF (VALIRENE Project)
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//     LICENCE AND COPYRIGHT/DISTRIBUTION CONDITIONS
//
//     (To be confirmed with Hugh Evans (ESA) and Daniel Heynderickx (DHC))
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//     DESCRIPTION
//
//     Part of an example/demonstration Python wrapper for some of the UNILIB
//     subroutines, specifically those needed for UNILIB example1.  It defines
//     the interfaces to the subroutines UT990, UM510, UT540, UM520, UM536, and
//     UM530, as well as the structs ZGEO, ZVEC and ZDAT.  The interfaces are
//     described as though called from C, which is what is required for
//     interfacing to Cython.  It is assumed that the gfortran compilation adds
//     only one underscore to the Fortran subroutine/function/struct name.
//     If interfaces are required to the UNILIB Fortran common blocks, similar
//     C-type declarations can be included.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ==============================================================================
#include "stdbool.h"

//#define NX140_CONST 16
//
//
// The following are required for access to the structs used by UNILIB:
//
// struct zxyz_
// {
//   double x, y, z;
// };

struct zgeo_
{
  double radius;
  double colat;
  double elong;
};

struct zvec_
{
  double dnrm;
  double rho;
  double theta;
  double phi;
};

// struct zlbl_
// {
//   double finv;
//   double fbmp;
//   double fkauf;
//   double flmi;
//   double falp0;
//   double fphi;
//   double ftim;
//   char   label[1];
//   bool   linv, lbmp, lkauf, llmi, lalp0, lphi, ltim;
// };

// struct zpnt_
// {
//   struct zgeo_ coord;
//   struct zvec_ b;
//   double rcurv;
// };

// struct zseg_
// {
//   struct zpnt_ beg;
//   double arcl;
//   double csalp;
//   double dtbnd;
//   double rkstp[3];
// };

// struct zind_
// {
//   int jbeg;
//   int jend;
//   int jmirpn;
//   int jmirps;
// };

// struct zfln_
// {
//   struct zpnt_ equat;
//   struct zgeo_ footpn;
//   struct zgeo_ footps;
//   struct zvec_ drift;
//   double dtdft;
//   int    kwest, keast;
//   struct zind_ ind;
// };

// struct zimf_
// {
//   char   label[32];
//   int    kinner;
//   int    norder;
//  double ccoeff[16][16];   // Redefined range from nx140 == 16
//   double ccoeff[NX140_CONST][NX140_CONST];
//   double colat, elong;
//   double gmmo;
//   double tzero, epoch, saarot;
// };

// struct zsun_
// {
//   double utdeg;
//   double gha;
//   struct zxyz_ dir;
// };

// struct zemf_
// {
//   double vdst;
//   double wdens, wvel;
//   double vkp;
//   double val;
//   double pdyn, bximf, byimf, bzimf, stdoff;
//   double g1, g2, w1, w2, w3, w4, w5, w6;
//   double trans[3];
//   double tilt;
//   char   label[20];
//   int    kouter;
//   int    ikp;
//   char   lbltns[4];
// };

struct zdat_
{
  double secs;
  double amjd;
  int    iyear, imonth, iday;
  int    ihour, imin, idummy;
};

// struct zatm_
// {
//   double ut;
//   double rzss;
//   double f107a, f107;
//   double apind[7];
//   double fkpx;
//   int    katm;
//   int    kion;
//   int    kyear, kday;
// };

//
// The following define the interfaces to some of the Fortran subroutines.
// (UT990, UM510, UT540, UM520, UM536, and UM530 - used for the UNILIB Ex.1)
//
extern void ut990_ (int*, int*, int*);
extern void um510_ (int*, double*, char*, int*, int*, int);
extern void ut540_ (struct zdat_* );
extern void ut545_ (struct zdat_* );
extern void um520_ (int*, double*, double*, char*, int*, int* , int);
extern void um522_ (double*, int*);
extern void um524_ (int*);
extern void um536_ (struct zgeo_*, struct zgeo_*);
extern void um530_ (struct zgeo_*, struct zvec_*, int*);
extern void ul220_ (struct zgeo_*, double*, int*, double*, double*, double*, double*, double*, double*, int*);
extern void ul225_ (double*, double*, int*, double*, double*, int*);

extern int  get_doy_ (int *, int *, int *);
extern void coord_trans1_ (int *, int *, int *, int *,
                           double *, double *, double *);

