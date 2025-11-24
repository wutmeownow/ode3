///
/// Starter template for second baseball problem
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

double const M_TO_FT = 3.28084; // conversion for meters to feet
double const MPH_TO_MPS = 0.44704; // conversion for mph to meters per second

using namespace std;

// model drag force from fitzpatrick with fD(v) = -F(v)|v|v * m
// with F(v) = a + b/(1+exp[(v-vd)/delta]) a,b,vd,delta constants
// magnus force is fM = S(v) * omega × v
// a good approximation of S(v) is B = S/m = 4.1x10^-4

struct Params {
  double g;   // acceleration [m/s^2]
  double a; // a,b,vd,delta constants for drag force from fitzpatrick
  double b;
  double vd;
  double delta;
  double B; // coefficient for magnus force
  double omega; // rotation of ball
  double cosPhi; // cosine of angle of rotation vector
  double sinPhi; // sine of angle of rotation vector
};

// first we define our vaiables
// x    = time
// y[0] = position along i (x) axis   ; f_ri = dri/dt => velocity along i axis  
// y[1] = velocity along i (x) axis   ; f_vi = dvi/dt => acceleration along i axis
// y[2] = position along j (z) axis   ; f_rj = drj/dt => velocity along j axis  
// y[3] = velocity along j (z) axis   ; f_vj = dvj/dt => acceleration along j axis
// y[4] = position along k (y) axis   ; f_rj = drj/dt => velocity along k axis  
// y[5] = velocity along k (y) axis   ; f_vj = dvj/dt => acceleration along k axis

// the function for the gsl solver holds the equations of motion of each coordinate
int func (double t, const double y[], double f[], void *params) {
  Params *p = (Params*)params;
  double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
  double F = p->a + p->b / (1 + exp((v - p->vd)/p->delta));
  f[0] = y[1];                         // f_ri
  f[1] = -F * v * y[1] + p->B * p->omega * (y[3]*p->sinPhi - y[5]*p->cosPhi);        // f_vi
  f[2] = y[3];                         // f_rj
  f[3] = -F * v * y[3] + p->B * p->omega * y[1] * p->sinPhi - p->g;    // f_vj
  f[4] = y[5];
  f[5] = -F * v * y[5] + p->B * p->omega * y[1] * p->cosPhi; // f_vk
  return GSL_SUCCESS;
}


int main(int argc, char **argv){

  // we have 6 initial conditions for this problem
  // y[0] = y[2] = y[4] = 0;  // init x,y,z
  // y[1] = v0*cos(theta0);   // vx  "x is line towards the plate
  // y[3] = 0;                // vy  "y" is measured as left/right divergence from line to plate
  // y[5] = v0*sin(theta0);   // vz  "z" is vertival measure
  // vector<double> y0(6);
  double y[6];

  bool showPlot=false;
  // pitches
  // slider ip=0
  // curve ip=1
  // screwball ip=2
  // fast ip=3
  int ip=1;    // default pitch
  int c;

  // set default parameters from fitzpatrick
  Params pars;
  pars.g=9.81;
  pars.a=0.0039;
  pars.b=0.0058;
  pars.vd=35; // m/s
  pars.delta=5 // m/s
  pars.B=4.1e-4;
  pars.omega=1800/60.; // rps


  while ((c = getopt (argc, argv, "p:n")) != -1) {
    switch (c) {
    case 'p':
      ip = atoi(optarg);
      break;
    case 'n':
      showPlot=false;
      break;
    }
  }

  double v0; // launch velocity in mph
  double phi; // rotation vector angle in deg
  TString title;
  if (ip==0){
    cout << "Setting up initial conditions for slider" << endl;
    title="Slider";
    phi=0.;
    v0=85.;
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
    title="Curveball";
    phi=45.;
    v0=85.;
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
    title="Screwball";
    phi=135.;
    v0=85.;
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
    title="Fastball";
    phi=225.;
    v0=95.;
  }

  // finish setting params
  pars.cosPhi = cos(phi*M_PI/180);
  pars.sinPhi = sin(phi*M_PI/180);

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  double x0=0.; // distance of release
  double z0=0;  // height of release [m]
  double y0=0.; // depth of release
  double theta0=1; // angle of velocity at release (degrees)
  double vx0=v0* cos(theta0*M_PI/180) * MPH_TO_MPS; // initial x vel
  double vz0=v0* sin(theta0*M_PI/180) * MPH_TO_MPS; // initial z vel
  double vy0=0.;  // initial y vel
  double xend=18.44;   // meters to plate
  double yend=0;    // tbd
  double zend=0;    // tbd
  double vxend=0;
  double vyend=0;
  double vzend=0;

  // write code here


  // to compare to the plots in Fitzpatrick, output your results in **feet**
  // do not change these lines
  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,x) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("********************************\n");

  // plot the trajectory.  See Fitzpatrick for plot details
  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

