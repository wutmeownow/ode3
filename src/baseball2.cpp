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
#include "TAxis.h"
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
  f[3] = -F * v * y[3] - p->B * p->omega * y[1] * p->sinPhi - p->g;    // f_vj
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
  pars.delta=5; // m/s
  pars.B=4.1e-4;
  pars.omega=2*M_PI*1800/60.; // rps


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
    title="slider";
    phi=0.;
    v0=85.;
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
    title="curveball";
    phi=45.;
    v0=85.;
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
    title="screwball";
    phi=135.;
    v0=85.;
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
    title="fastball";
    phi=225.;
    v0=95.;
  }

  // finish setting params
  pars.cosPhi = cos(phi*M_PI/180);
  pars.sinPhi = sin(phi*M_PI/180);
  void *p_par = (void*) &pars;

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  double x0=0.; // distance of release
  double z0=0;  // height of release [m]
  double y0=0.; // depth of release
  double theta0=1; // angle of velocity at release (degrees)
  double vx0=v0* cos(theta0*M_PI/180) * MPH_TO_MPS; // initial x vel m/s
  double vz0=v0* sin(theta0*M_PI/180) * MPH_TO_MPS; // initial z vel m/s
  double vy0=0.;  // initial y vel m/s
  double xend=18.44;   // meters to plate
  double yend=0;    // tbd
  double zend=0;    // tbd
  double vxend=0;
  double vyend=0;
  double vzend=0;

  // write code here
  // set initial conditions
  double y[6]={x0, vx0, z0, vz0, y0, vy0};
  double t = 0.0;
  double dt = 0.0001; // ms time step

  // canvas for graph
  TCanvas *c2 = new TCanvas("c2","baseball2 gsl",600,600);
  TGraph *tzvsx = new TGraph();
  tzvsx->SetTitle(title + ";x (ft);z (ft)");
  TGraph *tyvsx = new TGraph();
  tyvsx->SetTitle(title + ";x (ft);y (ft)");

  

  // Init the ODE solver, by defining the system of equations.
  // The second parameter below is the Jacobian (a function if provided to calcualte
  // partial derivatives of the functions).  This is not used for RK algorithms,
  // so we set the function pointer to 0.
  // The params array allows us to set parameters in the functions.
  // {func, jacoboan fcn, dimension, params (can be 0 if not used)}
  gsl_odeiv2_system sys = {func, 0, 6, p_par};

  // set gsl driver
  gsl_odeiv2_driver *drv =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-9, 1e-9, 0.0);


  // set initial points
  tzvsx->SetPoint(tzvsx->GetN(),y[0]* M_TO_FT,y[2]* M_TO_FT);
  tyvsx->SetPoint(tyvsx->GetN(),y[0]* M_TO_FT,y[4]* M_TO_FT);
  // integrate until x reaches xend
  while (y[0] < xend) {
    int status = gsl_odeiv2_driver_apply(drv, &t, t + dt, y);
    if (status != GSL_SUCCESS) { break; }
    tzvsx->SetPoint(tzvsx->GetN(),y[0]* M_TO_FT,y[2]* M_TO_FT);
    tyvsx->SetPoint(tyvsx->GetN(),y[0]* M_TO_FT,y[4]* M_TO_FT);
  }

  // graph
  // graph
  tzvsx->Draw("AC");
  tzvsx->GetXaxis()->SetLimits(0,xend*M_TO_FT);
  tzvsx->GetYaxis()->SetRangeUser(-4,2);
  tzvsx->Draw("AC");
  tyvsx->SetLineStyle(2);
  tyvsx->Draw("L SAME");
  c2->Print(title + "GSL.pdf");

  // set final parameters for comparison
  xend = y[0] * M_TO_FT;
  zend = y[2] * M_TO_FT;
  yend = y[4] * M_TO_FT;
  vxend = y[1] / MPH_TO_MPS;
  vzend = y[3] / MPH_TO_MPS;
  vyend = y[5] / MPH_TO_MPS;
  double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]) / MPH_TO_MPS;


  // to compare to the plots in Fitzpatrick, output your results in **feet**
  // do not change these lines
  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,z) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("(v) = (%lf)\n",v);
  printf("(t) = (%lf)\n",t);
  printf("********************************\n");

  // plot the trajectory.  See Fitzpatrick for plot details
  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

