///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
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

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
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
  f[0] = y[1];                         // f_ri
  f[1] = -(p->b * p->d + p->c * p->d * p->d * v) * y[1] / p->m;        // f_vi
  f[2] = y[3];                         // f_rj
  f[3] = -(p->b * p->d + p->c * p->d * p->d * v) * y[3] / p->m - p->g;    // f_vj
  f[4] = y[5];
  f[5] = -(p->b * p->d + p->c * p->d * p->d * v) * y[5] / p->m; // f_vk
  return GSL_SUCCESS;
}


int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.075;   
  pars.b=1.6e-4;  
  pars.c=0.25;
  void *p_par = (void*) &pars;

  double x0=0.; // distance of release
  double z0=1.4;             // height of release [m]
  double y0=0.; // depth of release
  double theta0=1;         // angle of velocity at release (degrees)
  double xend=18.5;       // meters to plate
  double yend = 0.;
  double zend = 0.9;
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays
  

  double vPitch = 0;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here

  // canvas for graph
  TCanvas *c2 = new TCanvas("c2","baseball gsl",600,600);
  TGraph *tvpitch = new TGraph();
  tvpitch->SetTitle(";iteration;V pitch (m/s)");

  // Init the ODE solver, by defining the system of equations.
  // The second parameter below is the Jacobian (a function if provided to calcualte
  // partial derivatives of the functions).  This is not used for RK algorithms,
  // so we set the function pointer to 0.
  // The params array allows us to set parameters in the functions.
  // {func, jacoboan fcn, dimension, params (can be 0 if not used)}
  gsl_odeiv2_system sys = {func, 0, 6, p_par};

  // Search range for pitch speed (reasonable MLB range)
  double v_low = 5.0;
  double v_high = 60.0;

  // adjustable search loop from chatgpt
  for (int iter = 0; iter < 40; iter++) {
    double v_try = 0.5*(v_low + v_high);

    // compute vx0, vz0 for this trial
    double vx0_try = v_try * cos(theta0*M_PI/180);
    double vz0_try = v_try * sin(theta0*M_PI/180);
    double vy0_try = 0.0;

    // reset GSL driver each iteration (required!)
    gsl_odeiv2_driver *drv =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,1e-9, 1e-9, 0.0);
    
    // reset initial conditions
    double t_local = 0.0;
    double y_local[6] = {x0, vx0_try, z0, vz0_try, y0, vy0_try};

    // integrate until x reaches xend or ball hits ground
    while (y_local[0] < xend && y_local[2] > 0.0) {
        int status = gsl_odeiv2_driver_apply(drv, &t_local, t_local + 0.01, y_local);
        if (status != GSL_SUCCESS) { break; }
    }

    double z_at_plate = y_local[2];

    gsl_odeiv2_driver_free(drv);

    // If we hit ground before reaching the plate, treat as too slow
    if (y_local[0] < xend) {
        v_low = v_try;
        tvpitch->SetPoint(iter,iter,0.5*(v_low + v_high));
        continue;
    }

    // Use whether z is above or below the target to adjust search
    if (z_at_plate > zend) {
        // too high -> throw slower
        v_high = v_try;
    } else {
        // too low -> throw harder
        v_low = v_try;
    }
    tvpitch->SetPoint(iter,iter,0.5*(v_low + v_high));
  }

  // vpitch is average of lower and higher velocity
  vPitch = 0.5*(v_low + v_high);

  // graph
  tvpitch->Draw("AC");
  c2->Print("baseballGSL.pdf");

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

