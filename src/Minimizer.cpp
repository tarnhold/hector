/*! \file   Minimizer.cpp
 *  \author Machiel Bos
 *
 * Main class which selects the minization process to be used and provides
 * the subroutine to compute the confidence interval of the estimated
 * parameters by evaluating the Fisher information matrix numerically.
 *
 * This class is only called by the main program and I guess it is only
 * called once because the program is only meant to analyse one station.
 * Therefore, there is no need to make this class a Singleton.
 *
 *  This script is part of Hector 1.9
 *
 *  Hector is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  Hector is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Hector. If not, see <http://www.gnu.org/licenses/>
 *
 */
//==============================================================================
  #include "Minimizer.h"
  #include <cmath>
  #include <iostream>
  #include <iomanip>
  #include <ostream>
  #include <cstdlib>
  #include <cstdio>
  #include <sys/types.h>
  #include <unistd.h>

//  #define DEBUG
//==============================================================================
// Subroutines
//==============================================================================

//---!!---------------------
  Minimizer::Minimizer(void)
//---!!---------------------
  {
    NoiseModel  &noisemodel = NoiseModel::getInstance();
    Control     &control=Control::getInstance();

    Nparam = noisemodel.get_Nparam();
    param  = new double[Nparam];
    error  = new double[Nparam];

    using namespace std;
    //--- Do we randomise first gues of noise parameter values?
    try {
      randomise_first_guess = control.get_bool("RandomiseFirstGuess");
    }
    catch (exception &e) {
      randomise_first_guess = false;
    }
  }



//---!!----------------------
  Minimizer::~Minimizer(void)
//---!!----------------------
  {
    if (param!=NULL)  delete[] param;
    if (error!=NULL)  delete[] error;
  }



/*! Wrapper to likelihood function
 *
 * The function does not require additional parameters because I work
 * with singletons which get their information independently when they are
 * created.
 */
//---------------------------
  double my_f(double *param_)
//---------------------------
  {
     Likelihood    &likelihood=Likelihood::getInstance();

     return likelihood.compute(param_);
  }
    

 
/*! For given values of p, q, x and H, compute xhat and AR(p),d and MA(q)
 *  parameters.
 * 
 *  Sometimes I don't want to include fractional noise in the computations.
 *  To allow for this, FI is only included if the value of 'd' is not zero
 *  at entry.
 */
//---------------------------
  void Minimizer::solve(void)
//---------------------------
  {
    int           i;
    double        sigma_eta,ln_det_C,*theta=NULL;
    NoiseModel    &noisemodel  = NoiseModel::getInstance();
    Likelihood    &likelihood  = Likelihood::getInstance();

    //--- Nelder & Mead variables
    double        *step,*start,*xmin,ynewlo,reqmin=1.0e-15;
    int           kconvge=10,kcount=1000,ifault,numres,icount;
    long          seed;
    double        size,w;

    //--- Allocate memory for arrays
    step  = new double[Nparam];
    start = new double[Nparam];
    xmin  = new double[Nparam];

    using namespace std;
    //--- Check if simple ordinary least-squares is required
    if (Nparam==0) {
      //--- no Minimizer required, simple OLS
      cout << "performing Ordinary Least-Squares" << endl;
      likelihood.compute_LeastSquares(param);
      //--- Even if no noise parameters are estimated, it is good to show
      //    which noise model was chosen and what its amplitude is. The
      //    error array is thus not filled but this is okay because nothing
      //    is shown.
      sigma_eta = likelihood.get_sigma_eta();
      noisemodel.show(param,error,sigma_eta);
      likelihood.compute_L_and_ICs(param);
      show_L_and_ICs();
      cout << endl;
    } else {    
      //--- Find minimum using Nelder-Mead Simplex method ---------

      //--- Initiate random generator
      boost::random::mt19937 rng;
      boost::random::uniform_real_distribution<double> gen(0.0, 0.25);

      //--- Starting point 
      for (i=0;i<Nparam;i++) {
        if (randomise_first_guess==true) {
          w = gen(rng);
          start[i] = 0.04 + w;
        } else {
          start[i] = 0.10;
        }
        step[i] = 0.3;
      }

      //--- Call Nelder & Mead
      nelmin(my_f,Nparam,start,xmin,&ynewlo,reqmin,step,kconvge,kcount,
						   &icount,&numres,&ifault);
      printf ("%5d",icount);
      for (i=0;i<Nparam;i++) {
        printf("%11.5f ",xmin[i]);
      }
      cout << "\n";
      cout << "  Return code IFAULT = " << ifault << "\n";
      cout << "\n";
      cout << "  Estimate of minimizing value X*:\n";
      cout << "  F(X*) = " << ynewlo << "\n";
      cout << "\n";
      cout << "  Number of iterations = " << icount << "\n";
      cout << "  Number of restarts =   " << numres << "\n";

      //--- save minimum
      for (i=0;i<Nparam;i++)   param[i] = xmin[i];
  
      //--- Compute confidence intervals for noise parameters. This sets
      //    the values for the error-array used later in noisemode->show.
      compute_confidence_intervals(param);

      //--- Show results in a way a human can understand
      sigma_eta = likelihood.get_sigma_eta();
      noisemodel.show(param,error,sigma_eta);
    }

    //--- Show results of least-squares
    likelihood.show_leastsquares();

    delete[] start;
    delete[] step;
    delete[] xmin;
  }



/*! To compute the second derivative of the likelihood near the minimum,
 *  I need to create vectors that are slightly displaced from the
 *  vector, containing noise parameters, at the minimum. I need to 
 *  do that sometimes in two directions. Therefore, in this subroutine
 *  columns k0 and k1 are displaced  by s0 and s1 respectively.
 *
 * \param[in]  k0   : index of vector entry that needs changing
 * \param[in]  k1   : index of vector entry that needs changing
 * \param[in]  s0   : amount of displacement in entry k0
 * \param[in]  s1   : amount of displacement in entry k1
 * \param[in]  param: array of optimal noise parameters
 * \param[out] X    : array of nearly optimal noise parameters
 */
//-----------------------------------------------------------------------
  void Minimizer::fill_X(int k0, double s0, int k1, double s1, double *X)
//-----------------------------------------------------------------------
  {
    int  i;

    using namespace std;
    //--- sanity check
    if (k0<0 || k0>=Nparam || k1<0 || k1>=Nparam) {
      cerr << "k0=" << k0 << ", k1=" << k1 << " : not right!" << endl;
      exit(EXIT_FAILURE);
    }
    
    //--- Start filling array X
    for (i=0;i<Nparam;i++) {
      X[i]=param[i];
      if (i==k0) X[i] += s0;
      if (i==k1) X[i] += s1;
    }
  }



/* Show L_min, BIC and AIC
 * Nparam is a class parameter
 */
//------------------------------------
  void Minimizer::show_L_and_ICs(void)
//------------------------------------
  {
    using namespace std;
    int             k,n;
    double          *H,ln_L,ln_det_I,AIC,BIC,BIC_tp,BIC_c;
    DesignMatrix    *designmatrix = DesignMatrix::getInstance();
    Likelihood      &likelihood = Likelihood::getInstance();
    JSON            &json = JSON::getInstance();
 
    //--- How many parameters are estimated? 
    designmatrix->get_H(n,&H);
    k = Nparam + n + 1;

    ln_L     = likelihood.get_ln_L();
    ln_det_I = likelihood.get_ln_det_I();
    AIC      = likelihood.get_AIC();
    BIC      = likelihood.get_BIC();
    BIC_tp   = likelihood.get_BIC_tp();
    BIC_c    = likelihood.get_BIC_c();

    cout << endl << "Likelihood value" << endl << "--------------------"
         << endl;
    cout << fixed << setprecision(3);
    cout << "min log(L)=" << ln_L << endl;
    cout << "k         =" << n << " + " << Nparam << " + 1 = " << k << endl;
    cout << "AIC       =" << AIC << endl;
    cout << "BIC       =" << BIC << endl;
    cout << "BIC_tp    =" << BIC_tp << endl;
    cout << "BIC_c     =" << BIC_c << endl;
    cout << "ln_det_I  =" << ln_det_I << endl;

    //--- Add to JSON file
    json.write_double("ln_L",ln_L);
    json.write_double("AIC",AIC);
    json.write_double("BIC",BIC);
    json.write_double("BIC_tp",BIC_tp);
    json.write_double("BIC_c",BIC_c);
    json.write_double("ln_det_I",ln_det_I);

    //--- free memory
    delete[] H;
  }



/*! Simply compute curvature
 */
//---------------------------------------------
  void Minimizer::compute_inv_Fisher(double *C)
//---------------------------------------------
  {
    int           i,j,info,k;
    char          Up='L';
    const double  ds=1.0e-3,TINY=1.0e-6;
    double        lnL_min,lnL[4],*X=NULL;
    Likelihood    &likelihood   = Likelihood::getInstance();
    NoiseModel    &noisemodel   = NoiseModel::getInstance();

    using namespace std;
    //--- Create array to hold parameters
    X = new double[Nparam];

    lnL_min = likelihood.compute(param);
    for (i=0;i<Nparam;i++) {
      for (j=i;j<Nparam;j++) {
        //--- diagonal terms
        if (j==i) {
          //--- Investigate both sides of the minimum
          fill_X(i, ds, i, 0.0, X); 
          if (noisemodel.compute_penalty(X)>0.0) throw "out of range";
          lnL[0] = likelihood.compute(X);
          fill_X(i,-ds, i, 0.0, X); 
          if (noisemodel.compute_penalty(X)>0.0) throw "out of range";
          lnL[1] = likelihood.compute(X);
          C[i + Nparam*i] = (lnL[0]+lnL[1] - 2.0*lnL_min)/pow(ds,2.0);
        } else {
          fill_X(i, 0.5*ds, j, 0.5*ds, X); 
          if (noisemodel.compute_penalty(X)>0.0) throw "out of range";
          lnL[0] = likelihood.compute(X);
          fill_X(i,-0.5*ds, j, 0.5*ds, X); 
          if (noisemodel.compute_penalty(X)>0.0) throw "out of range";
          lnL[1] = likelihood.compute(X);
          fill_X(i, 0.5*ds, j,-0.5*ds, X); 
          if (noisemodel.compute_penalty(X)>0.0) throw "out of range";
          lnL[2] = likelihood.compute(X);
          fill_X(i,-0.5*ds, j,-0.5*ds, X); 
          if (noisemodel.compute_penalty(X)>0.0) throw "out of range";
          lnL[3] = likelihood.compute(X);
          C[Nparam*j + i] = (lnL[0]-lnL[1]-lnL[2]+lnL[3])/pow(ds,2.0);
        }
      }
    }

    //--- Invert matrix
    dpotrf_(&Up,&Nparam,C,&Nparam,&info);
    dpotri_(&Up,&Nparam,C,&Nparam,&info);

    //--- Put optimum back in memory in Likelihood class
    lnL_min = likelihood.compute(param);

    //--- free up memory
    if (X!=NULL) delete[] X;
  }



/*! Compute some confidence limits for the estimated noise parameters. I
 *  assume that the computed noise parameters correspond to the maximum
 *  likelihood value. 
 * 
 * \param[i] param : vector with estimated parameters (passed to show_L_and_ICs)
 */
//-----------------------------------------------------------
  void Minimizer::compute_confidence_intervals(double *param)
//-----------------------------------------------------------
  {
    int         i;
    double      *C=NULL;
    Likelihood  &likelihood=Likelihood::getInstance();

    using namespace std;
    //--- Compute covariance matrix
    C = new double[Nparam*Nparam];
    memset(C,0,(Nparam*Nparam)*sizeof(double));

    try {
      compute_inv_Fisher(C);
    }
    catch(const char *str) {
      cout << "Cannot compute Fisher matrix: " << str << '\n';
      for (i=0;i<Nparam;i++) C[i*Nparam + i]=0.0;
    }
       
    //--- compute error array
    for (i=0;i<Nparam;i++) {
      error[i] = sqrt(C[i*Nparam + i]);
    }

    //--- The way I modify the 'param' variable to compute the Fisher 
    //    information matrix in compute_inv_Fisher is dangerous because I rely
    //    on the fact that the last call to likelihood->compute was done
    //    with the optimal values for 'param' but anyway. Now I call 
    //    likelihood->compute(param) last so that the other classes copy 
    //    the optimal value of the estimated parameters into their instances.
    likelihood.compute_L_and_ICs(param);
    show_L_and_ICs();

    //--- free memory
    if (C!=NULL) delete[] C;
  }



/*! Copy values of parameter param to param_
 */
//-----------------------------------------
  void Minimizer::get_param(double *param_)
//-----------------------------------------
  {
    int    i;

    for (i=0;i<Nparam;i++) {
      param_[i] = param[i];
    }
  }
