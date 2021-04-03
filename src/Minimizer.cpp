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
 * \date 13/1/2012  Santa Clara
 */
//==============================================================================
  #include "Minimizer.h"
  #include <cmath>
  #include <iostream>
  #include <iomanip>
  #include <ostream>
  #include <cstdlib>
  #include <cstdio>

//  #define DEBUG
//==============================================================================
// Subroutines
//==============================================================================

//---!!---------------------
  Minimizer::Minimizer(void)
//---!!---------------------
  {
    NoiseModel  *noisemodel = NoiseModel::getInstance();

    using namespace std;
    Nparam = noisemodel->get_Nparam();
    param  = new double[Nparam];
    error  = new double[Nparam];
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
//---------------------------------------------
  double my_f(const gsl_vector *v, void *dummy)
//---------------------------------------------
  {
     Likelihood    *likelihood=Likelihood::getInstance();

     return likelihood->compute(v->data);
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
    double        sigma_eta,lndeterminant,*theta=NULL,FTOL=1.0e-4;
    NoiseModel    *noisemodel  = NoiseModel::getInstance();
    Likelihood    *likelihood  = Likelihood::getInstance();

    //--- GSL minimizer variables
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer            *s = NULL;
    gsl_vector                         *ss, *x;
    gsl_multimin_function              minex_func;
    size_t                             iter = 0;
    int                                status;
    double                             size;

    using namespace std;
    //--- Check if simple ordinary least-squares is required
    if (Nparam==0) {
      //--- no Minimizer required, simple OLS
      cout << "performing Ordinary Least-Squares" << endl;
      likelihood->compute_LeastSquares(param,lndeterminant,sigma_eta);
      show_BIC();
      cout << endl;
    } else {    
      //--- Find minimum using Nelder-Mead Simplex method ---------

      //--- Starting point 
      x = gsl_vector_alloc (Nparam);
      for (i=0;i<Nparam;i++) {
        gsl_vector_set (x, i, 0.1);
      }
     
      //--- Set initial step sizes to 0.3 
      ss = gsl_vector_alloc (Nparam);
      gsl_vector_set_all (ss, 0.3);
     
      //--- Initialize method and iterate 
      minex_func.n = Nparam;
      minex_func.f = &my_f;
      minex_func.params = NULL;
     
      s = gsl_multimin_fminimizer_alloc (T, Nparam);
      gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

      do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
           
        if (status) 
          break;
     
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, FTOL);
     
        if (status == GSL_SUCCESS) {
          printf ("converged to minimum at\n");
        }
     
        printf ("%5zd",iter);
        for (i=0;i<Nparam;i++) {
          printf("%11.5f ",gsl_vector_get(s->x,i));
        }
        printf("  f()= %10.6f  size=%.3f\n",s->fval, size);
      } while (status == GSL_CONTINUE && iter < 1500);

      //--- Sanity check    
      if (iter==500) {
        cerr << "Did not converge!!" << endl;
        exit(EXIT_FAILURE);
      }
 
      //--- save minimum
      for (i=0;i<Nparam;i++) param[i] = gsl_vector_get(s->x,i);
  
      gsl_vector_free(x);
      gsl_vector_free(ss);
      gsl_multimin_fminimizer_free (s);

      //--- Compute confidence intervals for noise parameters
      compute_confidence_intervals();

      //--- Show results in a way a human can understand
      noisemodel->show(param,error);
    }

    //--- Show results of least-squares
    likelihood->show_leastsquares();
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
 */
//------------------------------
  void Minimizer::show_BIC(void)
//------------------------------
  {
    int           m;
    double        lnL_min;
    Likelihood    *likelihood   = Likelihood::getInstance();
    Observations  *observations = Observations::getInstance();
   
    using namespace std; 
    //--- Compute log(likelihood) at minimum
    lnL_min = likelihood->compute(param);
    cout << endl << "Likelihood value" << endl << "--------------------"
         << endl;
    //--- Remember, I use -L to find minimum instead of maximum
    m = observations->number_of_observations()-observations->number_of_gaps();
    cout << fixed << setprecision(3);
    cout << "min log(L)=" << -lnL_min << endl;
    cout << "AIC       =" << 2.0*(Nparam + lnL_min) << endl;
    cout << "BIC       =" << Nparam*log(m) + 2.0*lnL_min << endl;
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
    Likelihood    *likelihood   = Likelihood::getInstance();
    NoiseModel    *noisemodel   = NoiseModel::getInstance();

    using namespace std;
    //--- Create array to hold parameters
    X = new double[Nparam];

    //--- Compute log(likelihood) at minimum
    lnL_min = likelihood->compute(param);
    show_BIC();

    for (i=0;i<Nparam;i++) {
      for (j=i;j<Nparam;j++) {
        //--- diagonal terms
        if (j==i) {
          //--- Investigate both sides of the minimum
          fill_X(i, ds, i, 0.0, X); 
          if (noisemodel->compute_penalty(X)>0.0) throw "out of range";
          lnL[0] = likelihood->compute(X);
          fill_X(i,-ds, i, 0.0, X); 
          if (noisemodel->compute_penalty(X)>0.0) throw "out of range";
          lnL[1] = likelihood->compute(X);
          C[i + Nparam*i] = (lnL[0]+lnL[1] - 2.0*lnL_min)/pow(ds,2.0);
        } else {
          fill_X(i, 0.5*ds, j, 0.5*ds, X); 
          if (noisemodel->compute_penalty(X)>0.0) throw "out of range";
          lnL[0] = likelihood->compute(X);
          fill_X(i,-0.5*ds, j, 0.5*ds, X); 
          if (noisemodel->compute_penalty(X)>0.0) throw "out of range";
          lnL[1] = likelihood->compute(X);
          fill_X(i, 0.5*ds, j,-0.5*ds, X); 
          if (noisemodel->compute_penalty(X)>0.0) throw "out of range";
          lnL[2] = likelihood->compute(X);
          fill_X(i,-0.5*ds, j,-0.5*ds, X); 
          if (noisemodel->compute_penalty(X)>0.0) throw "out of range";
          lnL[3] = likelihood->compute(X);
          C[Nparam*j + i] = (lnL[0]-lnL[1]-lnL[2]+lnL[3])/pow(ds,2.0);
        }
      }
    }

    //--- Invert matrix
#ifdef DEBUG
    likelihood->show_matrix("C",C,Nparam,Nparam);
#endif
    clapack_dpotrf(CblasColMajor,CblasUpper,Nparam,C,Nparam);
    clapack_dpotri(CblasColMajor,CblasUpper,Nparam,C,Nparam);
#ifdef DEBUG
    likelihood->show_matrix("C",C,Nparam,Nparam);
#endif

    //--- free up memory
    if (X!=NULL) delete[] X;
  }



/*! Compute some confidence limits for the estimated noise parameters. I
 *  assume that the computed noise parameters correspond to the maximum
 *  likelihood value. 
 */
//--------------------------------------------------
  void Minimizer::compute_confidence_intervals(void)
//--------------------------------------------------
  {
    int      i;
    double   *C=NULL;

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

    //--- free memory
    if (C!=NULL) delete[] C;
  }
