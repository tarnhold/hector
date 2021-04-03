/*! \file    NoiseModel.cpp
 *  \author  Machiel Bos
 *
 * Main interface to various Noise model implementations. It allows the
 * combination of various noise models.
 *
 * \date 13/ 1/2012  Santa Clara
 * \date  7/10/2012  Santa Clara
 * \date 16/ 1/2013  Santa Clara
 */
//=============================================================================
  #include "NoiseModel.h"
  #include "White.h"
  #include "Powerlaw.h"
  #include "PowerlawApprox.h"
  #include "ARFIMA.h"
  #include "GenGaussMarkov.h"
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstring>

//  #define DEBUG

//=============================================================================
// Subroutines
//=============================================================================
//

//++++++ Singleton stuff +++++++++++++++++
bool NoiseModel::instanceFlag = false;
NoiseModel* NoiseModel::singleton = NULL;
//++++++++++++++++++++++++++++++++++++++++


//--!!------------------------
  NoiseModel::NoiseModel(void)
//--!!------------------------
  {
    int      i;
    Control  *control = Control::getInstance();

    using namespace std;
    noisemodel = new char* [5];
    try {
      control->get_name_list("NoiseModels",noisemodel,Nmodels);
    } 
    catch (const char* str) {
      cerr << str << endl;
      exit(EXIT_FAILURE);
    }
    if (Nmodels<1) {
      cerr << "No noise model was specified!" << endl;
      exit(EXIT_FAILURE);
    } else if (Nmodels>5) {
      cerr << "Only up to 5 noise models can be specified!" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Allocate memory space to store each individual noise model
    Nparam         = Nmodels-1;
    modelIndv      = new NoiseModelBaseClass* [Nmodels];
    NparamIndv     = new int[Nmodels];
    phi            = new double[Nparam];
    fraction_fixed = new double[Nparam];
   
    //--- Avoid problems, set h to NULL 
    h = NULL;
 
    //--- create a generator chosen by the environment variable GSL_RNG_TYPE
    gsl_rng_env_setup();
     
    //--- Initiate random generator
    T_random = gsl_rng_default;
    r_random = gsl_rng_alloc (T_random);

    printf ("generator type: %s\n", gsl_rng_name (r_random));
    printf ("seed = %lu\n", gsl_rng_default_seed);
    printf ("first value = %lu\n", gsl_rng_get (r_random));

    //--- Create all necessary classes
    for (i=0;i<Nmodels;i++) {
      if (strcmp(noisemodel[i],"White")==0) {
        modelIndv[i] = new White;
      } else if (strcmp(noisemodel[i],"Powerlaw")==0) {
        modelIndv[i] = new Powerlaw();
      } else if (strcmp(noisemodel[i],"Flicker")==0) {
        modelIndv[i] = new Powerlaw(0.5);
      } else if (strcmp(noisemodel[i],"RandomWalk")==0) {
        modelIndv[i] = new Powerlaw(1.0);
      } else if (strcmp(noisemodel[i],"PowerlawApprox")==0) {
        modelIndv[i] = new PowerlawApprox;
      } else if (strcmp(noisemodel[i],"FlickerApprox")==0) {
        modelIndv[i] = new PowerlawApprox(0.5);
      } else if (strcmp(noisemodel[i],"ARFIMA")==0) {
        modelIndv[i] = new ARFIMA();
      } else if (strcmp(noisemodel[i],"ARMA")==0) {
        modelIndv[i] = new ARFIMA(0.0);
      } else if (strcmp(noisemodel[i],"GGM")==0) {
        modelIndv[i] = new GenGaussMarkov;
      } else {
        cerr << "Unknown noise model: " << noisemodel[i] << endl;
        exit(EXIT_FAILURE);
      }
      NparamIndv[i] = modelIndv[i]->get_Nparam();
      Nparam       += NparamIndv[i];
    }
#ifdef DEBUG
     for (i=0;i<Nmodels;i++) {
       cout << noisemodel[i] << ", Nparam=" << NparamIndv[i] << endl;
     }
     cout << "Nparam total = " << Nparam << endl;
#endif
  }   
  


//--!!-------------------------
  NoiseModel::~NoiseModel(void)
//--!!-------------------------
  {
    int    i;

    if (h!=NULL) { //-- This implies also F_h != NULL
      delete[] h;
      delete[] w;
      fftw_free (F_h);
      fftw_free (F_w);
      fftw_destroy_plan ( plan_forward );
      fftw_destroy_plan ( plan_backward );
    }
    if (Nmodels>0) {
      for (i=0;i<Nmodels;i++) {
        delete[] noisemodel;
      }
      delete[] modelIndv;
      delete[] NparamIndv;
      delete[] phi;
      delete[] fraction_fixed;
    }
    gsl_rng_free (r_random);
  }



/*! For each noise model, get the first column of the covariance matrix,
 *  scale each with phi_i and sum the total. 
 *  The 'param' array is managed by the Minimizer class, 'gamma_x' by the 
 *  Likelihood class. I don't know if this is good C++ style but anyway.
 */
//----------------------------------------------------------------------
  void NoiseModel::get_covariance(double *param, int m, double *gamma_x)
//----------------------------------------------------------------------
  {
    int      i,j,k;
    double   fraction,*gamma_xIndv;

    using namespace std;
#ifdef DEBUG
    cout << "m=" << m << endl;
    for (i=0;i<Nparam;i++) cout << "param[" << i << "]= " << param[i] << endl;
#endif
    gamma_xIndv = new double[m];
    memset(gamma_x,0.0,m*sizeof(double));
    k = Nmodels-1;
    for (i=0;i<Nmodels;i++) {
      modelIndv[i]->get_covariance(&param[k],m,gamma_xIndv);
      
      //--- Compute fraction
      fraction = 1.0;
      if (i<Nmodels-1) {
        for (j=0;j<Nmodels-1;j++) { 
          if (j==i) fraction *= (1.0-param[j]);
          else      fraction *= param[j];
        } 
      } else {
        for (j=0;j<Nmodels-1;j++) fraction *= param[j]; 
      }
#ifdef DEBUG
      cout << "i=" << i << ", fraction=" << fraction << endl;
#endif

      //--- Scale covariance vector with fraction and add it to the total 
      cblas_daxpy(m,fraction,gamma_xIndv,1,gamma_x,1);

      //--- move pointer in param-array
      k += NparamIndv[i];
    }

    delete[] gamma_xIndv;
#ifdef DEBUG
    for (i=0;i<10;i++) {
      cout << "gamma_x, i=" << i << " : " << gamma_x[i] << endl;
    }
#endif
  }



/* Show the values of the parameters determined by Minimizer. The
 * first Nmodels-1 parameters are my phi's, the others are actual NoiseModel
 * parameters.
 */
//---------------------------------------------------
  void NoiseModel::show(double *param, double *error)
//---------------------------------------------------
  {
    int     i,j,k;
    double  fraction;

    using namespace std;
    k = Nmodels-1;
    for (i=0;i<Nmodels;i++) {
      
      //--- Compute fraction and show it
      fraction = 1.0;
      if (i<Nmodels-1) {
        for (j=0;j<Nmodels-1;j++) { 
          if (j==i) fraction *= (1.0-param[j]);
          else      fraction *= param[j];
        } 
      } else {
        for (j=0;j<Nmodels-1;j++) fraction *= param[j]; 
      }
      cout << endl << noisemodel[i] << ": fraction=" << fraction << endl;

      //--- Show value noise parameters
      modelIndv[i]->show(&param[k],&error[k]);

      //--- move pointer in param-array
      k += NparamIndv[i];
    }
  }



/* Each noise model has its own penalty function but in addition I don't let
 * my phi values go outside the [0:1] range.
 */
//-------------------------------------------------
  double NoiseModel::compute_penalty(double *param)
//-------------------------------------------------
  {
    int      i,k;
    double   penalty=0.0,LARGE=1.0e8;

    for (i=0;i<Nmodels-1;i++) {
      if (param[i]<0.0) {
        penalty += (0.0-param[i])*LARGE;
        param[i] = 0.0;
      } else if (param[i]>1.0) {
        penalty += (param[i]-1.0)*LARGE;
        param[i] = 1.0;
      }
    }
    k = Nmodels-1;
    for (i=0;i<Nmodels;i++) {
      penalty += modelIndv[i]->compute_penalty(&param[k]);

      //--- move pointer in param-array
      k += NparamIndv[i];
    }

    return penalty;
  }



/* To make a power spectral density plot, the values of param are needed.
 * In the future I could do something fancy by storing these values somewhere
 * on file and reading them automatically. Until then: manual input.
 */
//--------------------------------
  void NoiseModel::setup_PSD(void)
//--------------------------------
  {
    int    i;

    using namespace std;
    //--- Ask user about the fractions
    for (i=0;i<Nmodels;i++) {
      cout << "Enter fraction for model " << noisemodel[i] << ": ";
      cin >> fraction_fixed[i];
    }
    //--- Setup each noise model
    for (i=0;i<Nmodels;i++) {
      cout << endl << noisemodel[i] << ":" << endl;
      modelIndv[i]->setup_PSD();
    }
  }



/* For the power spectral density plot, created by the program ModelSpectrum,
 * I need to know for each model the value of G. These are summed after
 * scaling by their fraction.
 */
//-------------------------------------------
  double NoiseModel::compute_G(double lambda)
//-------------------------------------------
  {
    int      i,j;
    double   G=0.0;

    using namespace std;
    for (i=0;i<Nmodels;i++) {
      G += fraction_fixed[i]*modelIndv[i]->compute_G(lambda);
    }

    return G;
  }



/* To make a power spectral density plot, the values of param are needed.
 * In the future I could do something fancy by storing these values somewhere
 * on file and reading them automatically. Until then: manual input.
 */
//----------------------------------------
  void NoiseModel::setup_MonteCarlo(int m)
//----------------------------------------
  {
    int    i,j;

    using namespace std;
    //--- Ask user about the driving white noise and the fractions
    cout << "Enter sigma of driving white noise: ";
    cin >> sigma_fixed;
    for (i=0;i<Nmodels;i++) {
      cout << "Enter fraction for model " << noisemodel[i] << ": ";
      cin >> fraction_fixed[i];
    }
 
    //--- Allocate memory to hold impulse function and its Fourier transform
    //    for each noise model.
    h   = new double [2*m];
    w   = new double [2*m];
    F_h = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nmodels*(m+1));
    F_w = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(m+1));

    //--- we need to create a forward and backward plan
    plan_forward  = fftw_plan_dft_r2c_1d ( 2*m,   w, F_w, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d ( 2*m, F_w,   w, FFTW_ESTIMATE);

    //--- Setup each noise model and already compute the Fourier transform of h
    for (i=0;i<Nmodels;i++) {
      cout << endl << noisemodel[i] << ":" << endl;
      modelIndv[i]->compute_impulse_response(m,h);
      memset(&h[m],0,m*sizeof(double)); //--- zero padding of last m entries
      fftw_execute_dft_r2c(plan_forward,h,&F_h[i*(m+1)]);
#ifdef DEBUG
      for (j=0;j<m+1;j++) cout << "j=" << j << ", " << F_h[j][0] << ", "
					 	    << F_h[j][1] << endl;
#endif
    }
    cout << "Setup completed..." << endl;
  }



/*! Create a synthetic noise time-series using the given noise properties
 */
//-----------------------------------------------
  void NoiseModel::create_noise(int m, double *y)
//-----------------------------------------------
  {
    int     i,j,k;
    double  Scale,re,im;

    //--- Inverse FFT needs to be scaled
    Scale  = 1.0/static_cast<double>(2*m); 

    //--- Nicely start with zeros
    memset(y,0,m*sizeof(double));

    using namespace std;
    for (i=0;i<Nmodels;i++) {
      memset(&w[m],0,m*sizeof(double)); //--- zero padding of last m entries
      for (j=0;j<m;j++) w[j] = gsl_ran_gaussian(r_random,fraction_fixed[i]);
      fftw_execute_dft_r2c(plan_forward,w,F_w);

      //--- Compute F_h*F_w in frequency domain
      k = i*(m+1);
      for (j=0;j<m+1;j++) {
        re = F_h[k+j][0]*F_w[j][0]-F_h[k+j][1]*F_w[j][1];
        im = F_h[k+j][1]*F_w[j][0]+F_h[k+j][0]*F_w[j][1];
        F_w[j][0] = re;
        F_w[j][1] = im;
      }
 
      //--- Convert result of convolution back to time domain
      fftw_execute_dft_c2r(plan_backward,F_w,w);
      //for (j=0;j<2*m;j++) cout << "> j=" << j << ", " << w[j] << endl;

      //--- add to vector y
      cblas_daxpy(m,1.0,w,1,y,1);
      //for (j=0;j<m;j++) cout << "y --> j=" << j << ", " << w[j] << endl;
    }

    //--- Finally, apply white noise sigma
    cblas_dscal(m,sigma_fixed*Scale,y,1);
  }


  
/* The noise model is used so many times by other classes that it makes
 * sense to turn this into a singleton
 */
//-----------------------------------------
  NoiseModel* NoiseModel::getInstance(void)
//-----------------------------------------
  {
    if (instanceFlag==false) {
      instanceFlag=true;
      singleton = new NoiseModel();
    }
    return singleton;
  }
