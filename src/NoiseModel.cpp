/*! \file    NoiseModel.cpp
 *  \author  Machiel Bos
 *
 * Main interface to various Noise model implementations. It allows the
 * combination of various noise models.
 *
 * \date 13/ 1/2012  Santa Clara
 * \date  7/10/2012  Santa Clara
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
    control->get_name_list("NoiseModels",noisemodel,Nmodels);
    if (Nmodels<1) {
      cerr << "No noise model was specified!" << endl;
      exit(EXIT_FAILURE);
    } else if (Nmodels>5) {
      cerr << "Only up to 5 noise models can be specified!" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Allocate memory space to store each individual noise model
    Nparam       = Nmodels-1;
    modelIndv    = new NoiseModelBaseClass* [Nmodels];
    NparamIndv   = new int[Nmodels];
    phi          = new double[Nparam];
    fraction_PSD = new double[Nparam];
    
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

    if (Nmodels>0) {
      for (i=0;i<Nmodels;i++) {
        delete[] noisemodel;
      }
      delete[] modelIndv;
      delete[] NparamIndv;
      delete[] phi;
      delete[] fraction_PSD;
    }
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
      cin >> fraction_PSD[i];
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
      G += fraction_PSD[i]*modelIndv[i]->compute_G(lambda);
    }

    return G;
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
