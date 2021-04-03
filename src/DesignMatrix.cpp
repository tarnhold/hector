/*! \file   DesignMatrix.cpp
 *  \author Machiel Bos
 *
 * This class provides the Design Matrix H that is needed in the 
 * Least-Squares part of the computation of the Likelihood.
 *
 * \date 13/1/2012   Coimbra Library
 */
//==============================================================================
  #include "DesignMatrix.h"
  #include <iostream>
  #include <ostream>
  #include <iomanip>
  #include <cmath>
  #include <cstdlib>

  #define TINY 1.0e-5
//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================

//++++++ Singleton stuff +++++++++++++++++++++
bool DesignMatrix::instanceFlag = false;
DesignMatrix* DesignMatrix::singleton = NULL;
//++++++++++++++++++++++++++++++++++++++++++++


//---!!------------------------------------------------
  DesignMatrix::DesignMatrix(void) : tpi(8.0*atan(1.0))
//---!!------------------------------------------------
  {
    using namespace std;
    bool            estimate_offsets;
    int             i,j,k,l;
    double          *t,*x;
    char            **periodic_signals;
    Observations    *observations=Observations::getInstance();
    Control         *control=Control::getInstance();

    control->get_string("PhysicalUnit",unit);
    first_difference = control->get_bool("firstdifference");
    //--- determine size of matrices
    if (first_difference==true) {
      n = 1;
    } else {
      n = 2; // offset + trend
    }

    //--- Do we need to include seasonal signal
    seasonal_signal      = control->get_bool("seasonalsignal");
    halfseasonal_signal  = control->get_bool("halfseasonalsignal");
    if (seasonal_signal==true)     n +=2;
    if (halfseasonal_signal==true) n +=2;

    //--- An additional number of 20 periodic signals may be added
    periodic_signals = new char* [20];
    control->get_name_list("periodicsignals",periodic_signals,
							n_periodic_signals);
    if (n_periodic_signals>20) {
      cerr << "Only up to 20 periodic signals are allowed!" << endl;
      exit(EXIT_FAILURE);
    }
    n += 2*n_periodic_signals;
    if (n_periodic_signals>0) {
      periods = new double[n_periodic_signals];
      for (i=0;i<n_periodic_signals;i++) {
        sscanf(periodic_signals[i],"%lf",&periods[i]);
      }
    } else {
      periods = NULL;
    }

    observations->get_values(m,&t,&x);
    Ngaps = observations->number_of_gaps();
    dt = 1.0/(24.0*3600.0*observations->get_fs()); // unit is MJD days!!!
#ifdef DEBUG
    cout << "m=" << m << endl;
    cout << "dt=" << dt << endl;
#endif

    //--- Maybe we need to include offsets too
    estimate_offsets = control->get_bool("estimateoffsets");
    if (first_difference==false && estimate_offsets==true) {
      observations->get_offsets(offsets);
      n_offsets  = offsets.size();
      n         += n_offsets;
    } else {
      n_offsets  = 0;
    }

    //--- Allocate memory
    try {
      H = new double[m*n];
      if (Ngaps>0)  F = new double[Ngaps*m];
      else          F = NULL;
    }
    catch(bad_alloc) {
      cout << "Need more memory to store Design Matrix" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Construct the design matrix H
    th     = 0.5*(t[0] + t[m-1]);
    for (i=0;i<m;i++) {
      j = 1;
      if (first_difference==false) {
        H[i] = 1.0;
        H[i+j*m] = t[i]-th;                                     j++;
        if (seasonal_signal==true) {
          H[i+j*m] = cos(tpi*(t[i]-51544.0)/365.25);            j++;
          H[i+j*m] = sin(tpi*(t[i]-51544.0)/365.25);            j++;
        }
        if (halfseasonal_signal==true) {
          H[i+j*m] = cos(2.0*tpi*(t[i]-51544.0)/365.25);        j++;
          H[i+j*m] = sin(2.0*tpi*(t[i]-51544.0)/365.25);        j++;
        }
        for (k=0;k<n_periodic_signals;k++) {
          H[i+j*m] = cos(tpi*(t[i]-51544.0)/periods[k]);        j++;
          H[i+j*m] = sin(tpi*(t[i]-51544.0)/periods[k]);        j++;
        }
        for (l=0;l<n_offsets;l++) {
          for (k=0;k<m;k++) {
            if (t[k]+TINY>offsets[l]) H[k + (j+l)*m] = 1.0;
          }
        }
        j += n_offsets;
      } else {
        //--- Now treat the first differenced data case
        H[i] = dt;
        if (seasonal_signal==true) {
          H[i+j*m] = cos(tpi*(t[i]-51544.0)/365.25) -
                          cos(tpi*((t[i]-dt)-51544.0)/365.25);  j++;
          H[i+j*m] = sin(tpi*(t[i]-51544.0)/365.25) -
                          sin(tpi*((t[i]-dt)-51544.0)/365.25);  j++;
        }
        if (halfseasonal_signal==true) {
          H[i+j*m] = cos(2.0*tpi*(t[i]-51544.0)/365.25) -
                          cos(2.0*tpi*((t[i]-dt)-51544.0)/365.25);  j++;
          H[i+j*m] = sin(2.0*tpi*(t[i]-51544.0)/365.25) -
                          sin(2.0*tpi*((t[i]-dt)-51544.0)/365.25);  j++;
        }
        for (k=0;k<n_periodic_signals;k++) {
          H[i+j*m] = cos(tpi*(t[i]-51544.0)/periods[k]) -
			  cos(tpi*((t[i]-dt)-51544.0)/periods[k]);  j++;
          H[i+j*m] = sin(tpi*(t[i]-51544.0)/periods[k]) -
			  sin(tpi*((t[i]-dt)-51544.0)/periods[k]);  j++;
        }
      }
    }

    //--- Construct selection matrix F
    if (Ngaps>0) {
      //--- Make zero
      //cblas_dscal(Ngaps*m,0.0,F,1); !!! Does not work NaN*0.0=NAN !!!!
      memset(F,0,Ngaps*m*sizeof(double));

      //--- Put the ones in the right places
      j=0;
      for (i=0;i<m;i++) {
        if (isnan(x[i])) {
          F[i + j*m] = 1.0;
          j++;
        }
      }
    } else {
      F = NULL;
    }

    //--- free memory
    for (i=0;i<n_periodic_signals;i++) delete[] periodic_signals[i];
    delete[] periodic_signals;
  }



/*! Free up memory
 */
//---!!----------------------------
  DesignMatrix::~DesignMatrix(void)
//---!!----------------------------
  {
    if (H!=NULL)        delete[] H;
    if (F!=NULL)        delete[] F;
    if (periods!=NULL)  delete[] periods;
  }



/*! Copy designmatrix H to a new place in memory and return its pointer
 *
 * \param[out]   n_   number of columns
 * \param[out]   H_   design matrix
 */
//----------------------------------------------
  void DesignMatrix::get_H(int& n_, double **H_)
//----------------------------------------------
  {
    n_  = n;
    *H_ = new double[n*m];
    cblas_dcopy(n*m,H,1,*H_,1);
  }



/*! For my gap stuff, I need a special selection matrix 'F'
 *
 * \param[out]  F_   my special matrix to deal with missing data
 */
//-------------------------------------
  void DesignMatrix::get_F(double **F_)
//-------------------------------------
  {
    using namespace std;
    if (Ngaps>0) {
      *F_ = new double[Ngaps*m];
      cblas_dcopy(Ngaps*m,F,1,*F_,1);
    } else {
      cerr << "DesignMatrix: Warning! asking for F but Ngaps=0!" << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! show meaning of estimated least-squares parameters on screen
 *
 * \param[in]  theta    vector containing estimated parameters of LS
 * \param[in]  error    vector with error of values in vector theta
 */
//-------------------------------------------------------------
  void DesignMatrix::show_results(double *theta, double *error)
//-------------------------------------------------------------
  {
    int      i,j;
    double   ds = 365.25;

    using namespace std;
    i=0;
    cout << fixed << setprecision(6);
    if (first_difference==false) {
      cout << "bias : " << theta[0] << " +/- " << error[0] << " " << unit 
           << " (at MJD=" << th << ")" << endl;
      cout << "trend: " << theta[1]*ds << " +/- " << error[1]*ds
           << " " << unit << "/year" << endl;
      i += 2;
    } else {
      cout << "trend: " << theta[0]*ds << " +/- " << error[0]*ds
           << " " << unit << "/year" << endl;
      i += 1;
    }

    if (seasonal_signal==true) {
      cout << "cos yearly : " << theta[i] <<" +/- "<<error[i]<<" "<<unit<<endl;
      i++;
      cout << "sin yearly : " << theta[i] <<" +/- "<<error[i]<<" "<<unit<<endl;
      i++;
    }

    if (halfseasonal_signal==true) {
      cout << "cos hyearly : " << theta[i] <<" +/- "<<error[i]<<" "<<unit<<endl;
      i++;
      cout << "sin hyearly : " << theta[i] <<" +/- "<<error[i]<<" "<<unit<<endl;
      i++;
    }
    for (j=0;j<n_periodic_signals;j++) {
      printf("cos %7.2lf : %lf +/- %lf %s\n",periods[j],theta[i],
							error[i],unit); i++;
      printf("sin %7.2lf : %lf +/- %lf %s\n",periods[j],theta[i],
							error[i],unit); i++;
    }
    for (j=0;j<n_offsets;j++) {
      printf("offset at %10.4lf : %lf +/- %lf %s\n",offsets[j],
					       theta[i],error[i],unit); i++;
    }
  }



/*! Compute xhat
 *
 * \param[in]  theta   vector with estimated parameters from LS
 */
//----------------------------------------------------
  void DesignMatrix::compute_xhat(const double *theta)
//----------------------------------------------------
  {
    double        xhat[m];
    Observations  *observations=Observations::getInstance();
   
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
					     m,1,n,1.0,H,m,theta,n,0.0,xhat,m);
    observations->set_xhat(xhat);
  }



//--!!-----------------------------------------
  DesignMatrix* DesignMatrix::getInstance(void)
//--!!-----------------------------------------
  {
    if (instanceFlag==false) {
      singleton = new DesignMatrix();
      instanceFlag=true;
    } 
    return singleton;
  }



/*! For this singleton, the instance can be manually destroyed */
//----------------------------------------
  void DesignMatrix::destroyInstance(void)
//----------------------------------------
  {
    if (instanceFlag==true) {
      delete singleton;
      instanceFlag = false;
    }
  }

