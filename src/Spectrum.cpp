/*! \file    Spectrum.cpp
 *  \author  Machiel Bos
 * 
 *  This class contains some useful subroutines for computing the
 *  Welch power spectrum.
 *
 *  Reference:
 *  ----------
 *  Spectral Analysis and Filter Theory in Applied Geophysics by
 *  Burkhard Buttkus (2000), Springer.
 *
 *  This script is part of Hector 1.8
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
  #include <iostream>
  #include <ostream>
  #include <cmath>
  #include <cstring>
  #include <cstdlib>
  #include "Spectrum.h"
  #include "Control.h"

//  #define DEBUG
//  #define TIMING

//==============================================================================
// Subroutines
//==============================================================================


//--!!---------------------------------------------------------------------
  Spectrum::Spectrum(int n_, int segments_, double fs_) : pi(atan(1.0)*4.0)
//--!!---------------------------------------------------------------------
  {
    using namespace std;
    Control      &control=Control::getInstance();
    string       name;
    int          i,l;
    double       U,s;

    //--- Store values of n and segments in class
    n = n_;
    segments = segments_;
    fs = fs_;

    //--- Which window function needs to be applied?
    try {
      control.get_string("WindowFunction",name);
      if (name.compare("Parzen")==0) {
        windowfunction = &Spectrum::Parzen;
      } else if (name.compare("Hann")==0) {
        windowfunction = &Spectrum::Hann;
      } else {
        cerr << "Unknown WindowFunction: " << name << endl;
        exit(EXIT_FAILURE);
      }
      cout << "window function : " << name << endl;
    }
    catch (exception &e) {
      cout << e.what() << endl;
      cout << "Parzen window function is used." << endl;
      windowfunction = &Spectrum::Parzen;
    }
   
    //--- Fraction which determines how much the window function is applied    
    try {
      fraction = control.get_double("Fraction");
      if (fraction<0.0 || fraction>1.0) {
        cerr << "fraction must lie between 0 and 1, not:" << fraction << endl;
        exit(EXIT_FAILURE);
      }
      cout << "window function fraction: " << fraction << endl;
    }
    catch (exception &e) {
      cout << e.what() << endl;
      cout << "A fraction of 0.1 is used." << endl;
      fraction = 0.1; // first and last 10% is put through the window function
    }

    //--- Read data from file. Note that the Observation class fills the
    //    missing data with NaN's.
    dt   = 1.0/fs;

    L     = n/segments;
    K     = 2*segments-1;  // Total number of segments used in computation
                           // using half overlap.
    N     = L*segments;

    cout << "Number of data points n : " << n << endl;
    cout << "Number of data used   N : " << N << endl;
    cout << "Number of segments    K : " << K << endl;
    cout << "Length of segments    L : " << L << endl;

    //--- Compute threshold filter
    m = static_cast<int>(fraction*static_cast<double>(L/2));

    //--- Prepare FFT transformation
    try {
      Y     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (L/2+1));
      f     = new double[L/2+1];
      dummy = new double[L];
    } catch(bad_alloc) {
      cerr << "Not enough memory to hold all data...." << endl;
      exit(EXIT_FAILURE);
    }
    plan_forward = fftw_plan_dft_r2c_1d(L,dummy,Y,FFTW_ESTIMATE);

    //--- Compute frequencies
    for (i=0;i<L/2+1;i++) {
      f[i] = static_cast<double>(i)/(dt*static_cast<double>(L));
    }

    //--- Compute how much the window function alters the scale
    U = 0.0;
    for (i=0;i<L;i++) {
      //--- Compute fraction of window function between L/2-m and L/2
      //    I put it here for maximum speed although it is not very clear.
      l = i-L/2;       // array goes from 0..L window func from -L/2 to L/2
      if (l<0) l = -l; // symmetric window function, only look as l>0
      if (l>=L/2-m) {  // touced by filter
        if (m>0) {
          s = static_cast<double>(l-(L/2-m))/static_cast<double>(m);
        } else {
          s = 0.0;
        }
        //--- fraction s has been computed, call window function
        U += pow((*this.*windowfunction)(s),2.0);
      } else {
        //--- No need for any window function! Weight is simply 1
        U += 1.0;
      } 
    }
    U /= static_cast<double>(L); // follow equation 9.96
    scale = dt/(U*static_cast<double>(L));
    cout << "U : " << U << endl;
    cout << "dt: " << dt << endl;
    cout << "scale for G to get Amplitude (mm): " <<
                                        1000.0 * sqrt(2.0/(L*dt)) << endl;
  }



//--!!---------------------
  Spectrum::~Spectrum(void)
//--!!---------------------
  {
    delete[] dummy;
    delete[] f;
    fftw_destroy_plan (plan_forward);
    fftw_free (Y);
  }



/*! apply Parzen window function 
 *
 * \param[in] s  fraction between 0 and 1 (1 is end of filter)
 *
 * \returns{value of window function}
 */
//---------------------------------
  double Spectrum::Parzen(double s)
//---------------------------------
  {
    double  w;

    if (s<0.5) {
      w = 1.0 - 6.0*s*s + 6.0*s*s*s;
    } else {
      s = 1.0 - s;  // function is: 2.0*(1-s)^3
      w = 2.0*s*s*s;
    }
    return w;
  }



/*! apply Hann window function
 *
 * \param[in] s  fraction between 0 and 1 (1 is end of filter)
 *
 * \returns{value of window function}
 */
//-------------------------------
  double Spectrum::Hann(double s)
//-------------------------------
  {
    double          w;

    w = 0.5*(1.0 + cos(pi*s));
    
    return w;
  }



/*! To improve the PSD estimate, the mean must be removed.
 *
 * \param[in]   x     segment array with observations values
 * \param[in]   n     length of segment
 */
//--------------------------------------------
  void Spectrum::remove_mean(double *x, int n)
//--------------------------------------------
  {
    int     i,j;
    double  mean;

    mean=0.0;
    j   =0;
    for (i=0;i<n;i++) {
      if (!std::isnan(x[i])) {
        mean += x[i];
        j++;
      }
    }
    mean /= static_cast<double>(j);
    for (i=0;i<n;i++) if (!std::isnan(x[i])) x[i] -= mean;
  }
 

 
/*!Compute the spectrum using Welch averaging process and a window
 * function.
 *
 * \param[out] N        number of frequencies
 * \param[out] f        array filled with frenquencies (Hz)
 * \param[out] G        array filled with spectrum (mm^2/Hz)
 * \param[in]  segment  number of segments, without counting overlaps
 */
//--------------------------------------------
  double Spectrum::welch(double *y, double *G)
//--------------------------------------------
  {
    using namespace std;
    int             i,j,k,l,skipped_segments,ms1,ms2,ms3;
    double          Variance_Gf,Re,Im;
    double          percentage_gaps,s;
    clock_t         dtime1,dtime2,dsum2=0.0,dsum3=0.0;


    //--- For Andre Rodrigues, do some timing
    dtime1 = clock();

    for (i=0;i<=L/2;i++)  G[i] = 0.0;
    Variance_Gf = 0.0;
    skipped_segments = 0;
    for (j=0;j<K;j++) {

      dtime2 = clock();
      //--- Fill dummy array to process next segment
      for (k=0,i=0;i<L;i++) {
        if (std::isnan(y[j*L/2+i])) {
          dummy[i] = 0.0;
          k++;
        } else {

          //--- The following lines to compute fraction s were first included
          //    in the subroutines of Parzen and Hann. However, for maximum
          //    speed, I moved them here.
          l = i-L/2;
          if (l<0) l = -l;
          if (l>=L/2-m) { // touced by filter
            if (m>0) {
              s = static_cast<double>(l-(L/2-m))/static_cast<double>(m);
            } else {
              s = 0.0;
            }
            dummy[i] = y[j*L/2+i]* (*this.*windowfunction)(s);
       
          //--- No need to call window funtion. One function call avoided!
          } else {
            dummy[i] = y[j*L/2+i];
          }
        }
      }
      dsum3 += clock() - dtime2;

      percentage_gaps = 100.0*static_cast<double>(k)/static_cast<double>(L);

      // I begin to irritate me about the need to interpolate the data and
      // the horrible misfit between modelled spectrum and periodograms when
      // there are a lot of gaps. Now I try to skip segments which have too
      // many gaps which reduce the real power of the PSD.
      if (percentage_gaps>50.0) {
        skipped_segments++;
        cout << "Skipping segment " << j << ", gaps=" <<percentage_gaps<< endl;
      } else {
        dtime2 = clock();
        fftw_execute_dft_r2c(plan_forward,dummy,Y);
        dsum2 += clock() - dtime2;

        G[0] += scale*(Y[0][0]*Y[0][0] + Y[0][1]*Y[0][1]);
        if (L%2==1) {
          G[L/2] += scale*(Y[L/2][0]*Y[L/2][0] + Y[L/2][1]*Y[L/2][1]);
        } else {
          G[L/2] += 2.0*scale*(Y[L/2][0]*Y[L/2][0] + Y[L/2][1]*Y[L/2][1]);
        }
        for (i=1;i<(L/2);i++) {
          G[i] += 2.0*scale*(Y[i][0]*Y[i][0] + Y[i][1]*Y[i][1]);
        }
      }
    }

    //--- Divide total to get averaged PSD, which is more accurate
    Variance_Gf = 0.0;
    for (i=0;i<L/2+1;i++) {
      if (skipped_segments<K) {
        G[i] /= static_cast<double>(K-skipped_segments);
      } else {
        G[i]  = 0.0;
      }
      Variance_Gf+=G[i]/(dt*static_cast<double>(L));
    }
  
    dtime1 = clock() - dtime1;
    ms1 = (float(dtime1))/ CLOCKS_PER_SEC  * 1000;
    ms2 = (float(dsum2))/ CLOCKS_PER_SEC  * 1000;
    ms3 = (float(dsum3))/ CLOCKS_PER_SEC  * 1000;
#ifdef TIMING
    cout << "** Timing - Welch     : " << ms1 << endl;
    cout << "** Timing - FFT       : " << ms2 << endl;
    cout << "** Timing - Window    : " << ms3 << endl;
#endif

    //--- Done! just return variance_Gf as bonus
    return Variance_Gf;
  }



//-----------------------------------
  void Spectrum::get_freq(double *f_)
//-----------------------------------
  {
    int  i;
   
    for (i=0;i<L/2+1;i++) f_[i] = f[i];
  }
