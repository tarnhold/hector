/*! \file    EstimateSpectrum.cpp
 *  \author  Machiel Bos
 *  \version 1.0
 * 
 *  A program to estimate the Power Spectral Density of a data set. The data 
 *  set consists out of two columns (MJD + data) or three (MJD, data and 
 *  model). If there are three columns, the PSD is computed for the residuals 
 *  (data-model). It uses the Welch periodogram method and Parzen window 
 *  function. The equation numbers and book pages refer to:
 *
 *  Reference:
 *  ----------
 *  Spectral Analysis and Filter Theory in Applied Geophysics by
 *  Burkhard Buttkus (2000), Springer.
 *
 * \date 16/ 9/2004  Vila Nova de Gaia
 * \date 13/ 7/2005  Vila Nova de Gaia
 * \date  6/ 7/2009  CIIMAR, Porto
 * \date 25/ 7/2011  Coimbra
 * \date 12/11/2012  Santa Clara
 */
//==============================================================================
  #include <iostream>
  #include <iomanip>
  #include <ostream>
  #include <istream>
  #include <fstream>
  #include <cmath>
  #include <vector>
  #include <cstring>
  #include <string>
  #include <cstdlib>
  #include <cstdlib>
  #include <fftw3.h>
  #include "Observations.h"

//  #define DEBUG

//==============================================================================
// Global variables
//==============================================================================

    char           filename[20];
    double         dt;


//==============================================================================
// Subroutines
//==============================================================================


/*! apply Parzen window function [-M..M]
 * fraction is how much of the data is used in the window function. Buttus
 * says 10% is normal. Thus Parzen filter applies [-M,-0.9M] and [0.9M,M]
 * with 1 in between and zeros outside.
 *
 * \param[in] i          index between 0 and M (location)
 * \param[in] M          half-width of window (-M:M)
 * \param[in] fraction   the fraction inside the window which does not alter
 *                       the data. Only 1-fraction gets altered.
 * \returns{value of window function}
 */
//--------------------------------------------
  double Parzen(int i, int M, double fraction)
//--------------------------------------------
  {
    int             m;
    double          w;
    static double   pi = atan(1.0)*4.0,s;

    if (i>=M) return 0.0; // outside range always zero
    if (i<0) i = -i;      // only symmetric window functions treated here
    m = static_cast<int>(fraction*static_cast<double>(M));
    if (i<M-m) { // untouced by filter
      return 1.0;
    } else { // filter applies
      if (m>0) {
        //--- i ranges from M-m to M
        s = static_cast<double>(i-(M-m))/static_cast<double>(m);
      } else {
        s = 0;
      }
    }
    if (i<M/2) {
      w = 1.0 - 6.0*s*s + 6.0*s*s*s;
    } else {
      s = 1.0 - s;  // function is: 2.0*(1-s)^3
      w = 2.0*s*s*s;
    }
    return w;
  }



/*! apply Hann window function [-M..M]
 * fraction is how much of the data is used in the window function. Buttus
 * says 10% is normal. Thus Parzen filter applies [-M,-0.9M] and [0.9M,M]
 * with 1 in between and zeros outside.
 *
 * \param[in] i          index between 0 and M (location)
 * \param[in] M          half-width of window (-M:M)
 * \param[in] fraction   the fraction inside the window which does not alter
 *                       the data. Only 1-fraction gets altered.
 * \returns{value of window function}
 */
//------------------------------------------
  double Hann(int i, int M, double fraction)
//------------------------------------------
  {
    int             m;
    double          w;
    static double   pi = atan(1.0)*4.0,s;

    if (i>=M) return 0.0; // outside range always zero
    if (i<0) i = -i;      // only symmetric window functions treated here
    m = static_cast<int>(fraction*static_cast<double>(M));
    if (i<M-m) { // untouced by filter
      return 1.0;
    } else { // filter applies
      if (m>0) {
        //--- i ranges from M-m to M
        s = static_cast<double>(i-(M-m))/static_cast<double>(m);
      } else {
        s = 0;
      }
    }
    w = 0.5*(1.0 + cos(pi*s));
    
    return w;
  }



/*! To improve the PSD estimate, the mean must be removed.
 *
 * \param[in]   x     segment array with observations values
 * \param[in]   n     length of segment
 */
//----------------------------------
  void remove_mean(double *x, int n)
//----------------------------------
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
//----------------------------------------------------------------------
  void compute_periodogram(int& N, double **f, double **G, int segments)
//----------------------------------------------------------------------
  {
    using namespace std;
    Observations    &data=Observations::getInstance();
    int             i,j,L,K,n,k,skipped_segments;
    double          Variance_xt,Variance_Gf,U,scale,Re,Im,fraction;
    double          *t,*y,*dummy,MJD,obs,mod,dt,freq[2],percentage_gaps;
    double const    deg = 45.0/atan(1.0);
    char            line[80];
    fstream         fp;
    FILE            *fp_res;
    string          name;
    fftw_plan       plan_forward;
    fftw_complex    *Y=NULL;
    double         (*windowfunction)(int, int, double);
    Control        &control=Control::getInstance();

    //--- Which window function needs to be applied?
    try {
      control.get_string("WindowFunction",name);
      if (name.compare("Parzen")==0) {
        windowfunction = Parzen; 
      } else if (name.compare("Hann")==0) {
        windowfunction = Hann;
      } else {
        cerr << "Unknown WindowFunction: " << name << endl;
        exit(EXIT_FAILURE);
      }
      cout << "window function : " << name << endl;
    }
    catch (exception &e) {
      cout << e.what() << endl;
      cout << "Parzen window function is used." << endl;
      windowfunction = Parzen; 
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
    dt   = 1.0/data.get_fs();
    data.get_values(n,&t,&y); // only arrays t and y are needed from now on.

    //--- First, compute the segment length: L= n/segments
    //    FFTW does not require a length that is  a power of 2
    L     = n/segments; 
    K     = 2*segments-1;  // Total number of segments used in computation
                           // using half overlap.
    N     = L*segments;    
    dummy = new double[L]; // Its segment is put through window filter and put
                           // into the 'dummy' array and afterwards FFT-ed.

    cout << "Number of data points n : " << N << endl;
    cout << "Number of data used   N : " << N << endl;
    cout << "Number of segments    K : " << K << endl;
    cout << "Length of segments    L : " << L << endl;

    //--- Prepare FFT transformation
    try {
      Y  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (L/2+1));
      *G = new double[L/2+1]; // spectrum (periodogram) amplitude
      *f = new double[L/2+1]; // frequencies, used for plotting
    } catch(bad_alloc) {
      cerr << "Not enough memory to hold all data...." << endl;
      exit(EXIT_FAILURE);
    }
    plan_forward = fftw_plan_dft_r2c_1d(L,dummy,Y,FFTW_ESTIMATE);

    //--- Compute how much the window function alters the scale
    U = 0.0; 
    for (i=0;i<L;i++) {
      U += pow((*windowfunction)(i-L/2,L/2,fraction),2.0);
    }
    U /= static_cast<double>(L); // follow equation 9.96
    scale = dt/(U*static_cast<double>(L));
    cout << "U : " << U << endl;
    cout << "dt: " << dt << endl;
    cout << "scale for G to get Amplitude (mm): " << 
 					1000.0 * sqrt(2.0/(L*dt)) << endl;

    //--- Remove Mean to improve results at low frequencies
    //remove_mean(y,N);

#ifdef DEBUG
    fp_res = fopen("residuals.dat","w");
    for (i=0;i<N;i++) fprintf(fp_res,"%e\n",y[i]);
    fclose(fp_res);
#endif

    Variance_xt = 0.0; // Total variance, page 167, SAaFT
    for (j=0,i=0;i<N;i++) {
      if (!isnan(y[i])) { // Skip NaN's
        Variance_xt += y[i]*y[i];
        j++;
      }
    }
    //--- Avoid catastrophes
    if (j>1) Variance_xt /= static_cast<double>(j-1);
    else     Variance_xt  = 0.0;

    cout << "Total variance in signal (time domain): " << Variance_xt << endl;
    for (i=0;i<=L/2;i++)  (*G)[i] = 0.0;
    Variance_Gf = 0.0;
    skipped_segments = 0;
    for (j=0;j<K;j++) {

      //--- Fill dummy array to process next segment
      for (k=0,i=0;i<L;i++) {
        if (isnan(y[j*L/2+i])) {
          dummy[i] = 0.0;
          k++;
        } else {
          dummy[i] = y[j*L/2+i]* (*windowfunction)(i-L/2,L/2,fraction);
        }
      }

      percentage_gaps = 100.0*static_cast<double>(k)/static_cast<double>(L);

      // I begin to irritate me about the need to interpolate the data and
      // the horrible misfit between modelled spectrum and periodograms when
      // there are a lot of gaps. Now I try to skip segments which have too
      // many gaps which reduce the real power of the PSD.
      if (percentage_gaps>50.0) {
        skipped_segments++;
        cout << "Skipping segment " << j << ", gaps=" <<percentage_gaps<< endl;
      } else {
        //--- Remove mean of segment
        // remove_mean(dummy,L);
        fftw_execute_dft_r2c(plan_forward,dummy,Y);
        for (i=0;i<(L/2+1);i++) {
          (*G)[i] += 2.0*scale*(Y[i][0]*Y[i][0] + Y[i][1]*Y[i][1]);
        }
      }
    }

    //--- Divide total to get averaged PSD, which is more accurate
    Variance_Gf = 0.0;
    for (i=0;i<L/2+1;i++) {
      if (skipped_segments<K) {
        (*G)[i] /= static_cast<double>(K-skipped_segments);
      } else {
        (*G)[i]  = 0.0;
      }
      if (i==0 || i==L/2)  Variance_Gf+=(*G)[i]/(2.0*dt*static_cast<double>(L));
      else                 Variance_Gf+=(*G)[i]/(dt*static_cast<double>(L));

      (*f)[i] = static_cast<double>(i)/(dt*static_cast<double>(L));
    }
    freq[0] = (*f)[1];
    freq[1] = (*f)[L/2];
    cout << "Total variance in signal (spectrum)   : " << Variance_Gf << endl;
    cout << "freq0: " << scientific << freq[0] << endl;
    cout << "freq1: " << scientific << freq[1] << endl;
    //--- note: energy under plot = \int_0^{f_N}G df ~ G/(2*dt*L/2)
    //          Therefore, variance computed above takes into account
    //          negative frequencies. I don't need those so multiply by 2

    N = L/2+1;

    //--- Clean up mess
    delete[] y;
    delete[] t;
    delete[] dummy;
    fftw_destroy_plan (plan_forward);
    fftw_free (Y);
  }



/*!Compute the spectrum using Welch averaging process and Parzen window
 * function.
 *
 * \param[in]  N        number of frequencies
 * \param[in]  f        array filled with frenquencies (Hz)
 * \param[in]  G        array filled with spectrum (mm^2/Hz)
 */
//---------------------------------------------------
  void write_periodogram(int N, double *f, double *G)
//---------------------------------------------------
  {
    using namespace std;
    fstream       fp;
    int           i;
    string        filename;
    Control       &control=Control::getInstance();

    try {
      control.get_string("OutputFile",filename);
    } catch (exception &e) {
      filename = "estimatespectrum.out";
    }
    cout << "--> " << filename << endl;
    fp.open(filename.c_str(),ios::out);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    //--- zero frequency is normally skipped.
    for (i=1;i<N;i++) fp << scientific << f[i] << "   " << G[i] << endl;
    fp.close();
  }



/*! Make a file that contains the Welch periodogram for the given station name
 *
 * \param[in] segments : number of divisions of the time-series, without
 *                       counting the overlaps (50%).
 */
//-----------------------------------------
  void make_Welch_periodogram(int segments)
//-----------------------------------------
  {
    double    *G=NULL,*f=NULL;
    int       i,N;

    compute_periodogram(N,&f,&G,segments);
    write_periodogram(N,f,G);
    delete[] f;
    delete[] G;
  }



/*! Check if character array contains a integer
 */
//-------------------------
  bool IsInteger(char *str)
//-------------------------
  {
    int   i=0;

    using namespace std;
    while (str[i]!='\0') { 
      if (isdigit(str[i])==false) return false;
      i++;
    }
    if (i>0) return true; else return false;
  }

       

//==============================================================================
// Main program
//==============================================================================

  int main(int argc, char *argv[])
  {
    int       segments=4;

    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      Control &control = Control::getInstance("estimatespectrum.ctl");
    } else {
      if (IsInteger(argv[1])==true) {
        segments = atoi(argv[1]);
        if (segments<=0) {
          cerr << "Unacceptable number of segments: " << segments << endl;
          exit(EXIT_FAILURE);
        }
        if (argc==2)      
		Control &control= Control::getInstance("estimatespectrum.ctl");
        else if (argc==3) 
                Control &control= Control::getInstance(argv[2]); 
        else {
          cerr << "estimatespectrum [number of segments (not counting "
               << "overlaps)] [controlfile.ctl]" << endl; 
          exit(EXIT_FAILURE);
        }
      } else {
        if (argc==2) Control &control = Control::getInstance(argv[1]); 
        else {
          cerr << "estimatespectrum [number of segments (not counting "
               << "overlaps)] [controlfile.ctl]" << endl; 
          exit(EXIT_FAILURE);
        }
      }
    }

    //--- Start estimatespectrum
    cout << endl
         << "************************************" << endl
         << "    estimatespectrum, version " << VERSION << endl
         << "************************************" << endl;

    make_Welch_periodogram(segments);

    return EXIT_SUCCESS;
  }
