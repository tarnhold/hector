/*! \file    EstimateSpectrum.cpp
 *  \author  Machiel Bos
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
  #include <ctime>
  #include "Observations.h"
  #include "Spectrum.h"

//  #define DEBUG
//  #define TIMING

//==============================================================================
// Subroutines
//==============================================================================


 
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
    clock_t         dtime;
    int             ms;
    
    dtime = clock();
    Observations    &data=Observations::getInstance();
    dtime = clock() - dtime;
    ms = (float(dtime))/ CLOCKS_PER_SEC  * 1000;
#ifdef TIMING
    cout << "** Timing - Read file : " << ms << endl;
#endif
    Spectrum        *spectrum;
    int             i,j,L,n;
    double          Variance_xt,Variance_Gf,dt,fs;
    double          *t,*y,freq[2];
    fstream         fp;
    FILE            *fp_res;

    //--- Read data from file. Note that the Observation class fills the
    //    missing data with NaN's.
    fs   = data.get_fs();
    dt   = 1.0/fs;
    data.get_values(n,&t,&y); // only arrays t and y are needed from now on.
    L    = n/segments;
    N    = L/2+1;

    //--- Initialise Spectrum class
    spectrum = new Spectrum(n,segments,fs);

    try {
      *G = new double[N]; // spectrum (periodogram) amplitude
      *f = new double[N]; // frequencies, used for plotting
    } catch(bad_alloc) {
      cerr << "Not enough memory to hold all data...." << endl;
      exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    fp_res = fopen("residuals.dat","w");
    for (i=0;i<n;i++) fprintf(fp_res,"%e\n",y[i]);
    fclose(fp_res);
#endif

    Variance_xt = 0.0; // Total variance, page 167, SAaFT
    for (j=0,i=0;i<n;i++) {
      if (!std::isnan(y[i])) { // Skip NaN's
        Variance_xt += y[i]*y[i];
        j++;
      }
    }
    //--- Avoid catastrophes
    if (j>1) Variance_xt /= static_cast<double>(j-1);
    else     Variance_xt  = 0.0;

    cout << "Total variance in signal (time domain): " << Variance_xt << endl;

    //--- Use Welch method to compute power-spectrum
    Variance_Gf = spectrum->welch(y,*G);

    //--- Get frequencies
    spectrum->get_freq(*f);

    freq[0] = (*f)[1];
    freq[1] = (*f)[L/2];
    cout << "Total variance in signal (spectrum)   : " << Variance_Gf << endl;
    cout << "freq0: " << scientific << freq[0] << endl;
    cout << "freq1: " << scientific << freq[1] << endl;
    //--- note: energy under plot = \int_0^{f_N}G df ~ G/(2*dt*L/2)
    //          Therefore, variance computed above takes into account
    //          negative frequencies. I don't need those so multiply by 2

    //--- Clean up mess
    delete[] y;
    delete[] t;
    delete   spectrum;
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
    int       i,N,ms;
    clock_t   dt;


    using namespace std;
    compute_periodogram(N,&f,&G,segments);

    dt = clock();
    write_periodogram(N,f,G);
    dt = clock() -dt;
    ms = (float(dt))/ CLOCKS_PER_SEC  * 1000;
#ifdef TIMING
    cout << "** Timing - Write file: " << ms << endl;
#endif


     
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
