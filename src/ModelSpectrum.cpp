/*! \file    ModelSpectrum.cpp
 *  \author  Machiel Bos
 *
 * This program computes the modelled spectrum for given noise parameters.
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
  #include <cmath>
  #include <cstdlib>
  #include <cstdio>
  #include <iostream>
  #include <ostream>
  #include <vector>
  #include "ModelSpectrum.h"

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================

//---------------------------------------------------
  static int compare (const void * a, const void * b)
//---------------------------------------------------
  {
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;  
  }



//---!!--------------------------------------------------
  ModelSpectrum::ModelSpectrum(void) : pi(4.0*atan(1.0)),
                                       tpi(2.0*pi)
//---!!--------------------------------------------------
  {
    using namespace std;
    int               i,j;
    double            T,*theta;
    complex<double>   a,b,c,w;
    Control           &control=Control::getInstance();
    NoiseModel        &noisemodel = NoiseModel::getInstance();

    //--- Check for Monte Carlo derived confidence intervals
    try {
      montecarlo  = control.get_bool("MonteCarloConfidence");
    }
    catch (exception &e) {
      montecarlo = false;
    }

    //--- Do we need Monte Carlo?
    if (montecarlo==true) {
      try {
        n_simulations = control.get_int("NumberOfSimulations");
        m             = control.get_int("NumberOfPoints");
        segments      = control.get_int("NumberOfSegments");
      }
      catch (exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
      }
    }

    //--- Need to know the standard deviation of epsilon (the white noise)
    cout << "Enter the standard deviation of the driving noise: ";
    cin >> sigma;
    noisemodel.set_sigma_fixed(sigma);

    //--- For correct plots, I need to know the sampling period
    cout << "Enter the sampling period in hours: ";
    cin >> T;
    fs = 1.0/(T*3600.0);

    //--- make sure area under PSD figure is equal to variance
    scale = 2.0*pow(sigma,2.0)/fs; //--- no negative frequencies (2x)
  }





//---------------------------------
  void ModelSpectrum::compute(void)
//---------------------------------
  {
    using namespace std;
    const int   N=500;
    int         i,choice,Nparam;
    double      dlambda,lambda,*G,*f,s,freq[2],*params_fixed=NULL;
    fstream     fp;
    string      filename;
    Control     &control=Control::getInstance();
    NoiseModel  &noisemodel = NoiseModel::getInstance();

    //--- Setup Noisemodel parameters
    Nparam = noisemodel.get_Nparam();
    if (Nparam>0) params_fixed = new double[Nparam];
    noisemodel.set_noise_parameters(params_fixed);

    //--- create space to hold PSD
    G = new double[N+1];
    f = new double[N+1];

    cout << "1) Linear or 2) logarithmic scaling of frequency?: ";
    cin >> choice;
    if (choice==1) {
      for (i=0;i<=N;i++) {
        f[i] = 0.5*static_cast<double>(i)/static_cast<double>(N)*fs;
        G[i] = scale*noisemodel.compute_G(tpi*f[i]/fs);
      }
    } else if (choice==2) {
      cout << "Enter freq0 and freq1: ";
      cin >> freq[0] >> freq[1];
      cout << "freq0 : " << freq[0] << endl;
      cout << "freq1 : " << freq[1] << endl;
      freq[0] = log(freq[0]);
      freq[1] = log(freq[1]);
      for (i=0;i<=N;i++) {
        s    = static_cast<double>(i)/static_cast<double>(N);
        f[i] = exp((1.0-s)*freq[0] + s*freq[1]);
        G[i] = scale*noisemodel.compute_G(tpi*f[i]/fs);
      }
    } else {
      cerr << "Unknown choice" << endl;
      exit(EXIT_FAILURE);
    }

    //--- write results to file
    try {
      control.get_string("OutputFile",filename);
    } catch (exception &e) {
      filename = "modelspectrum.out";
    }
    cout << "--> " << filename << endl;
    fp.open(filename.c_str(),ios::out);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    for (i=0;i<=N;i++) {
      fp << scientific << f[i] << "   " << G[i] << endl;
    }
    fp.close();

    //--- Do we need confidence levels?
    if (montecarlo==true) {
      noisemodel.setup_MonteCarlo(m);
      compute_confidence_levels();
    }

    //--- Clean up mess
    delete[] G;
    delete[] f;
    if (params_fixed!=NULL) delete[] params_fixed;
  }
    


//---------------------------------------------------
  void ModelSpectrum::compute_confidence_levels(void)
//---------------------------------------------------
  {
    const int   N=500;
    int         n2,i,ms,j,L;
    double      *f,*G,*y,*G_5,*G_50,*G_95,*dummy,f_int;
    double      G_5_int,G_50_int,G_95_int,s,f0,f1;
    Spectrum    spectrum(m,segments,fs);
    NoiseModel  &noisemodel=NoiseModel::getInstance();
    clock_t     dt;
    FILE        *fp;

    using namespace std;
    //--- number of frequencies in the FFT
    L  = m/segments;
    n2 = L/2 + 1;
   
    //--- Allocate memory
    try {
      y    = new double[m];
      f    = new double[n2];
      G    = new double[n2*n_simulations];
      G_5  = new double[n2];
      G_50 = new double[n2];
      G_95 = new double[n2];
      dummy= new double[n_simulations];
    }
    catch (bad_alloc()) {
      cerr << "Need more memory for f and G...." << endl;
      exit(EXIT_FAILURE);
    }

    dt = clock();
    //--- create synthetic noise and run welch over it
    for (i=0;i<n_simulations;i++) {
      //--- create noise
      noisemodel.create_noise(m,y);
      spectrum.welch(y,&G[i*n2]);
    }
    dt = clock() - dt;
    ms = (float(dt))/ CLOCKS_PER_SEC  * 1000;
    cout << "** Timing - All simul : " << ms << endl;

    //--- Compute mean and std of G
    dt = clock();
    for (j=0;j<n2;j++) {
      cblas_dcopy(n_simulations,&G[j],n2,dummy,1);
      qsort(dummy,n_simulations,sizeof(double),compare);
      G_5[j]  = dummy[ 5*(n_simulations/100)]; 
      G_50[j] = dummy[50*(n_simulations/100)]; 
      G_95[j] = dummy[95*(n_simulations/100)]; 
    }
    dt = clock() - dt;
    ms = (float(dt))/ CLOCKS_PER_SEC  * 1000;
    cout << "** Timing - percentiles  : " << ms << endl;
    cout << "m=" << m << ", segments=" << segments << " fs=" << fs << endl;

    spectrum.get_freq(f);
    fp = fopen("modelspectrum_percentiles.out","w");
    j  = 1;
    f0 = log(f[1]);
    f1 = log(f[n2-1]);
    for (i=0;i<N;i++) { 
      s  = static_cast<double>(i)/static_cast<double>(N);
      f_int = exp((1.0-s)*f0 + s*f1);
      while (f[j]<f_int and j<n2-1) j++;
      s = (f_int - f[j-1])/(f[j]-f[j-1]);
      G_5_int  = G_5[j-1]*(1-s)  + G_5[j]*(s);
      G_50_int = G_50[j-1]*(1-s) + G_50[j]*(s);
      G_95_int = G_95[j-1]*(1-s) + G_95[j]*(s);
      fprintf(fp,"%e %e %e %e\n",f_int,G_5_int,G_50_int,G_95_int);
    } 
    fclose(fp);

    //--- Clean up mess
    delete[] y;
    delete[] f;
    delete[] G_5;
    delete[] G_50;
    delete[] G_95;
    delete[] dummy;
    delete[] G;
  }


//==============================================================================
// Main program
//==============================================================================

  int main(int argc, char* argv[]) 
  {
    time_t    start,end;

    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      Control &control = Control::getInstance("modelspectrum.ctl");
    } else if (argc==2) {
      Control &control = Control::getInstance(argv[1]);
    } else {
      cerr << "correct usage: modelspectrum [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    cout << endl
         << "************************************" << endl
         << "    modelspectrum, version " << VERSION <<endl
         << "************************************" << endl;

    //--- start modelspectrum
    ModelSpectrum PSD;

    time(&start);
    PSD.compute();
    time(&end);

    cout << "ModelSpectrum took " << difftime(end,start) << " s" << endl;

    return EXIT_SUCCESS;
  }

