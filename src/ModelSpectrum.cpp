/*! \file    ModelSpectrum.cpp
 *  \author  Machiel Bos
 *  \version 1.0
 *
 * For my sea-level trend work I need to make plots of the spectral density
 * for the estimated noise model versus observed spectral density. This
 * program provides a subroutine for making the modelled spectrum.
 *
 * \date  22/7/2011  CIIMAR,  Porto
 * \date   5/6/2012  Santa Clara
 */
//==============================================================================
  #include <cmath>
  #include <cstdlib>
  #include <cstdio>
  #include <iostream>
  #include <ostream>
  #include "ModelSpectrum.h"

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================

//---!!--------------------------------------------------
  ModelSpectrum::ModelSpectrum(void) : pi(4.0*atan(1.0)),
                                       tpi(2.0*pi)
//---!!--------------------------------------------------
  {
    using namespace std;
    int               i,j;
    double            T,*theta;
    complex<double>   a,b,c,w;

    //--- Need to know the standard deviation of epsilon (the white noise)
    cout << "Enter the standard deviation of the driving noise: ";
    cin >> sigma;

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
    NoiseModel  &noisemodel=NoiseModel::getInstance();
    fstream     fp;
    string      filename;
    Control     &control=Control::getInstance();

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

    //--- Clean up mess
    delete[] G;
    delete[] f;
    if (params_fixed!=NULL) delete[] params_fixed;
  }
    

//==============================================================================
// Main program
//==============================================================================

  int main(int argc, char* argv[]) 
  {
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
         << "    modelspectrum, version " << VERSION << endl
         << "************************************" << endl << endl;

    //--- start modelspectrum
    ModelSpectrum PSD;

    PSD.compute();

    return EXIT_SUCCESS;
  }

