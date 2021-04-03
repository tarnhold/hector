/*! \file   SimulateNoise.cpp
 *  \author Machiel Bos
 *
 * Implementation of the Hosking (1984) method.
 *
 * \date 15/1/2013
 */
//==============================================================================
  #include "SimulateNoise.h"
  #include <iostream>
  #include <ostream>
  #include <cstring>

//==============================================================================
// Subroutines
//==============================================================================

/*! Class constructor
 */
//---!!-----------------------------
  SimulateNoise::SimulateNoise(void)
//---!!-----------------------------
  {
    Control       *control=Control::getInstance();
    NoiseModel    *noisemodel=NoiseModel::getInstance();
    int           i,j,k,n_simulations,m,ms;
    char          directory[80],label[30],filename[120];
    double        *y,*MJD,dt,ts;
    FILE          *fp;

    using namespace std;
    //--- Which noise models must be used
    try {
      control->get_string("SimulationDir",directory);
      control->get_string("SimulationLabel",label);
      n_simulations = control->get_int("NumberOfSimulations");
      m             = control->get_int("NumberOfPoints");
      dt            = control->get_double("SamplingPeriod");
      ms            = control->get_int("TimeNoiseStart"); 
    }
    catch (const char* str) {
      cerr << str << endl;
      exit(EXIT_FAILURE);
    }

    //--- Prepare NoiseModel impulse responses
    noisemodel->setup_MonteCarlo(ms+m);

    //--- Store results in array y
    MJD = new double[m];
    y   = new double[ms+m];

    //--- Already create time array
    MJD[0] = 51544.0; // 1 January 2000
    for (i=1;i<m;i++) {
      MJD[i] = MJD[i-1] + dt;
    }

    //--- Run all simulations
    for (i=0;i<n_simulations;i++) {
      cout << "simulation run:" << i << endl;
      //--- Open file to store time-series
      k = strlen(directory);
      if (directory[k-1]=='/') {
        sprintf(filename,"%s%s_%d.mom",directory,label,i);
      } else {
        sprintf(filename,"%s/%s_%d.mom",directory,label,i);
      }
      fp = fopen(filename,"w");
      if (fp==NULL) {
        cerr << "Could not open" << filename << endl;
        exit(EXIT_FAILURE);
      }

      //--- Create the synthetic noise
      noisemodel->create_noise(ms+m,y);
      cout << "noise created" << endl;

      //--- write results to file
      for (j=ms;j<(ms+m);j++) {
        fprintf(fp,"%lf  %lf\n",MJD[j-ms],y[j]);
      }

      //--- close file
      fclose(fp);
    }

    //--- Free memory
    delete[] MJD;
    delete[] y;
  }



/*!
 */
//---!!------------------------------
  SimulateNoise::~SimulateNoise(void)
//---!!------------------------------
  {
  }



//==============================================================================
// Main program
//==============================================================================

  int main(int argc, char* argv[])
  {
    Control         *control;

    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      control = Control::getInstance("simulatenoise.ctl");
    } else if (argc==2) {
      control = Control::getInstance(argv[1]);
    } else {
      cerr << "correct usage: simulatenoise [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    cout << endl
         << "************************************" << endl
         << "    simulatenoise, version " << VERSION << endl
         << "************************************" << endl;

    //--- Now it's save to initiate the SimulateNoise class
    SimulateNoise    simulate;

    return EXIT_SUCCESS;
  }
