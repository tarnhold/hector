/*! \file   SimulateNoise.cpp
 *  \author Machiel Bos
 *
 * Implementation of the Hosking (1984) method.
 *
 *  This script is part of Hector 1.7.2
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
    using namespace std;
    Control          &control=Control::getInstance();
    NoiseModel       &noisemodel=NoiseModel::getInstance();
    int              i,j,k,n_simulations,m,ms;
    string           directory,label,filename;
    double           *y,*MJD,dt,ts;
    FILE             *fp;
    stringstream     ss;
  
    //--- Which noise models must be used
    try {
      control.get_string("SimulationDir",directory);
      if (directory.at(directory.length()-1)!='/') {
        directory += "/";
      }
      control.get_string("SimulationLabel",label);
      n_simulations = control.get_int("NumberOfSimulations");
      m             = control.get_int("NumberOfPoints");
      dt            = control.get_double("SamplingPeriod");
      ms            = control.get_int("TimeNoiseStart"); 
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    //--- Prepare NoiseModel impulse responses
    noisemodel.setup_MonteCarlo(ms+m);

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
      ss.str("");
      ss << i;
      filename = directory + label + "_" + ss.str() + ".mom";
      fp = fopen(filename.c_str(),"w");
      if (fp==NULL) {
        cerr << "Could not open" << filename << endl;
        exit(EXIT_FAILURE);
      }

      //--- Create the synthetic noise
      noisemodel.create_noise(ms+m,y);
      cout << "noise created" << endl;

      //--- write results to file
      fprintf(fp,"# sampling period %lf\n",dt);
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
    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      Control &control = Control::getInstance("simulatenoise.ctl");
    } else if (argc==2) {
      Control &control = Control::getInstance(argv[1]);
    } else {
      cerr << "correct usage: simulatenoise [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    cout << endl
         << "************************************" << endl
         << "    simulatenoise, version " << VERSION << endl 
         << "." << SUBVERSION
         << "************************************" << endl;

    //--- Now it's save to initiate the SimulateNoise class
    SimulateNoise    simulate;

    return EXIT_SUCCESS;
  }
