/*! \filename FindOffset.cpp
 *  \author   Machiel Bos
 *
 * Perform a single cycle of offset detection:
 *  1) Compute noise and SLT model parameters
 *  2) Compute for each epoch the -2*log(L) - 2*log(f_s) - 2*log(f_theta)
 *  3) Save results to file and show lowest value on screen
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
  #include "Control.h"
  #include "Observations.h"
  #include "DesignMatrix.h"
  #include "Likelihood.h"
  #include "Minimizer.h"
  #include "JSON.h"
  #include <iostream>
  #include <iomanip>
  #include <ostream>
  #include <fstream>
  #include <ctime>
  #include <cstdlib>

//==============================================================================
// Main program
//==============================================================================


  int main(int argc, char *argv[])
  {
    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      Control &control = Control::getInstance("findoffset.ctl");
    } else if (argc==2) {
      Control &control = Control::getInstance(argv[1]);
    } else {
      cerr << "correct usage: findoffset [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Start estimatetrend
    cout << endl
         << "************************************" << endl
         << "    findoffset, version " << VERSION << endl
         << "************************************" << endl;

    //--- Here the rest of the variable declarations occur
    fstream        fp;
    const double   TINY=1.0e-6;
    JSON           &json = JSON::getInstance("findoffset.json");
    Minimizer      minimizer;
    Likelihood     &likelihood=Likelihood::getInstance();
    DesignMatrix   *designmatrix = DesignMatrix::getInstance();
    Observations   &observations = Observations::getInstance();
    time_t         start,end;
    int            i,Nparam,m,n=-1,i_min;
    double         dt,*t,*x,*param=NULL,t_min,BIC_c_min,*BIC_c=NULL;
    string         true_param[10];
    Control        &control = Control::getInstance();

    //--- Are we using a-priori noise parameter values?
    try {
      control.get_name_list("TrueParams",true_param,n);
      param = new double[n];
      for (i=0;i<n;i++) {
        param[i] = atof(true_param[i].c_str());
        cout << "param " << i << ") =" << param[i] << endl;
      }
      //--- We also need to know BIC_c without any extra offset being estimated
      //    Assume user is smart enough to provide valid parameter values
      //    and directly call compute_LeastSquares and not compute.
      likelihood.compute_LeastSquares(param);
      likelihood.show_leastsquares();
      likelihood.compute_L_and_ICs(param);
      cout << "BIC_c    =" << likelihood.get_BIC_c() << endl;
    } catch (exception &e) {
    
      time(&start);
      //--- Use MLE to estimate noise parameters. This will also show on the
      //    screen the BIC_c for the case no extra offset is being estimated,
      //    besides trend value.   
      minimizer.solve();

      //--- Remember noise parameter values
      Nparam = minimizer.get_Nparam();
      param = new double[Nparam];
      minimizer.get_param(param);
      time(&end);
      cout << "Total computing time: " << difftime(end,start) << " sec" << endl;
    }

    //--- Add another column for the new offset. This implies creating a new
    //    design matrix H and Likelihood instance so: reset them. 
    observations.get_values(m,&t,&x);
    dt = observations.get_fs()*24.0*3600.0;
    observations.add_offset(t[0]+dt);
    designmatrix->resetInstance();
    designmatrix=DesignMatrix::getInstance();
    likelihood.reset_method();

    //--- Allocate memory for BIC_c array
    BIC_c = new double[m];
    memset(BIC_c,0.0,m*sizeof(double));

    //--- Likelihood class has been reset. Make sure compute_LeastSquares is
    //    called to fill A1, A2, y1, y2, etc...
    //    To call 'compute' instead of 'compute_LeastSquares' is a little safer
    //    in case the parameter values are out of acceptable range.
    likelihood.compute(param);
   
    //--- Now compute BIC_c for all possible offset locations 
    fp.open("findoffset.out",ios::out);
    time(&start);
    likelihood.compute_BIC_cs(BIC_c);
  
    i_min     = 1; 
    BIC_c_min = BIC_c[1]; 
    t_min     = t[1];
    if (BIC_c[1]<9.0e99)
      fp << fixed << setprecision(3) << t[1] << "  "
       				   << setprecision(6) << BIC_c[1] << endl;
    for (i=2;i<m;i++) {
      if (BIC_c[i]<BIC_c_min) {
        i_min     = i;
        t_min     = t[i];
        BIC_c_min = BIC_c[i];
      }
      if (BIC_c[i]<9.0e99)
        fp << fixed << setprecision(3) << t[i] << "  "
         			     << setprecision(6) << BIC_c[i] << endl;
    }
    time(&end);
    fp.close();
    cout << "FindOffset i_min =" << i_min << endl;
    cout << "FindOffset MJD   =" << t_min << endl;
    cout << "FindOffset BIC_c =" << BIC_c_min << endl;
    cout << "Total computing time: " << difftime(end,start) << " sec" << endl;
   
    //--- Free memory
    if (BIC_c!=NULL)  delete[] BIC_c;
    if (param!=NULL)  delete[] param;
 
    return EXIT_SUCCESS;
  }

