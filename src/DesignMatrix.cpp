/*! \file   DesignMatrix.cpp
 *  \author Machiel Bos
 *
 * This class provides the Design Matrix H that is needed in the 
 * Least-Squares part of the computation of the Likelihood.
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
  #include "DesignMatrix.h"
  #include "Calendar.h"
  #include <iostream>
  #include <ostream>
  #include <fstream>
  #include <iomanip>
  #include <cmath>
  #include <cstdlib>
  #include <vector>

  #define TINY 1.0e-5
//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================

//++++++ Singleton stuff +++++++++++++++++++++
DesignMatrix* DesignMatrix::singleton = NULL;
//++++++++++++++++++++++++++++++++++++++++++++


//---!!-------------------------------------------------
  DesignMatrix::DesignMatrix(void) : tpi(8.0*atan(1.0)),
			        pi(0.5*tpi), median(0.0)
//---!!-------------------------------------------------
  {
    using namespace std;
    fstream                   fp;
    ofstream                  myfile;
    bool                      estimate_offsets,estimate_postseismic;
    bool                      estimate_sse;
    int                       i,j,k,l,year,month,day,hour,minute,Nnumbers=0;
    size_t                    i_multi;
    double                    *x,second,MJD,T,J,phi,s,T0,T1,coeff;
    string                    periodic_signals[40],numbers[40];
    Observations              &observations=Observations::getInstance();
    Control                   &control=Control::getInstance();
    Calendar                  calendar;
    string                    str,line,filename;
    stringstream              fs_line;
    string                    whitespaces (" \t\f\v\n\r");
    size_t                    found;
    vector<double>            row;
    vector<vector <double> >  multivariate;

    //--- Determine some general properties
    observations.get_values(m,&t,&x);
    Ngaps = observations.number_of_gaps();
    median = observations.get_median();
    dt = 1.0/(24.0*3600.0*observations.get_fs()); // unit is MJD days!!!
    t0 = t[0];
#ifdef DEBUG
    cout << "m=" << m << endl;
    cout << "dt=" << dt << endl;
#endif

    try {
      control.get_string("PhysicalUnit",unit);
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    //--- Are we estimating one trend or several?    
    try {
      estimate_multitrend = control.get_bool("estimatemultitrend");
    }
    catch (exception &e) {
      estimate_multitrend = false;
    }

    //--- How many segments or what's the degree of the polynomial
    if (estimate_multitrend==true) {
      observations.get_breaks(breaks);
      n_breaks = breaks.size();
      n        = 2+n_breaks;
    } else {
      //--- Janusz and Anna want to investigate higher degrees
      try {
        degree_polynomial = control.get_int("DegreePolynomial");
      }
      catch (exception &e) {
        cout << "No Polynomial degree set, using offset + linear trend" << endl;
        degree_polynomial = 1;
      }
      if (degree_polynomial<0 || degree_polynomial>6) {
        cerr << "Only polynomial degrees  between 0 and 6 are allowed" << endl;
        exit(EXIT_FAILURE);
      }
      n = degree_polynomial+1;
    }

    //--- Do we need to include seasonal signal
    try {
      varying_seasonal = control.get_bool("estimatevaryingseasonal");
    }
    catch (exception &e) {
      varying_seasonal = false;
    }
    if (varying_seasonal==false) {
      try {
        seasonal_signal      = control.get_bool("seasonalsignal");
        halfseasonal_signal  = control.get_bool("halfseasonalsignal");
      } 
      catch (exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
      } 
    } else {
      seasonal_signal = halfseasonal_signal = false;
    }

    //--- Varying seasonal signals?
    if (varying_seasonal==true) {
      try {
        varyingseasonal_N = control.get_int("varyingseasonal_N");
      } 
      catch (exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
      } 
      cout << "varyingseasonal_N = " << varyingseasonal_N << endl;
      index_seasonal = n;
      n += 4*(1+varyingseasonal_N);
    } else {
      if (seasonal_signal==true)     n +=2;
      if (halfseasonal_signal==true) n +=2;
    }

    //--- An additional number of 40 periodic signals may be added
    //    if keyword is not found, then that is okay too.
    try {
      control.get_name_list("periodicsignals",periodic_signals,
							n_periodic_signals);
    }
    catch (exception &e) {
      cout << "No extra periodic signal is included." << endl;
      n_periodic_signals=0;
    }
    if (n_periodic_signals>40) {
      cerr << "Only up to 40 periodic signals are allowed!" << endl;
      exit(EXIT_FAILURE);
    }
    n += 2*n_periodic_signals;
    if (n_periodic_signals>0) {
      periods = new double[n_periodic_signals];
      for (i=0;i<n_periodic_signals;i++) {
        periods[i] = atof(periodic_signals[i].c_str());
        if (periods[i]==0.0) {
          cerr << "Cannot have a 0 day period or Hector could not interpret"
               << ": " << periodic_signals[i] << endl;
          exit(EXIT_FAILURE);
        }
      }
    } else {
      periods = NULL;
    }

    //--- Maybe we need to include offsets too
    try {
      estimate_offsets = control.get_bool("estimateoffsets");
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
    if (estimate_offsets==true) {
      //--- For the prior it is convenient to know where the estimated
      //    sizes are stored in vector param.
      index_offset = n;

      observations.get_offsets(offsets);
      n_offsets  = offsets.size();
      n         += n_offsets;
    } else {
      n_offsets  = 0;
    }

    //--- Check also for postseicmic relaxation modelling
    try {
      estimate_postseismic = control.get_bool("estimatepostseismic");
    }
    catch (exception &e) {
      estimate_postseismic=false;
    }
    if (estimate_postseismic==true) {
      observations.get_postseismiclog(postseismiclog);
      n_postseismiclog  = postseismiclog.size();
      n                += n_postseismiclog;
      observations.get_postseismicexp(postseismicexp);
      n_postseismicexp  = postseismicexp.size();
      n                += n_postseismicexp;
    } else {
      n_postseismiclog = 0;
      n_postseismicexp = 0;
    }

    //--- Check for slow slip event modelling
    try {
      estimate_sse = control.get_bool("estimateslowslipevent");
    }
    catch (exception &e) {
      estimate_sse=false;
    }
    if (estimate_sse==true) {
      observations.get_ssetanh(ssetanh);
      n_ssetanh  = ssetanh.size();
      n         += n_ssetanh;
    } else {
      n_ssetanh  = 0;
    }

    //--- Do we need to include a geophysical signals in matrix H?
    try {
      estimate_multivariate = control.get_bool("estimatemultivariate");
    }
    catch (exception &e) {
      estimate_multivariate=false;
    }
    if (estimate_multivariate==false) {
      n_channels=0; //--- Remember that we have no extra channels!
    } else {
      try {
        control.get_string("MultiVariateFile",filename);
      } catch (exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
      }
      //--- Open multivariate file
      fp.open(filename.c_str(),ios::in);
      if (!fp.is_open()) {
        cerr << "Could not open " << filename << endl;
        exit(EXIT_FAILURE);
      }
      //--- Read file
      while (getline(fp,line)) {  // read the whole line into string 'line'
        row.clear();
        fs_line.clear();
        //--- Remove trailing spaces since >> chokes on them
        found = line.find_last_not_of(whitespaces);
        if (found!=std::string::npos) {
          line.erase(found+1);
        }
        fs_line << line;   // copy line into the string-stream

        while (fs_line.good()==true && fs_line.eof()!=true) {
          fs_line >> str;
          if (str.length()>0) {
            row.push_back(atof(str.c_str()));
          }
        }
        //--- Save row
        multivariate.push_back(row);
      }
      fp.close();

      //--- Sanity check
      if (static_cast<size_t>(m)>multivariate.size()) {
        cerr << "There are " << m << " observations but only " 
             << multivariate.size()
             << " rows in the multivariate file!" << endl;
        exit(EXIT_FAILURE);
      }
      n_channels = row.size()-1;
      n += n_channels;
    }

    //**** End of getting information about what to put in Design Matrix ****


    //--- Allocate memory
    try {
      H = new double[m*n];
      memset(H,0,(m*n)*sizeof(double));
      if (Ngaps>0)  F = new double[Ngaps*m];
      else          F = NULL;
    }
    catch(const bad_alloc &) {
      cerr << "Need more memory to store Design Matrix" << endl;
      cerr << "m=" << m << ", n=" << n << ", Ngaps=" << Ngaps << endl;
      exit(EXIT_FAILURE);
    }

    //--- See if there is a reference Epoch we should use for th
    try {
      control.get_name_list("ReferenceEpoch",numbers,Nnumbers);
    } 
    catch (exception &e) {
#ifdef DEBUG
      cerr << e.what() << endl;
#endif
    }
    if (Nnumbers==0) { // No ReferenceEpoch found, use middle of time-series
      th = 0.5*(t[0] + t[m-1]);
    } else if (Nnumbers==3) {
      year  = atoi(numbers[0].c_str());
      month = atoi(numbers[1].c_str());
      day   = atoi(numbers[2].c_str());
      hour  = 0;
      minute= 0;
      second= 0.0;
      th    = calendar.compute_MJD(year,month,day,hour,minute,second);
    } else {
      cerr << "Correct usage: ReferenceEpoch year month day" << endl;
      cout << "Nnumbers=" << Nnumbers << endl;
      exit(EXIT_FAILURE);
    }

    //--- Construct the design matrix H
    i_multi = 0;
    for (i=0;i<m;i++) {
      if (estimate_multitrend==true) {
        //--- multi-trend 
        H[i] = 1.0;  //--- nominal bias
        j = 0;
        while (j<n_breaks && t[i]+TINY>breaks[j]) {
          H[i+(1+j)*m] = 0.0;
          j++;
        }
        if (j==0) {
          H[i+m] = t[i]-t0;
        } else {
          H[i+(1+j)*m] = t[i]-breaks[j-1];
        }
        j=2+n_breaks;
          
      } else {
        //--- polynomial trend
        for (j=0,J=0.0;j<=degree_polynomial;j++,J+=1.0) {
          H[i+j*m] = pow(t[i]-th,J);
        }
        j = degree_polynomial+1;
      }

      //--- Seasonal signal
      if (varying_seasonal==true) {
        phi = tpi*(t[i]-51544.0)/365.25;
        
        //--- Define basis function
        s  = 2.0*(t[i]-t0)/(t[m-1]-t[0]);
        T0 = 1.0;
        T1 = s;
        for (k=0;k<=varyingseasonal_N;k++) {
          if (k==0) {
            T = T0;
          } else if (k==1) {
            T = T1;
          } else {
            T  = 2.0*s*T1 - T0;
            T0 = T1;
            T1 = T;
          }
          H[i+(j+k*4+0)*m] = T*cos(phi);
          H[i+(j+k*4+1)*m] = T*sin(phi);
          H[i+(j+k*4+2)*m] = T*cos(2.0*phi);
          H[i+(j+k*4+3)*m] = T*sin(2.0*phi);
        }
        j += (1+varyingseasonal_N)*4;
      } else {
        //--- Constant seasonal signal
        if (seasonal_signal==true) {
          H[i+j*m] = cos(tpi*(t[i]-51544.0)/365.25);            j++;
          H[i+j*m] = sin(tpi*(t[i]-51544.0)/365.25);            j++;
        }
        if (halfseasonal_signal==true) {
          H[i+j*m] = cos(2.0*tpi*(t[i]-51544.0)/365.25);        j++;
          H[i+j*m] = sin(2.0*tpi*(t[i]-51544.0)/365.25);        j++;
        }
      }

      //--- Periodic signals
      for (k=0;k<n_periodic_signals;k++) {
        H[i+j*m] = cos(tpi*(t[i]-51544.0)/periods[k]);        j++;
        H[i+j*m] = sin(tpi*(t[i]-51544.0)/periods[k]);        j++;
      }

      //--- Offsets
      for (l=0;l<n_offsets;l++) {
        if (t[i]+TINY>offsets[l]) H[i + (j+l)*m] = 1.0;
        else                      H[i + (j+l)*m] = 0.0;
      }
      j += n_offsets;

      //--- Post-seismic stuff
      for (l=0;l<n_postseismiclog;l++) {
        MJD = postseismiclog[l].MJD;
        T   = postseismiclog[l].T;
        if (t[i]+TINY>MJD) H[i + (j+l)*m] = log1p((t[i]-MJD)/T);
      }
      j += n_postseismiclog;
      for (l=0;l<n_postseismicexp;l++) {
        MJD = postseismicexp[l].MJD;
        T   = postseismicexp[l].T;
        if (t[i]+TINY>MJD) H[i + (j+l)*m] = 1.0 - exp(-(t[i]-MJD)/T);
      }
      j += n_postseismicexp;

      //--- Slow Slip event
      for (l=0;l<n_ssetanh;l++) {
        MJD = ssetanh[l].MJD;
        T   = ssetanh[l].T;
        H[i + (j+l)*m] = 0.5 * ( tanh((t[i]-MJD)/T) - 1.0);
      }
      j += n_ssetanh;

      if (estimate_multivariate==true) {
        row = multivariate[i_multi];
        //--- Skip rows until first observation, afterwards this is not needed
        while (t[i]-row[0]>0.001 && i_multi<multivariate.size()-1) {
          i_multi++;
          row = multivariate[i_multi];
        }
        //--- Have we somehow already reached the end?
        if (i_multi==multivariate.size()-1) {
          cerr << "Multivariate time series is too short!" << endl;
          exit(EXIT_FAILURE);

        //--- Problem with time step: not exact same epoch in both series
        } else if (fabs(row[0]-t[i])>TINY) {
          cerr << "Epoch of multivariate file is: " << row[0] << endl;
          cerr << "Time in time series is: " << t[i] << endl;
          cerr << "This should not happen!" << endl;
          exit(EXIT_FAILURE);
        }

        //--- row is now for t[i] -> copy columns 
        for (l=0;l<n_channels;l++) H[i+(j+l)*m] = row[l+1];

        //--- Already move to next row
        i_multi++;

        //--- Move n_channel columns to the right
        j += n_channels; 
      }
    }

#ifdef DEBUG
    myfile.open("H.dat", ios::out);
    for (i=0;i<m;i++) {
      if (isnan(x[i])==false) {
        for (j=0;j<n;j++) myfile << H[i + j*m] << "  ";
        myfile << endl;
      }
    }
    myfile.close();
#endif

    //--- Construct selection matrix F
    if (Ngaps>0) {
      //--- Make zero
      //cblas_dscal(Ngaps*m,0.0,F,1); !!! Does not work NaN*0.0=NAN !!!!
      memset(F,0,Ngaps*m*sizeof(double));

      //--- Put the ones in the right places
      j=0;
      for (i=0;i<m;i++) {
        if (std::isnan(x[i])) {
          F[i + j*m] = 1.0;
          j++;
        }
      }
    } else {
      F = NULL;
    }

    //--- For numerical stability, it is advantageous to orthogonalise
    //    the basis functions for the seasonal signal using Gram-Schmidt
    //    algorithm
    if (varying_seasonal==true) {
      for (i=0;i<(1+varyingseasonal_N)*4;i++) {
        for (j=0;j<i;j++) {
          k = m*(index_seasonal + i);
          l = m*(index_seasonal + j);
          coeff = -cblas_ddot(m,&H[k],1,&H[l],1)/cblas_ddot(m,&H[l],1,&H[l],1);
          cblas_daxpy(m,coeff,&H[l],1,&H[k],1); 
        }
      }  
    }

    //--- Free some memory
    delete[] x;
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
    if (t!=NULL)        delete[] t;
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



/*! Compute correct total amplitude and uncertainty using Rice distribution
 */
//------------------------------------------------------------------------
  void DesignMatrix::compute_Amp(double Ac, double As, 
			  double sigma_in, double& Amp, double& sigma_out)
//------------------------------------------------------------------------
  {
    double   nu,L12;

    nu        = sqrt(Ac*Ac + As*As);
    L12       = gsl_sf_hyperg_1F1(-0.5,1.0,-pow(nu/sigma_in,2.0)/2.0);

    Amp       = sqrt(pi/2.0)*sigma_in*L12;
    sigma_out = sqrt(2.0*pow(sigma_in,2.0) + pow(nu,2.0) - 
				pi*pow(sigma_in,2.0)/2.0*pow(L12,2.0));
  }



/*! Compute phase value and its uncertainty
 */
//------------------------------------------------------------------------
  void DesignMatrix::compute_Pha(double Ac, double As, double sigma_Ac,
                          double sigma_As, double& Pha, double& sigma_out)
//------------------------------------------------------------------------
  {
    Pha       = atan2(As, Ac);
    sigma_out = sqrt(pow(Ac/(pow(As,2.0)+pow(Ac,2.0))*sigma_As, 2.0) +
                    pow(-As/(pow(As,2.0)+pow(Ac,2.0))*sigma_Ac, 2.0));
  }


namespace {
//------------------------------------------------------------------------
  std::string values_formatted(const std::string& name, double value,
                                    double error, const std::string &unit)
//------------------------------------------------------------------------
  {
    using namespace std;
    stringstream ss;
    ss << setw(36) << name << " : " << fixed
       << setw(13) << setprecision(3) << value << " +/- "
       << setw(7)  << setprecision(3) << error << " "
       << unit;
    return ss.str();
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
    using namespace std;
    const double   rad=tpi/360.0;
    int            i,j,year,month,day,hour,minute;
    size_t         k;
    double         ds = 365.25,second,corr_pos = 0.0,corr_vel = 0.0,Amp,Pha,sigma_out;
    double         *xhat,*theta_dummy;
    fstream        fp;
    Calendar       calendar;
    ostringstream  buff;

    i=0;

    if (estimate_multitrend==true) {
      cout << values_formatted("bias", theta[0]+median, error[0], unit) << endl;
      for (j=0;j<=n_breaks;j++) {
        cout << values_formatted("trend", theta[1+j]*ds, error[1+j]*ds, unit)
             << "/year" << endl;
      }
      i += 2+n_breaks;

    } else {
      calendar.compute_date(th,year,month,day,hour,minute,second);
      cout << values_formatted("bias", theta[0]+median, error[0], unit)
           << " (at ";
      cout.fill('0');
      cout << setw(4) << year << "/"
           << setw(2) << month << "/"
           << setw(2) << day << ", "
           << setw(2) << hour << ":"
           << setw(2) << minute << ":"
           << setw(6) << setprecision(3) << second << ")" << endl;
      cout.fill(' ');

      if (degree_polynomial>0) {
        cout << values_formatted("trend", theta[1]*ds, error[1]*ds, unit)
             << "/year" << endl;
      }
      if (degree_polynomial>1) {
        cout << values_formatted("quadratic (half acceleration)",
                                theta[2]*ds*ds, error[2]*ds*ds, unit)
             << "/year^2" << endl;
      }
      for (j=3;j<=degree_polynomial;j++) {
        buff << "degree " << j;
        cout << values_formatted(buff.str(), theta[j]*pow(ds,1.0*j),
                                      error[j]*pow(ds,1.0*j), unit)
             << "/year^" << j << endl;
        buff.str("");
      }
      i += degree_polynomial + 1;
    }

    //--- Seasonal signal
    if (varying_seasonal==true) {
      theta_dummy = new double[n];
      xhat = new double[m];
      for (j=0;j<n;j++) theta_dummy[j]=0.0;
      for (j=0;j<(1+varyingseasonal_N)*4;j++) {
	theta_dummy[index_seasonal+j]=theta[index_seasonal+j];
        //cout << theta_dummy[index_seasonal+i] << endl;
      }
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                     m,1,n,1.0,H,m,theta_dummy,n,0.0,xhat,m); 
      cout.precision(5);
      fp.open("varying_seasonal.mom",ios::out);
      if (!fp.is_open()) {
        cerr << "Could not open varying_seasonal.mom" << endl;
        exit(EXIT_FAILURE);
      }
      //--- Save the data to the output file
      for (j=0;j<m;j++) {
        fp << fixed << t[j] << "  " << xhat[j] << endl;
      }
      fp.close();

      //--- free memory
      delete[] theta_dummy;
      delete[] xhat;

      //--- Move to next parameter
      i += (1+varyingseasonal_N)*4;

    } else {
      if (seasonal_signal==true) {
        cout << values_formatted("cos yearly", theta[i], error[i], unit) << endl;
        i++;
        cout << values_formatted("sin yearly", theta[i], error[i], unit) << endl;
        i++;
        compute_Amp(theta[i-2],theta[i-1],0.5*(error[i-2]+error[i-1]),
							   Amp,sigma_out);
        cout << values_formatted("Amp yearly", Amp, sigma_out, unit) << endl;
        compute_Pha(theta[i-2],theta[i-1],error[i-2],error[i-1],Pha,sigma_out);
        cout << values_formatted("Pha yearly", Pha/rad, sigma_out/rad, "degrees") << endl;
      }

      if (halfseasonal_signal==true) {
        cout << values_formatted("cos hyearly", theta[i], error[i], unit) << endl;
        i++;
        cout << values_formatted("sin hyearly", theta[i], error[i], unit) << endl;
        i++;
        compute_Amp(theta[i-2],theta[i-1],0.5*(error[i-2]+error[i-1]),
							   Amp,sigma_out);
        cout << values_formatted("Amp hyearly", Amp, sigma_out, unit) << endl;
        compute_Pha(theta[i-2],theta[i-1],error[i-2],error[i-1],Pha,sigma_out);
        cout << values_formatted("Pha hyearly", Pha/rad, sigma_out/rad, "degrees") << endl;
      }
    }

    //--- Periodic signals
    for (j=0;j<n_periodic_signals;j++) {
      buff << "cos (T=" << fixed << setw(8) << setprecision(2) << periods[j] << ")";
      cout << values_formatted(buff.str(), theta[i], error[i], unit) << endl;
      buff.str("");
      i++;
      buff << "sin (T=" << fixed << setw(8) << setprecision(2) << periods[j] << ")";
      cout << values_formatted(buff.str(), theta[i], error[i], unit) << endl;
      buff.str("");
      i++;

      compute_Amp(theta[i-2],theta[i-1],0.5*(error[i-2]+error[i-1]),
							   Amp,sigma_out);
      buff << "amp (T=" << fixed << setw(8) << setprecision(2) << periods[j] << ")";
      cout << values_formatted(buff.str(), Amp, sigma_out, unit) << endl;
      buff.str("");
      compute_Pha(theta[i-2],theta[i-1],error[i-2],error[i-1],Pha,sigma_out);
      buff << "pha (T=" << fixed << setw(8) << setprecision(2) << periods[j] << ")";
      cout << values_formatted(buff.str(), Pha/rad, sigma_out/rad, "degrees") << endl;
      buff.str("");
    }
    for (j=0;j<n_offsets;j++) {
      // correction values are set to non-zero at break points if
      // there is a multitrend estimation
      corr_pos = 0.0;
      corr_vel = 0.0;

      if (estimate_multitrend && breaks.size()>0) {
        for (k=0;k<breaks.size();k++) {
          // calculate offset correction value at break point
          if (fabs(offsets[j]-breaks[k])<TINY) {
            if (k==0) {
              corr_pos = theta[1]*(breaks[k]-t0);
              corr_vel = sqrt(pow(error[i],2.0) + 
			   pow(error[1]*(breaks[k]-t0),2.0)) - error[i];
            } else {
              corr_pos = theta[1+k]*(breaks[k]-breaks[k-1]);
              corr_vel = sqrt(pow(error[i],2.0) + pow(error[1+k]*(breaks[k] 
						-breaks[k-1]),2.0)) - error[i];
            }
          }
        }
      }
      buff << "offset at " << fixed << setw(10) << setprecision(4) << offsets[j];
      cout << values_formatted(buff.str(), theta[i]-corr_pos, error[i]+corr_vel, unit)
           << endl;
      buff.str("");
      i++;
    }
    for (j=0;j<n_postseismiclog;j++) {
      buff << "log relax at " << setw(10) << setprecision(4) << postseismiclog[j].MJD
           << " (T=" << fixed << setw(8) << setprecision(2) << postseismiclog[j].T << ")";
      cout << values_formatted(buff.str(), theta[i], error[i], unit) << endl;
      buff.str("");
      i++;
    }
    for (j=0;j<n_postseismicexp;j++) {
      buff << "exp relax at " << setw(10) << setprecision(4) << postseismicexp[j].MJD
           << " (T=" << fixed << setw(8) << setprecision(2) << postseismicexp[j].T << ")";
      cout << values_formatted(buff.str(), theta[i], error[i], unit) << endl;
      buff.str("");
      i++;
    }
    for (j=0;j<n_ssetanh;j++) {
      printf("tanh sse at %10.4lf (T=%8.2lf) : %7.2lf +/- %5.2lf %s\n",
		ssetanh[j].MJD, ssetanh[j].T,
					theta[i],error[i],unit.c_str()); i++;
    }
    for (j=0;j<n_channels;j++) {
      buff << "scale factor of channel " << fixed << setw(1) << j+1;
      cout << values_formatted(buff.str(), theta[i], error[i], unit) << endl;
      buff.str("");
      i++;
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
    Observations  &observations=Observations::getInstance();
   
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
					     m,1,n,1.0,H,m,theta,n,0.0,xhat,m);
    observations.set_xhat(xhat);
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
