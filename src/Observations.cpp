/*! \file   Observations.cpp
 *  \author Machiel Bos
 *
 * I will work with various data formats. This forces me to treat the 
 * reading of observations a bit more generally than I normally do.
 *
 * Note that when reading the data, I first omit storing the gaps. This
 * helps later finding them back because of t_(i+1)!=t_i+dt. Otherwise
 * I needed to count the NaN's first to determine the width of the gap
 * before I could compute the slope. Nevertheless, once the data is
 * stored in Memory, I interpolate the data or fill them with NaN's.
 *
 * \date 27/7/2011  CIIMAR, Porto
 * \date 16/1/2012  Santa Clara
 */
//==============================================================================
  #include "Observations.h"
  #include <iostream>
  #include <ostream>
  #include <iomanip>
  #include <fstream>
  #include <string>
  #include <cstdlib>
  #include <cstdio>
  #include <cmath>

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================

//++++++ Singleton stuff +++++++++++++++++++++
bool Observations::instanceFlag = false;
Observations* Observations::singleton = NULL;
//++++++++++++++++++++++++++++++++++++++++++++


/* The filename should have an extension which will allow me to determine
 * the file type
 */
//---!!-------------------------------------------
  Observations::Observations(void) : NaN(sqrt(-1))
//---!!-------------------------------------------
  {
    using namespace std;
    string    fn;
    char      datafile[40],filename[120],directory[80];
    int       i,n,choice;
    Control   *control = Control::getInstance();

    //--- To ensure the right sampling frequency is read, fs=NaN at first
    fs=NaN;

    control->get_string("DataFile",datafile);
    control->get_string("DataDirectory",directory);
    scale_factor = control->get_double("ScaleFactor");
    strcpy(filename,directory);
    n = strlen(filename);
    if (filename[n-1]!='/') strcat(filename,"/");
    strcat(filename,datafile);
    fn = filename;
    i  = fn.find_last_of(".");
    strncpy(name,filename,i);
    strcpy(extension,&filename[i+1]);
    if (strcmp(extension,"mom")==0) {
      cout << "Data format: MJD, Observations, Model" << endl;
      read_observations = &Observations::read_mom;
    } else if (strcmp(extension,"enu")==0) {
      cout << "Data format: MJD, East, North & Up" << endl;
      read_observations = &Observations::read_enu;
    } else if (strcmp(extension,"neu")==0) {
      cout << "Data format: year fraction, North, East & Up" << endl;
      read_observations = &Observations::read_neu;
    } else if (strcmp(extension,"rlrdata")==0) {
      cout << "Data format: monthly PSMSL data" << endl;
      read_observations = &Observations::read_PSMSL_monthly;
    } else {
      cerr << "Cannot understand which type is : " << filename << endl;
      cout << "name=" << name << ", extension=" << extension << endl;
      exit(EXIT_FAILURE);
    }

    //--- Start with empty vectors
    t.clear();
    x.clear();
    xhat.clear();

    //--- read the data
    (*this.*read_observations)(filename);

    //--- For the handling of the gaps I need to know some things
    interpolate_data = control->get_bool("interpolate");
    first_difference = control->get_bool("firstdifference");

    //--- Need to interpolate the data?
    if (interpolate_data==true) {
      make_continuous(true);
    } else {
      //--- The Levinson & AmmarGrag method cannot deal with gaps in the
      //    time vector, only the FullCov method can. Therefore, fill gaps 
      //    with NaN's to make time-series continuous.
      make_continuous(false);
    }

    //--- Taking the first difference is necessary? 
    if (first_difference==true) {
      //--- Data is already continous
      take_first_difference();
    }

    //--- Now that I have read both the offset information and the data,
    //    I can remove trouble offsets (outside period, double counting)
    clean_offsets();

    //--- Count gaps    
    Ngaps=0;
    for (i=0;i<x.size();i++) {
      if (isnan(x[i])) Ngaps++;
    }

    //--- Give summary of read data
    cout << "Filename              : " << filename << endl;
    cout << "Number of observations: " << t.size() << endl;
    cout << "Percentage of gaps    : " << double(Ngaps)/t.size()*100.0 << endl;
  }


//-------------------------------------------------------
  long Observations::julday(int year, int month, int day)
//-------------------------------------------------------
/* Compute for given year, month and day the Julian day
 * 
 * Example: YEAR = 1970, MONTH = 1, DAY = 1, JD = 2440588. 
 * Reference: Fliegel, H. F. and van Flandern, T. C., Communications of 
 * the ACM, Vol. 11, No. 10 (October, 1968).
 *
 * http://www.usno.navy.mil/USNO/astronomical-applications/
 */
  {
    int       i,j,k;
    long int  jul;

    using namespace std;
    if (year<1801 || year>2099) {
      cerr << "year " << year << " is out of possible range!" << endl;
      exit(EXIT_FAILURE);
    }

    i = year;
    j = month;
    k = day;

    jul = k-32075+1461*(i+4800+(j-14)/12)/4+367*(j-2-(j-14)/12*12)/12 -
                                           3*((i+4900+(j-14)/12)/100)/4;
    return jul;
  }



/*! Read header information
 */
//------------------------------------------------------------------
  void Observations::read_header(std::fstream& fp, int component=-1)
//------------------------------------------------------------------
  {
    using namespace std;
    int            component_;
    char           line[80];
    double         MJD,T;

    offsets.clear();
    while (fp.peek() == '#') {
      fp.getline(line,80);
#ifdef DEBUG
      cout << "--->" << line << endl;
#endif
      if (strncmp("# offset",line,8)==0) {
        if (component==-1 && sscanf(&line[8],"%lf",&MJD)==1) {
          offsets.push_back(MJD);
        } else if (sscanf(&line[8],"%lf %d",&MJD,&component_)==2) {
          if (component==component_) offsets.push_back(MJD);
        } else {
          cerr << "Could not understand offset-line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      } else if (strncmp("# sampling period",line,17)==0) {
        if (sscanf(&line[17],"%lf",&T)==1) {
          if (!isnan(fs)) {
            if (fabs(fs*24.0*3600.0 - T)>0.001) {
              cerr << "Weird, fs was set but different to header info" << endl
                   << "T=" << T << ",  fs=" << fs << endl;
              exit(EXIT_FAILURE);
            }
          }
          fs = 1.0/(T*24.0*3600.0); // sampling frequency in Hz
        } else {
          cerr << "Could not understand sampling period-line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }


/*! Clean up any found offsets. I mean, check if they fall between
 *  start and end time of data and if there is not more than one
 *  offset in a sequence of missing data.
 */
//--------------------------------------
  void Observations::clean_offsets(void)
//--------------------------------------
  {
    using namespace std;
    const double  TINY=1.0e-5;
    bool          index_already_used;
    int           i,j,k,m=t.size();
    vector<int>   index;

    i=0;
    index.clear();
#ifdef DEBUG
    for (j=0;j<offsets.size();j++)
      cout << j << " : offset=" << offsets[j] << endl;
#endif
    while (i<offsets.size()) {
      if (offsets[i]>t[0] && offsets[i]<t[m-1]) {
        //--- find index
        j=0;
        while (t[j]<offsets[i] || (t[j]>offsets[i]-TINY && isnan(x[j]))) j++;
        //cout << "i=" << i << ", index=" << j << " t[j]=" << t[j] <<  endl;
        index_already_used = false;
        for (k=0;k<index.size();k++) {
          if (index[k]==j) {
            index_already_used = true;
          }
        }
        if (index_already_used==true) {
          cout << "offset " << offsets[i] << " is already used" << endl;
          offsets.erase(offsets.begin()+i);
        } else {
          index.push_back(j);
          //cout << index.size() << " : " << index[index.size()-1] << endl;
          i++;
        }
      } else {
        cout << "offset " << offsets[i] << " is outside time span" << endl;
        cout << t[0] << ", " << t[m-1] << endl;
        offsets.erase(offsets.begin()+i);
      }
    }
#ifdef DEBUG
    for (j=0;j<offsets.size();j++)
      cout << j << " : offset=" << offsets[j] << endl;
#endif

    //--- If first-difference of data is taken, then treat those observations
    //    as outliers.
    if (first_difference==true) {
      for (i=0;i<offsets.size();i++) {
        for (j=1;j<m;j++) {
          if (t[j-1]<offsets[i] && t[j]+TINY>offsets[i]) {
            set_one_x(j,NaN);
          }
        }
      }
    }
  }



/*! After removing the outliers, I need to create a new data file that includes
 *  the header information (sampling-period and offsets).
 */
//-------------------------------------------------
  void Observations::write_header(std::fstream& fp)
//-------------------------------------------------
  {
    int    i;

    using namespace std;
    fp << "# sampling period " << 1.0/(fs*24.0*3600.0) << endl;
    for (i=0;i<offsets.size();i++) {
      fp << "# offset " << offsets[i] << endl;
    }
  }



//------------------------------------------------------
  void Observations::read_PSMSL_monthly(char filename[])
//------------------------------------------------------
  {
    using namespace std;
    fstream   fp;
    string    fn(filename);
    double    fraction,MJD;
    int       MSL,flag,year,yearly,missing;
    char      line[80],complete_filename[120],fileout[120];
    long int  J;

    //--- Open file
    fs = 1.0/(30.4375*24.0*3600.0); // sampling frequency (Hz)
    fp.open(filename,ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << complete_filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Get first date. The rest is just adding sampling periods to MJD
    fp.getline(line,80);
    if (sscanf(line,"%lf;%d;%d;%d",&fraction,&MSL,&missing,&flag)==4) {
      year = static_cast<int>(floor(fraction));
      J    = julday(year,1,0);
      MJD  = static_cast<double>(J-2400001) + fmod(fraction,1.0)*365.25;
    } else {
      cerr << "Direct trouble reading the first line!" << endl;
      exit(EXIT_FAILURE);
    }
    //--- Read the whole file
    do {
      if (sscanf(line,"%lf;%d;%d;%d",&fraction,&MSL,&missing,&flag)==4) {
        if (flag==0 && MSL!=-99999) {
          t.push_back(MJD);
          x.push_back(static_cast<double>(scale_factor*MSL)); // scale data
        } else {
          cout << "found gap at: " << line << endl;
        }
        MJD += 30.4375; // I assume the PSMSL files have no gaps
      } else {
        cerr << "Unable to understand line: " << line;
        exit(EXIT_FAILURE);
      }
    } while (fp.getline(line,80)!=NULL);

    //--- Close file
    fp.close();
  }



//--------------------------------------------
  void Observations::read_mom(char filename[])
//--------------------------------------------
  {
    using namespace std;
    fstream    fp;
    char       line[80];
    double     MJD,obs,mod;

    //--- Open file
    fp.open(filename,ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Read header information if it exists
    read_header(fp);

    //--- Read file
    while (fp.getline(line,80)!=NULL) {
      if (sscanf(line,"%lf %lf %lf",&MJD,&obs,&mod)==3) {
#ifdef DEBUG
        cout << "Found x & xhat: " << obs << ", " << mod << endl;
#endif
        t.push_back(MJD);
        if (isnan(obs) || isnan(mod)) {
          cerr << "Found a NaN! : obs=" << obs << ", mod=" << mod << endl;
          exit(EXIT_FAILURE);
        } else {
          x.push_back(scale_factor*(obs-mod));
        }
      } else if (sscanf(line,"%lf %lf",&MJD,&obs)==2) {
#ifdef DEBUG
        cout << "Found x " << obs << endl;
#endif
        t.push_back(MJD);
        x.push_back(scale_factor*obs);
      } else {
        cerr << "Unable to understand line: " << line;
        exit(EXIT_FAILURE);
      }
    }

    //--- Close file
    fp.close();

    //--- Determine sampling period
    if (isnan(fs)) determine_fs();
  }



/*! See if I can determine the sampling period
 */
//-------------------------------------
  void Observations::determine_fs(void)
//-------------------------------------
  {
    int      i;
    double   dt;

    using namespace std;
    if (t.size()<2) {
      cerr << "Pointless exercise, I refuse..." << endl;
      exit(EXIT_FAILURE);
    } else {
      fs = 0.0;
      i  = 0;
      while (fs<1.0e-8 && i<10) {
        dt = t[i+1]-t[i];
        if      (fabs(48.0*dt-1.0)<1.0e-3)   fs = 1.0/(1800.0);
        else if (fabs(24.0*dt-1.0)<1.0e-3)   fs = 1.0/(3600.0);
        else if (fabs(dt-1.0)<1.0e-3)        fs = 1.0/(24.0*3600.0);
        else if (fabs(dt-7.0)<1.0e-3)        fs = 1.0/(7.0*24.0*3600.0);
        i++;
      }
      if (fs<1.0e-8) {
        cerr << "Cannot determine the sampling period! dt=" << dt << endl;
        exit(EXIT_FAILURE);
      } 
    }
#ifdef DEBUG
    cout << "fs=" << fs << endl;
#endif
  }



/*! Read my own brew GPS data format containing MJD, East, North and Up
 */
//--------------------------------------------
  void Observations::read_enu(char filename[])
//--------------------------------------------
  {
    using namespace std;
    fstream    fp;
    const double  TINY=1.0e-4;
    bool          first_line=true;
    char          line[80],component_name[20];
    int           component;
    double        MJD,obs[3],MJD_old;
    Control       *control = Control::getInstance();

    //--- Which component needs to be analysed?
    control->get_string("component",component_name);
    if      (strcmp(component_name,"East")==0)  component=0;
    else if (strcmp(component_name,"North")==0) component=1;
    else if (strcmp(component_name,"Up")==0)    component=2;
    else {
      cerr << "Unknown component name: " << component_name << endl;
      exit(EXIT_FAILURE);
    }

    //--- Open file
    fp.open(filename,ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Read header information if it exists
    read_header(fp,component);

    //--- Read file
    while (fp.getline(line,80)!=NULL) {
      if (sscanf(line,"%lf%lf%lf%lf",&MJD,&obs[0],&obs[1],&obs[2])==4) {
        if (first_line==true || MJD-MJD_old>TINY) {
          t.push_back(MJD);
          x.push_back(scale_factor*obs[component]);
          MJD_old = MJD;
          first_line=false;
        } else {
          cout << "Problem: " << line << endl;
        }
      } else {
        cerr << "Unable to understand line: " << line;
        exit(EXIT_FAILURE);
      }
    }

    //--- Close file
    fp.close();

    //--- Determine sampling period
    if (isnan(fs)) determine_fs();
  }



/* JPL, Simon and Guy Woppelman prefer North,East Up...
 */
//--------------------------------------------
  void Observations::read_neu(char filename[])
//--------------------------------------------
  {
    using namespace std;
    fstream       fp;
    const double  TINY=1.0e-4;
    bool          first_line=true,found_fs=false;
    char          line[120],component_name[20];
    int           decimal,remainder,component,i;
    double        MJD,obs[3],MJD_old,yearfraction,dt;
    Control       *control = Control::getInstance();

    //--- Which component needs to be analysed?
    control->get_string("component",component_name);
    if      (strcmp(component_name,"East")==0)  component=1;
    else if (strcmp(component_name,"North")==0) component=0;
    else if (strcmp(component_name,"Up")==0)    component=2;
    else {
      cerr << "Unknown component name: " << component_name << endl;
      exit(EXIT_FAILURE);
    }

    //--- Open file
    fp.open(filename,ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Read header information if it exists
    offsets.clear();
    while (fp.peek() == '#') {
      fp.getline(line,120);
#ifdef DEBUG
      cout << "--->" << line << endl;
#endif
      if (strncmp("#Offset",line,7)==0 || strncmp("# offset",line,8)==0) {
        if (sscanf(&line[8],"%lf %d",&yearfraction,&decimal)==2) {
          for (i=0;i<=component;i++) {
            remainder = decimal % 2;
            decimal  /= 2;
          }
          if (remainder==1) { 
            MJD = floor(365.25*(yearfraction-1970.0)+40587.0+0.1)-0.5; 
            offsets.push_back(MJD);
            cout << "offset at MJD=" << MJD << endl;
          }
        } else {
          cerr << "Could not understand offset-line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      }
    }

    //--- Read file
    while (fp.getline(line,120)!=NULL) {
      if (sscanf(line,"%lf%lf%lf%lf",&yearfraction,&obs[0],&obs[1],
							   &obs[2])==4) {
        MJD = floor(365.25*(yearfraction-1970.0)+40587.0+0.1)-0.5; 
        if (first_line==true || MJD-MJD_old>TINY) {
          t.push_back(MJD);
          x.push_back(scale_factor*obs[component]);
          if (found_fs==false && first_line==false) {
            dt = MJD-MJD_old;
            cout << "dt=" << dt << endl;
            if (fabs(dt-1.0)<0.1) {
              fs = 1.0/(24.0*3600.0); 
              found_fs=true;
            } else if (fabs(dt-7.0)<0.1) {
              fs = 1.0/(7.0*24.0*3600.0); 
              found_fs=true;
            }
          }  
          MJD_old = MJD;
          first_line=false;
        } else {
          cout << "Problem: " << line << endl;
        }
      } else {
        cerr << "Unable to understand line: " << line;
        exit(EXIT_FAILURE);
      }
    }

    //--- Close file
    fp.close();

    //--- Determine sampling period
    if (isnan(fs)) {
      cerr << "Could not determine sampling period!" << endl;
      exit(EXIT_FAILURE);
    } 
  }



/* For plotting purposes it might be interesting to save the observations in
 * MJD, Observations format
 */
//----------------------------------------
  void Observations::save_mom(bool header)
//----------------------------------------
  {
    using namespace std; 
    int      i;
    fstream  fp;
    char     filename[80];
    Control  *control=Control::getInstance();

    control->get_string("OutputFile",filename);

    cout.precision(5);
    cout << "--> " << filename << endl;
    fp.open(filename,ios::out);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    //--- In some cases, the header of the observation file needs to be
    //    copied to the output file
    if (header==true) write_header(fp);

    //--- Save the data to the output file
    for (i=0;i<t.size();i++) {
      if (xhat.size()==0) { 
        fp << fixed << t[i] << "  " << x[i] << endl;
      } else {
        if (!isnan(x[i])) {
          fp << fixed << t[i] << "  " << x[i] << "  " << xhat[i] << endl;
        }
      }
    }
    fp.close(); 
  }



/*! Copy computed xhat into the Observation class
 */
//------------------------------------------------
  void Observations::set_xhat(const double *xhat_)
//------------------------------------------------
  {
    int    i;

    using namespace std;
    for (i=0;i<x.size();i++) {
      xhat.push_back(xhat_[i]);
    }
  }



/* Although vectors are more beautiful creates than arrays, I can only
 * use arrays in LAPACK subroutines so will stick to simple arrays.
 */
//---------------------------------------------------------------
  void Observations::get_values(int& m, double **t_, double **x_)
//---------------------------------------------------------------
  { 
    int i;
   
    m = t.size();
    (*t_) = new double[m];
    (*x_) = new double[m];
    for (i=0;i<m;i++) {
      (*t_)[i] = t[i];
      (*x_)[i] = x[i];
    }
  }



/* For the removal of outliers I need to be able to set certain points to
 * NaN. It can also be used for other purposes.
 */
//-----------------------------------------------------
  void Observations::set_one_x(int index, double value)
//-----------------------------------------------------
  {
    x[index] = value;
    if (isnan(value)) Ngaps++;
  }



//----------------------------------------------------
  void Observations::make_continuous(bool interpolate)
//----------------------------------------------------
  {
    using namespace std;
    int             i,j;
    double          dt,slope;
    vector<double>  t_new,x_new;

#ifdef DEBUG
    cout << "before  m: " << t.size() << endl;
#endif

    //--- dt : sampling period (unit s -> day) 
    dt = 1.0/(24.0*3600.0*fs);

    //--- Create two new vectors
    t_new.clear();
    x_new.clear();

    //--- First point cannot be interpolated
    t_new.push_back(t[0]);
    x_new.push_back(x[0]);

    //--- Copy each row and interpolate when there's a gap
    j=0;
    cout.precision(4);
    for (i=1;i<t.size();i++) {
      if (fabs(t[i]-t[i-1]-dt)>0.1*dt) {
        slope = (x[i]-x[i-1])/(t[i]-t[i-1]);
        do {
          t_new.push_back(t_new[j] + dt);
          if (interpolate==true) {
            x_new.push_back(x[i-1] + slope*((t_new[j]+dt)-t[i-1]));
          } else {
            x_new.push_back(NaN);
          }
          j++;
        } while ((t[i]-t_new[j]-dt)>0.1*dt);
        if (fabs(t[i]-t_new[j]-dt)>0.1*dt) {
          cerr << "Oops, the data is not equally spaced over the missing data!"
               << endl << "MJD=" << t[i] << endl;
          exit(EXIT_FAILURE);
        }
      }
      t_new.push_back(t[i]);
      x_new.push_back(x[i]);
      j++;
    }

    //---- Old t and x values are no longer necessary
    t.clear();
    x.clear();

    //---- Copy t_new to t and x_new to x
    for (i=0;i<t_new.size();i++) t.push_back(t_new[i]);
    for (i=0;i<x_new.size();i++) x.push_back(x_new[i]);
#ifdef DEBUG
    cout << "after  m: " << t.size() << endl;
#endif
  }



/*! For the direct method there should be no NaN's in the data. These
 *  rows should be eliminated.
 */
//------------------------------------
  void Observations::remove_gaps(void)
//------------------------------------
  {
    using namespace std;
    int             i;
    vector<double>  t_new,x_new;

    //--- Create two new vectors, without NaN's
    t_new.clear();
    x_new.clear();

    for (i=0;i<t.size();i++) {
      if (!isnan(x[i])) {
        t_new.push_back(t[i]);
        x_new.push_back(x[i]);
      }
    }

    //---- Old t and x values are no longer necessary
    t.clear();
    x.clear();

    //---- Copy t_new to t and x_new to x
    for (i=0;i<t_new.size();i++) t.push_back(t_new[i]);
    for (i=0;i<x_new.size();i++) x.push_back(x_new[i]);
#ifdef DEBUG
    cout << "New size after removal of gaps: m=" << t.size() << endl;
#endif
    //--- Count gaps    
    Ngaps=0;
    for (i=0;i<x.size();i++) {
      if (isnan(x[i])) Ngaps++;
    }
  }



/* If the data is non-stationary, then taking the first-difference may
 * make it stationary. Note that you have to be carefull with gaps.
 */
//----------------------------------------------
  void Observations::take_first_difference(void)
//----------------------------------------------
  {
    using namespace std;
    int             i;
    vector<double>  t_new,x_new;

    //--- Start with empty new vectors
    t_new.clear();
    x_new.clear();

    //--- Take first difference
    for (i=1;i<t.size();i++) {
      t_new.push_back(t[i]);
      //--- I hope the following works when x[i] or x[i-1] is NaN.
      x_new.push_back(x[i]-x[i-1]);
    }

    //---- Old t and x values are no longer necessary
    t.clear();
    x.clear();

    //---- Copy t_new to t and x_new to x
    for (i=0;i<t_new.size();i++) t.push_back(t_new[i]);
    for (i=0;i<x_new.size();i++) x.push_back(x_new[i]);
  }



/*! Copy offsets to offsets_
 */
//-------------------------------------------------------------
  void Observations::get_offsets(std::vector<double>& offsets_)
//-------------------------------------------------------------
  {
    int   i;

    offsets_.clear();
    for (i=0;i<offsets.size();i++) {
      offsets_.push_back(offsets[i]);
    }
  }



//--!!-----------------------------------------
  Observations* Observations::getInstance(void)
//--!!-----------------------------------------
  {
    if (instanceFlag==false) {
      singleton = new Observations();
      instanceFlag=true;
    }
    return singleton;
  }      
