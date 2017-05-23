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
  #include "Observations.h"
  #include "Calendar.h"
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


/* The filename should have an extension which will allow me to determine
 * the file type
 */
//---!!-------------------------------------------
  Observations::Observations(void) : NaN(sqrt(-1))
//---!!-------------------------------------------
  {
    using namespace std;
    size_t    i;
    string    datafile,filename,directory;
    Control   &control = Control::getInstance();

    //--- To ensure the right sampling frequency is read, fs=NaN at first
    fs=NaN;

    //--- Avoid stupid errors
    offsets.clear();
    breaks.clear();
    postseismiclog.clear();
    postseismicexp.clear();

    //--- Get information on where to find the file with the observations
    try {
      control.get_string("DataFile",datafile);
      control.get_string("DataDirectory",directory);
      if (directory.at(directory.length()-1)!='/') {
        directory += "/";
      }
    } catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    //--- Is a scale factor different from 1 required?
    try {
      scale_factor = control.get_double("ScaleFactor");
    } catch (exception &e) {
      scale_factor = 1.0;
    }

    //--- Do we want to write empty mom records?
    try {
      write_empty_records = control.get_bool("WriteEmptyRecords");
    } catch (exception &e) {
      write_empty_records = false;
    }

    //--- construct filename
    filename = directory + datafile;
    
    //--- Try to guess the filename extension and split
    i  = filename.find_last_of(".");
    extension = filename.substr(i+1,filename.length()-i-1);
    
    if (extension.compare("mom")==0) {
      cout << "Data format: MJD, Observations, Model" << endl;
      read_observations = &Observations::read_mom;
    } else if (extension.compare("enu")==0) {
      cout << "Data format: MJD, East, North & Up" << endl;
      read_observations = &Observations::read_enu;
    } else if (extension.compare("neu")==0) {
      cout << "Data format: year fraction, North, East & Up" << endl;
      read_observations = &Observations::read_neu;
    } else if (extension.compare("rlrdata")==0) {
      cout << "Data format: monthly PSMSL data" << endl;
      read_observations = &Observations::read_PSMSL_monthly;
    } else {
      cerr << "Cannot understand which type is : " << filename << endl;
      cout << "extension=" << extension << endl;
      exit(EXIT_FAILURE);
    }

    //--- Start with empty vectors
    t.clear();
    x.clear();
    xhat.clear();

    //--- read the data
    (*this.*read_observations)(filename);

    //--- For the handling of the gaps I need to know some things
    try {
      interpolate_data = control.get_bool("interpolate");
    } catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    //--- Need to interpolate the data?
    if (interpolate_data==true) {
      make_continuous(true);
    } else {
      //--- The AmmarGrag method cannot deal with gaps in the
      //    time vector, only the FullCov method can. Therefore, fill gaps 
      //    with NaN's to make time-series continuous.
      make_continuous(false);
    }

    //--- Now that I have read both the offset information and the data,
    //    I can remove trouble offsets (outside period, double counting)
    clean_offsets();

    //--- Count gaps    
    Ngaps=0;
    for (i=0;i<x.size();i++) {
      if (std::isnan(x[i])) Ngaps++;
    }

    //--- Give summary of read data
    cout << "Filename              : " << filename << endl;
    cout << "Number of observations: " << t.size()-Ngaps << endl;
    cout << "Number of gaps        : " << Ngaps << endl;
    cout << "Percentage of gaps    : " << double(Ngaps)/t.size()*100.0 << endl;
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
    LogEntry       logdummy;
    ExpEntry       expdummy;
    TanhEntry      tanhdummy;

    breaks.clear();
    offsets.clear();
    postseismiclog.clear();
    postseismicexp.clear();
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
      } else if (strncmp("# break",line,7)==0) {
        if (sscanf(&line[8],"%lf",&MJD)==1) {
          if (breaks.size()==0 || 
		(breaks.size()>0 && breaks[breaks.size()-1]<MJD)) {
            breaks.push_back(MJD);
          } else {
            cerr << "Please list breaks in chronological order!" << endl;
            cerr << line << endl;
            exit(EXIT_FAILURE);
          }
        } else {
          cerr << "Could not understand break-line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      } else if (strncmp("# sampling period",line,17)==0) {
        if (sscanf(&line[17],"%lf",&T)==1) {
          if (!std::isnan(fs)) {
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
      } else if (strncmp("# log",line,5)==0) {
        if (sscanf(&line[5],"%lf %lf",&logdummy.MJD,&logdummy.T)==2) {
          postseismiclog.push_back(logdummy);
        } else {
          cerr << "Could not understand postseismic-log line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      } else if (strncmp("# exp",line,5)==0) {
        if (sscanf(&line[5],"%lf %lf",&expdummy.MJD,&expdummy.T)==2) {
          postseismicexp.push_back(expdummy);
        } else {
          cerr << "Could not understand postseismic-exp line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      } else if (strncmp("# tanh",line,6)==0) {
        if (sscanf(&line[6],"%lf %lf",&tanhdummy.MJD,&tanhdummy.T)==2) {
          ssetanh.push_back(tanhdummy);
        } else {
          cerr << "Could not understand SSE-tanh line:" << line << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }



/*! Read header information from other file
 */
//---------------------------------------------------------
  void Observations::read_external_header(std::fstream& fp)
//---------------------------------------------------------
  {
    using namespace std;
    int            component,day,month,year;
    char           line[80];
    double         MJD,comp[3];
    string         component_name;
    Calendar       calendar;
    Control        &control=Control::getInstance();

    //--- We need to know which component is requested.
    try {
      control.get_string("component",component_name);
    } catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
    if      (component_name.compare("East")==0)  component=0;
    else if (component_name.compare("North")==0) component=1;
    else if (component_name.compare("Up")==0)    component=2;
    else {
      cerr << "Unknown component name: " << component_name << endl;
      exit(EXIT_FAILURE);
    }

    //--- Clear offsets and start reading the lines
    offsets.clear();
    while (!fp.eof() ) {
      fp.getline(line,80);
      if (sscanf(line,"%d-%d-%d %lf %lf %lf",&day,&month,&year,&comp[0],
						    &comp[1],&comp[2])==6) {
        MJD = static_cast<double>(calendar.julday(year,month,day)-2400001);
        //--- Mom files don't come with component information so must
        //    look for the control file to know which component we need.
        cout << "MJD=" << MJD << ", " << comp[0] << ", " << comp[1] 
						<< ", " << comp[2] << endl;
        if (std::isnan(comp[component])) {
          offsets.push_back(MJD);
        }
      } else if (strlen(line)>0) {
        cerr << "Could not understand line:" << line << endl;
        exit(EXIT_FAILURE);
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
    const double   TINY=1.0e-5;
    bool           index_already_used;
    size_t         i,j,k,m=t.size();
    vector<size_t> index;

    i=0;
    index.clear();
#ifdef DEBUG
    for (j=0;j<offsets.size();j++)
      cout << j << " : offset=" << fixed << offsets[j] << endl;
#endif
    while (i<offsets.size()) {
      if (offsets[i]>t[0]-TINY && offsets[i]<t[m-1]+TINY) {
        //--- find index
        j=0;
        while (t[j]<offsets[i]-TINY || (t[j]>offsets[i]-TINY 
						&& std::isnan(x[j]))) j++;
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
        cout << fixed
             << "offset " << offsets[i] << " is outside time span ["
             << t[0] << ", " << t[m-1] << "]" << endl;
        offsets.erase(offsets.begin()+i);
      }
    }
#ifdef DEBUG
    for (j=0;j<offsets.size();j++)
      cout << j << " : offset=" << fixed << offsets[j] << endl;
#endif

  }



/*! After removing the outliers, I need to create a new data file that includes
 *  the header information (sampling-period and offsets).
 */
//-------------------------------------------------
  void Observations::write_header(std::fstream& fp)
//-------------------------------------------------
  {
    size_t     i;
    LogEntry   logdummy;
    ExpEntry   expdummy;
    TanhEntry  tanhdummy;

    using namespace std;
    fp << "# sampling period " << 1.0/(fs*24.0*3600.0) << endl;
    for (i=0;i<offsets.size();i++) {
      fp << "# offset " << offsets[i] << endl;
    }
    for (i=0;i<breaks.size();i++) {
      fp << "# break " << breaks[i] << endl;
    }
    for (i=0;i<postseismiclog.size();i++) {
      logdummy = postseismiclog[i];
      fp << "# log " << logdummy.MJD << "  " << logdummy.T << endl;
    }
    for (i=0;i<postseismicexp.size();i++) {
      expdummy = postseismicexp[i];
      fp << "# exp " << expdummy.MJD << "  " << expdummy.T << endl;
    }
    for (i=0;i<ssetanh.size();i++) {
      tanhdummy = ssetanh[i];
      fp << "# tanh " << tanhdummy.MJD << "  " << tanhdummy.T << endl;
    }
  }



//-----------------------------------------------------------
  void Observations::read_PSMSL_monthly(std::string filename)
//-----------------------------------------------------------
  {
    using namespace std;
    fstream   fp;
    double    fraction,MJD;
    int       MSL,flag,year,missing,month;
    char      line[80];

    //--- Open file
    fs = 1.0/(30.4375*24.0*3600.0); // sampling frequency (Hz)
    fp.open(filename.c_str(),ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Get first date. The rest is just adding sampling periods to MJD
    fp.getline(line,80);
    if (sscanf(line,"%lf;%d;%d;%d",&fraction,&MSL,&missing,&flag)==4) {
      year = static_cast<int>(floor(fraction));
      month= static_cast<int>(0.501 + (fraction-year)*12.0);
      MJD = 30.4375*(12*(year-1859) + (month-1.0)) + 59.0;
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
          //cout << "found gap at: " << line << endl;
        }
        MJD += 30.4375; // I assume the PSMSL files have no gaps
      } else {
        cerr << "Unable to understand line: " << line;
        exit(EXIT_FAILURE);
      }
      fp.getline(line,80);
    } while (!fp.eof());

    //--- Close file
    fp.close();
  }



//-------------------------------------------------
  void Observations::read_mom(std::string filename)
//-------------------------------------------------
  {
    using namespace std;
    fstream    fp,fp_offset;
    string     offset_filename;
    char       line[80];
    double     MJD,obs,mod;
    Control    &control=Control::getInstance();

    //--- Open file
    fp.open(filename.c_str(),ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Read header information if it exists
    read_header(fp);
    try {
      control.get_string("OffsetFile",offset_filename);
      fp_offset.open(offset_filename.c_str(),ios::in);
      if (!fp_offset.is_open()) {
        cerr << "Could not open " << offset_filename << endl;
        exit(EXIT_FAILURE);
      }
      read_external_header(fp_offset);
      fp_offset.close();
    }
    catch (exception &e) {
      //--- No OffsetFile found, use data file instead
    }

    //--- Read file
    fp.getline(line,80);
    while (!fp.eof()) {
      if (sscanf(line,"%lf %lf %lf",&MJD,&obs,&mod)==3) {
#ifdef DEBUG
        cout << "Found x & xhat: " << fixed << setprecision(5) << obs << ", " << mod << endl;
#endif
        t.push_back(MJD);
        if (std::isnan(obs) || std::isnan(mod)) {
          cerr << "Found a NaN! : obs=" << obs << ", mod=" << mod << endl;
          exit(EXIT_FAILURE);
        } else {
          x.push_back(scale_factor*(obs-mod));
        }
      } else if (sscanf(line,"%lf %lf",&MJD,&obs)==2) {
#ifdef DEBUG
        cout << "Found x " << fixed << setprecision(5) << obs << endl;
#endif
        if (!std::isnan(obs)) {
          t.push_back(MJD);
          x.push_back(scale_factor*obs);
        }
      } else if (sscanf(line,"%lf",&MJD)==1) {
#ifdef DEBUG
        cout << "Found empty MJD line, ignoring it" << endl;
#endif
      } else {
        cerr << "Unable to understand line: " << line << endl;
        exit(EXIT_FAILURE);
      }
      fp.getline(line,80);
    }

    //--- Close file
    fp.close();
    

    //--- Determine sampling period
    if (std::isnan(fs)) determine_fs();
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
        else if (fabs(dt-30.4375)<1.0e-3)    fs = 1.0/(30.4375*24.0*3600.0);
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
//-------------------------------------------------
  void Observations::read_enu(std::string filename)
//-------------------------------------------------
  {
    using namespace std;
    fstream       fp,fp_offset;
    const double  TINY=1.0e-4;
    bool          first_line=true;
    string        component_name,offset_filename;
    char          line[80];
    int           component;
    double        MJD,obs[3],MJD_old = 0.0;
    Control       &control = Control::getInstance();

    //--- Which component needs to be analysed?
    try {
      control.get_string("component",component_name);
    } catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
    if      (component_name.compare("East")==0)  component=0;
    else if (component_name.compare("North")==0) component=1;
    else if (component_name.compare("Up")==0)    component=2;
    else {
      cerr << "Unknown component name: " << component_name << endl;
      exit(EXIT_FAILURE);
    }

    //--- Open file
    fp.open(filename.c_str(),ios::in);
    if (!fp.is_open()) {
      cerr << "Could not open " << filename << endl;
      exit(EXIT_FAILURE);
    }

    //--- Read header information if it exists
    read_header(fp,component);
    try {
      control.get_string("OffsetFile",offset_filename);
      fp_offset.open(offset_filename.c_str(),ios::in);
      if (!fp_offset.is_open()) {
        cerr << "Could not open " << offset_filename << endl;
        exit(EXIT_FAILURE);
      }
      read_external_header(fp_offset);
      fp_offset.close();
    }
    catch (exception &e) {
      //--- No OffsetFile found, use data file instead
    }

    //--- Read file
    fp.getline(line,80);
    while (!fp.eof()) {
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
      fp.getline(line,80);
    }

    //--- Close file
    if (fp_offset.is_open()) fp_offset.close();
    fp.close();

    //--- Determine sampling period
    if (std::isnan(fs)) determine_fs();
  }



/* JPL, Simon and Guy Woppelman prefer North,East Up...
 */
//-------------------------------------------------
  void Observations::read_neu(std::string filename)
//-------------------------------------------------
  {
    using namespace std;
    fstream       fp;
    const double  TINY=1.0e-4;
    bool          first_line=true,found_fs=false;
    char          line[120];
    string        component_name;
    int           decimal,remainder=0,component,i;
    double        MJD,obs[3],MJD_old = 0.0,yearfraction,dt;
    Control       &control = Control::getInstance();

    //--- Which component needs to be analysed?
    try {
      control.get_string("component",component_name);
    } catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
    if      (component_name.compare("East")==0)  component=1;
    else if (component_name.compare("North")==0) component=0;
    else if (component_name.compare("Up")==0)    component=2;
    else {
      cerr << "Unknown component name: " << component_name << endl;
      exit(EXIT_FAILURE);
    }

    //--- Open file
    fp.open(filename.c_str(),ios::in);
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
    fp.getline(line,120);
    while (!fp.eof()) {
      if (sscanf(line,"%lf%lf%lf%lf",&yearfraction,&obs[0],&obs[1],
							   &obs[2])==4) {
        MJD = floor(365.25*(yearfraction-1970.0)+40587.0+0.1)-0.5; 
          cout << "# yearfraction=" << yearfraction-2000 << ", MJD_old=" << MJD_old << ", MJD=" << MJD << endl;
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
        cout << ">> yearfraction=" << yearfraction-2000 << ", MJD_old=" << MJD_old << ", MJD=" << MJD << endl;
      } else {
        cerr << "Unable to understand line: " << line;
        exit(EXIT_FAILURE);
      }
      fp.getline(line,120);
    }

    //--- Close file
    fp.close();

    //--- Determine sampling period
    if (std::isnan(fs)) {
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
    size_t   i;
    fstream  fp;
    string   filename;
    Control  &control=Control::getInstance();

    try {
      control.get_string("OutputFile",filename);
    } catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
    cout.precision(5);
    cout << "--> " << filename << endl;
    fp.open(filename.c_str(),ios::out);
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
        fp << setprecision(6) << fixed << t[i] << "  "
           << x[i] << endl;
       } else {
         if (!std::isnan(x[i])) {
          fp << setprecision(6) << fixed << t[i] << "  "
             << x[i] << "  " << xhat[i] << endl;
         } else {
           if (write_empty_records) {
             // just write epoch with no data
             fp << setprecision(6) << fixed << t[i] << endl;
          }
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
    size_t i;

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
    if (std::isnan(value)) Ngaps++;
  }



//----------------------------------------------------
  void Observations::make_continuous(bool interpolate)
//----------------------------------------------------
  {
    using namespace std;
    size_t          i,j;
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
          cerr << "Oops, missing data timing problem at MJD=" << t[i] << endl;
          cerr << "t_new="<<t_new[j]<< ", dt="<<dt << endl;
          cerr << "There could be various reasons for this" << endl;
          cerr << "1) There is no sampling information in the header and "
               << "Hector guessed wrongly" << endl;
          cerr << "  Solution: add \"# sampling period\" in header!" << endl;
          cerr << "2) time is not monotonically increasing" << endl;
          cerr << "3) time is an integer of dt" << endl;
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
    size_t          i;
    vector<double>  t_new,x_new;

    //--- Create two new vectors, without NaN's
    t_new.clear();
    x_new.clear();

    for (i=0;i<t.size();i++) {
      if (!std::isnan(x[i])) {
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
      if (std::isnan(x[i])) Ngaps++;
    }
  }



/*! Copy offsets to offsets_
 */
//-------------------------------------------------------------
  void Observations::get_offsets(std::vector<double>& offsets_)
//-------------------------------------------------------------
  {
    size_t i;

    offsets_.clear();
    for (i=0;i<offsets.size();i++) {
      offsets_.push_back(offsets[i]);
    }
  }



/*! Copy breaks to breaks_
 */
//-----------------------------------------------------------
  void Observations::get_breaks(std::vector<double>& breaks_)
//-----------------------------------------------------------
  {
    size_t i;

    breaks_.clear();
    for (i=0;i<breaks.size();i++) {
      breaks_.push_back(breaks[i]);
    }
  }



/*! Copy postseismiclog to postseismiclog_
 */
//-----------------------------------------------------------------------------
  void Observations::get_postseismiclog(std::vector<LogEntry>& postseismiclog_)
//-----------------------------------------------------------------------------
  {
    size_t i;

    postseismiclog_.clear();
    for (i=0;i<postseismiclog.size();i++) {
      postseismiclog_.push_back(postseismiclog[i]);
    }
  }



/*! Copy postseismicexp to postseismicexp_
 */
//-----------------------------------------------------------------------------
  void Observations::get_postseismicexp(std::vector<ExpEntry>& postseismicexp_)
//-----------------------------------------------------------------------------
  {
    size_t i;

    postseismicexp_.clear();
    for (i=0;i<postseismicexp.size();i++) {
      postseismicexp_.push_back(postseismicexp[i]);
    }
  }



/*! Copy ssetanh to ssetanh_ 
 */
//----------------------------------------------------------------
  void Observations::get_ssetanh(std::vector<TanhEntry>& ssetanh_)
//----------------------------------------------------------------
  {
    size_t i;

    ssetanh_.clear();
    for (i=0;i<ssetanh.size();i++) {
      ssetanh_.push_back(ssetanh[i]);
    }
  }



/*! Add offset
 */
//-----------------------------------------
  void Observations::add_offset(double MJD)
//-----------------------------------------
  {
    offsets.push_back(MJD);
  }



/*! Add offset
 */
//--------------------------------------------------------
  void Observations::change_offset(size_t column, double MJD)
//--------------------------------------------------------
  {
    using namespace std;
    if (column>=offsets.size()) {
      cerr << "Invalid column specification: " << column << endl;
      exit(EXIT_FAILURE);
    }
    offsets[column] = MJD;
  }


/*! To compute prior, I need first and last epoch
 */
//---------------------------------------------------
  void Observations::get_t0t1(double& t0, double& t1)
//---------------------------------------------------
  {
    t0 = t[0]; 
    t1 = t[t.size()-1]; 
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
