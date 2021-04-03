/*! \file    Control.cpp
 *  \author  Machiel Bos
 *
 * Facilitate the sharing of settings between classes
 * 
 * \date 12/1/2011  CIIMAR, Porto
 */
//============================================================================
  #include <cstdio>
  #include <cstdlib>
  #include <cstring>
  #include <string>
  #include <iostream>
  #include <ostream>
  #include <string>
  #include "Control.h"

//============================================================================
// Subroutines
//============================================================================


/*! Open control file 
 *
 * \param[in] file: filename of ctl-file, by default this is EstimateTrend.ctl 
 */
//--!!------------------------------
  Control::Control(std::string file)
//--!!------------------------------
  {
    using namespace std;
    ctl_file = file;
    open_file();
  }



/*! close the ctl files */
//--!!-------------------
  Control::~Control(void)
//--!!-------------------
  {
    close_file();
  }



/*! Rewind the file and read until line with line has been found and
 *  return. Otherwise exit program with warning.
 *
 * \param[in] label: name of the parameter
 * \return The file-pointer is now located just before the value */
//-------------------------------------------------
  void Control::find_label(const std::string label)
//-------------------------------------------------
  {
    using namespace std;
    string       str,dummy,error_message;

    fs.clear();
    fs.seekg (0,fs.beg); // start at beginning of Control file
    do {
      fs >> str;
      if (fs.eof()!=true) {
        if (str.compare(label)==0) {
          return; // file pointer now starts at right place
        } else {
          getline(fs,dummy); // skip unwanted values in rest of line
        }
      }
    } while (fs.eof()!=true);
    fs.clear();
    error_message = "Could not find label:" + label;
    throw(const_cast<const char *>(error_message.c_str()));
  }



/*! Find "label" and read a list of "n" strings into "value". Memory for 
 *  "value" must have been allocated before. 
 *
 * \param[in]  label[]: string with parameter name
 * \param[out] value: array of n *chars filled with values 
 * \param[out] n: the number of values read
 */
//--------------------------------------------------------------------------
  void Control::get_name_list(const std::string label, std::string *value, 
								     int& n)
//--------------------------------------------------------------------------
  {
    using namespace std;
    string        str,line;
    stringstream  fs_line;

    find_label(label); // point file pointer to the right location
    getline(fs,line);  // read the whole line into string 'line'
    fs_line << line;   // copy line into the string-stream
    n=0;
    while (fs_line.good()==true && fs_line.eof()!=true) {
      fs_line >> str;
      if (str.length()>0) {
        value[n] = str;
        n++;
      }
    }
  }



/*! Read a string
 *
 * \param[in]  label[]: name of parameter
 * \param{out] value[]: the read value (string)
 */
//---------------------------------------------------------------------
  void Control::get_string(const std::string label, std::string& value)
//---------------------------------------------------------------------
  {
    using namespace std;
    find_label(label); // point file pointer to the right location
    fs >> value;
    if (value.size()==0) {
      cerr << "Not found a string for label " << label << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Read a boolean parameter
 *
 * \param[in]  label: name of parameter (string)
 * \return boolean (yes=true, no=false)
 */
//-----------------------------------------------
  bool Control::get_bool(const std::string label)
//-----------------------------------------------
  {
    using namespace std;
    string   answer;

    find_label(label);
    fs >> answer;
    if      (answer.compare("Yes")==0 || answer.compare("yes")==0) return true;
    else if (answer.compare("No")==0  || answer.compare("no")==0)  return false;
    else {
      cerr << "Not a bool answer for label " << label << ":" << answer << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Read integer value
 * 
 * \param[in]  label: name of parameter (string)
 * \return integer that is read
 */
//---------------------------------------------
  int Control::get_int(const std::string label)
//---------------------------------------------
  {
    int   i;

    using namespace std;
    find_label(label);
    fs >> i;
    if (fs.fail()==true) {
      cerr << "Not an int answer for label " << label << endl;
      exit(EXIT_FAILURE);
    } else {
      return i;
    }
  }



/*! Read double value from file
 * 
 * \param[in]  label: name of parameter (string)
 * \return double that is read
 */
//---------------------------------------------------
  double Control::get_double(const std::string label)
//---------------------------------------------------
  {
    double   y;

    using namespace std;
    find_label(label);
    fs >> y;
    if (fs.fail()==true) {
      cerr << "Not an double answer for label " << label << endl;
      exit(EXIT_FAILURE);
    } else {
      return y;
    }
  }



/*! Open ctl_file (normally = EstimateTrend.ctl) and remember the file 
 *  pointer.
 */
//-----------------------------
  void Control::open_file(void)
//-----------------------------
  {
    using namespace std;
    fs.open(ctl_file.c_str(),fstream::in);
    if (fs.is_open()!=true) {
      cerr << "Could not open " << ctl_file << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! close the ctl file */
//----------------------------------
  void Control::close_file(void)
//----------------------------------
  {
    if (fs.is_open()==true)  fs.close();
  }



/* In my opinion this subroutine is the core trick of the singleton
 * implementation. A pointer is made static and points to the instance.
 */
//---------------------------------------------------
  Control* Control::getInstance(std::string filename)
//---------------------------------------------------
  {
    static Control singleton(filename);

    return &singleton;
  }
