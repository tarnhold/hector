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
  #include <fstream>
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
    strcpy(ctl_file,file.c_str());
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
//--------------------------------------------
  void Control::find_label(const char label[])
//--------------------------------------------
  {
    int   i,m;
    char  c,line[100],*retval;

    using namespace std;
    m  = strlen(label);
    rewind(fp_in);
    do {
      i=0; while((c=fgetc(fp_in))!=EOF && c!='\n' && c!=' '){line[i]=c;i++;} 
      line[i]='\0';
      if (i==m && strncmp(line,label,m)==0) return;       // if found...
      retval = fgets(line,100,fp_in); // else skip rest of line
    } while (c!=EOF);    
    cerr << "Could not find label: " << label << endl;
    exit(EXIT_FAILURE);
  }



/*! Find "label" and read a list of "n" strings into "value". Memory for 
 *  "value" must have been allocated before. 
 *
 * \param[in]  label[]: string with parameter name
 * \param[out] value: array of n *chars filled with values 
 * \param[out] n: the number of values read
 */
//---------------------------------------------------------------------
  void Control::get_name_list(const char label[], char **value, int& n)
//---------------------------------------------------------------------
  {
    int      i1,i2,m;
    char     line[150],*retval;

    find_label(label);
    retval = fgets(line,150,fp_in); // read the line containing the values
    m  = strlen(line)-1;   // -1 because of removing '\n' in count
    i1 = 0;
    while (i1<m && line[i1]==' ') i1++; // skip spaces
    i2 = i1;
    n  = 0; // sofar no single value has been read
    while (i2<m) { // look for values as long as the line contains characters
      while (i2<m && !(line[i2]==' ')) i2++; // find end of value
      value[n] = new char[i2-i1+1];
      strncpy(value[n],&line[i1],i2-i1);
      value[n][i2-i1]='\0';
      n++;
      i1=i2; // move i1 to end of value which is at i2
      while (i1<m && line[i1]==' ') i1++; // skip spaces
      i2=i1;
    }
  }



/*! Read a string
 *
 * \param[in]  label[]: name of parameter
 * \param{out] value[]: the read value (string)
 */
//----------------------------------------------------------
  void Control::get_string(const char label[], char value[])
//----------------------------------------------------------
  {
    int     i1,i2,m;
    char    line[100],*retval;

    using namespace std;
    find_label(label);
    retval = fgets(line,100,fp_in); // read the line containing the values
    m  = strlen(line);
    i1 = 0;
    while (i1<m && line[i1]==' ') i1++; // skip spaces
    if (i1==m) {
      cerr << "Not found a string for label " << label << endl;
      exit(EXIT_FAILURE);
    }
    i2 = i1;
    while (i2<m && !(line[i2]==' ' || line[i2]=='\n')) i2++; // find end 
    strncpy(value,&line[i1],i2-i1);
    value[i2-i1]='\0';
  }



/*! Read a boolean parameter
 *
 * \param[in]  label: name of parameter (string)
 * \return boolean (yes=true, no=false)
 */
//------------------------------------------
  bool Control::get_bool(const char label[])
//------------------------------------------
  {
    char   line[100],answer[100],*retval;

    using namespace std;
    find_label(label);
    retval = fgets(line,100,fp_in);
    sscanf(line,"%s",answer);
    if      (answer[0]=='y' || answer[0]=='Y') return true;
    else if (answer[0]=='n' || answer[0]=='N') return false;
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
//----------------------------------------
  int Control::get_int(const char label[])
//----------------------------------------
  {
    int   answer;

    using namespace std;
    find_label(label);
    if (fscanf(fp_in,"%d",&answer)>0) return answer;
    else {
      cerr << "Not an int answer for label " << label << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Read double value from file
 * 
 * \param[in]  label: name of parameter (string)
 * \return double that is read
 */
//----------------------------------------------
  double Control::get_double(const char label[])
//----------------------------------------------
  {
    double   answer;

    using namespace std;
    find_label(label);
    if (fscanf(fp_in,"%lf",&answer)>0) return answer;
    else {
      cerr<<"Not a double answer for label " << label << endl;
      exit(EXIT_FAILURE);
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
    fp_in = fopen(ctl_file,"rt");
    if (fp_in==NULL) {
      cerr << "Could not open " << ctl_file << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! close the ctl file */
//----------------------------------
  void Control::close_file(void)
//----------------------------------
  {
    if (fp_in!=NULL) {
      fclose(fp_in);
      fp_in = NULL;
    }
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
