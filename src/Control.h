/*! \file    Control.h
 *  \author  Machiel Bos
 *
 * A small class for reading the parameters that control the behaviour of
 * EstimateTrend. 
 *
 * \date 12/1/2012   CIIMAR, Porto
 */
//============================================================================

//============================================================================
// Class definition
//============================================================================

#ifndef __CONTROL
  #define __CONTROL
  #include <exception>
  #include <stdexcept>
  #include <cstdio>
  #include <string>
  #include <fstream>
  #include <sstream>

  /*! Following the Namelist tradition of FORTRAN, this class lets you define
   *  parameter values in a ASCII file which is read at run-time. 
   */
  class Control 
  {
    private:
      std::fstream   fs;
      std::string    ctl_file;
      void           find_label(const std::string label);
      void           open_file(void);
      void           close_file(void);

    public:
      Control(std::string file);
      ~Control(void);
      static Control* getInstance(std::string filename = "estimatetrend.ctl");
      void     get_name_list(const std::string label, std::string *value, 
									int& n);
      void     get_string(const std::string label, std::string& value);
      bool     get_bool(const std::string label);
      int      get_int(const std::string label);
      double   get_double(const std::string label);
  };

#endif
