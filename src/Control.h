/*! \file    Control.h
 *  \author  Machiel Bos
 *
 * Header file of Control.cpp
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
      void           find_label(const std::string& label);
      void           open_file(void);
      void           close_file(void);

      explicit Control(const std::string& file);
      ~Control(void);

    public:
      static Control& getInstance(std::string filename = "estimatetrend.ctl");
      void     get_name_list(const std::string& label, std::string *value,
									int& n);
      void     get_string(const std::string& label, std::string& value);
      bool     get_bool(const std::string& label);
      int      get_int(const std::string& label);
      double   get_double(const std::string& label);
  };

#endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
