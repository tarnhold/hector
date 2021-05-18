/*! \file    JSON.h
 *  \author  Machiel Bos
 *
 * Header file of JSON.cpp
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
//============================================================================

//============================================================================
// Class definition
//============================================================================

#ifndef __JSON
  #define __JSON
  #include <cstdio>
  #include <string>
  #include <fstream>
  #include <vector>
  #include "Control.h"

  class JSON 
  {
    private:
      std::fstream   fp;
      int            indent;
      bool           json_required;

      JSON(std::string file);
      ~JSON(void);

    public:
      static   JSON& getInstance(std::string filename = "estimatetrend.json");
      void     start_dictionary(const std::string label);
      void     end_dictionary(bool last);
      void     write_int(const std::string label, int i, bool last=false);
      void     write_double(const std::string label, double y, bool last=false);
      void     write_double_list(const std::string label, int n, double *y,
							       bool last=false);
      void     write_chararray_list(const std::string label, 
				   std::vector<std::string> s, bool last=false);
							      
  };

#endif
