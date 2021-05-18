/*! \filename JSON.cpp
 *  \author   Machiel Bos
 *
 * To run Hector in batch mode and/or include it in other scripts requires
 * that the output is parsed. The printout on the screen is not the most
 * convenient way to extract all the results. Therefore, this class provides
 * some conenient subroutines to add values to a JSON file.
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
 * \date 9/4/2020  Santa Clara (Coimbra)
 */
//==============================================================================
  #include "JSON.h"
  #include <iostream>
  #include <ostream>
  #include <fstream>

//==============================================================================
// Subroutines
//==============================================================================

//---!!---------------------------
  JSON::JSON(std::string filename)
//---!!---------------------------
  {
    Control &control = Control::getInstance();

    using namespace std;
    try {
      json_required = control.get_bool("JSON");
    }
    catch (exception e) {
      json_required = false;
    }

    if (json_required==true) { 
      fp.open(filename, ios::out);
      if (fp.is_open()==true) {
        fp << "{" << endl;
        indent = 2;
      } else {
        cerr << "Could not open file " << filename << endl;
        exit(EXIT_FAILURE);
      }
    }
  }



//---!!------------
  JSON::~JSON(void)
//---!!------------
  {
    using namespace std;
    if (json_required==true) {
      fp << "}" << endl;
      fp.close();
    }
  }



//----------------------------------------------------
  void JSON::start_dictionary(const std::string label)
//----------------------------------------------------
  {
    int    j;

    using namespace std;
    if (json_required==true) {
      for (j=0;j<indent;j++) fp << " ";
      fp << "\"" << label << "\" : {" << endl;
      indent += 2; 
    }
  }



//------------------------------------
  void JSON::end_dictionary(bool last)
//------------------------------------
  {
    int  j;

    using namespace std;
    if (json_required==true) {
      indent -= 2; 
      //--- remove last written comma
      for (j=0;j<indent;j++) fp << " ";
      if (last==true)  fp << "}" << endl;
      else             fp << "}," << endl;
    }
  }



//---------------------------------------------------------------
  void JSON::write_int(const std::string label, int i, bool last)
//---------------------------------------------------------------
  {
    int    j;

    using namespace std;
    if (json_required==true) {
      for (j=0;j<indent;j++) fp << " ";
      if (last==true) fp << "\"" << label << "\" : " << i << endl;
      else            fp << "\"" << label << "\" : " << i << "," << endl;
    }
  }



//---------------------------------------------------------------------
  void JSON::write_double(const std::string label, double y, bool last)
//---------------------------------------------------------------------
  {
    int    j;

    using namespace std;
    if (json_required==true) {
      for (j=0;j<indent;j++) fp << " ";
      if (last==true) fp << "\"" << label << "\" : " << y << endl;
      else            fp << "\"" << label << "\" : " << y << "," << endl;
    }
  }



//-------------------------------------------------------------------------
  void JSON::write_double_list(const std::string label, int n, double *y,
								 bool last)
//-------------------------------------------------------------------------
  {
    int    j;

    using namespace std;
    if (json_required==true) {
      for (j=0;j<indent;j++) fp << " ";
      fp << "\"" << label << "\" : [";
      for (j=0;j<n;j++) {
        fp  << y[j];
        if (j<n-1) fp << ",";
      }
      if (last==false)  fp << "]," << endl;
      else              fp << "]" << endl;
    }
  }



//---------------------------------------------------------------------------
  void JSON::write_chararray_list(const std::string label, 
			               std::vector<std::string> s, bool last)
//---------------------------------------------------------------------------
  {
    int    j,n=s.size();

    using namespace std;
    if (json_required==true) {
      for (j=0;j<indent;j++) fp << " ";
      fp << "\"" << label << "\" : [";
      for (j=0;j<n;j++) {
        fp  << "\"" << s[j] << "\"";
        if (j<n-1) fp << ",";
      }
      if (last==false)  fp << "]," << endl;
      else              fp << "]" << endl;
    }
  }


/* In my opinion this subroutine is the core trick of the singleton
 * implementation. A pointer is made static and points to the instance.
 * This is Meyers singleton.
 */
//---------------------------------------------
  JSON& JSON::getInstance(std::string filename)
//---------------------------------------------
  {
    static JSON singleton(filename);

    return singleton;
  }
