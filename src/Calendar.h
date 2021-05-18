/*! \file    Calendar.h
 *  \author  Machiel Bos
 *
 * Header file for Calendar.cpp
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

  #ifndef __CALENDAR
    #define __CALENDAR
 
    class Calendar
    {
      private:
        void      caldat(long int jul, int& year, int& month, int& day);

      public:
        long int  julday(int year, int month, int day);
        void      compute_date(double MJD, int& year, int& month, int& day,
				       int& hour, int& minute, double& second);
        double    compute_MJD(int year, int month, int day, int hour, 
						    int minute, double second);
        void      MJD_to_ISO8601(double MJD, std::string& date);
    };

  #endif
