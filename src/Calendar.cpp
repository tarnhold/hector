/*! \file    Calendar.cpp
 *  \author  Machiel Bos
 *
 * Since conversion between MJD and date is needed in DesignMatrix.cpp, I 
 * decided to create a new class with these subroutines. 
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
  #include <cmath>
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include "Calendar.h"

//==============================================================================
// Subroutines
//==============================================================================

namespace
{
  /*
   * Round double to decimal place given by precision.
   */
  double round_double(double value, double precision)
  {
    return floor((value * pow(10, precision) + 0.5)) / pow(10, precision);
  }
}

//--------------------------------------------------------------------
  void Calendar::caldat(long int jul, int& year, int& month, int& day)
/*--------------------------------------------------------------------
 * Compute calendar date for given julian date.
 * 
 * Example: YEAR = 1970, MONTH = 1, DAY = 1, JD = 2440588. 
 * Reference: Fliegel, H. F. and van Flandern, T. C., Communications of 
 * the ACM, Vol. 11, No. 10 (October, 1968).
 *
 * http://www.usno.navy.mil/USNO/astronomical-applications/
 */
  {
    long int   l,n,i,j,k; 

    using namespace std;
    l = jul+68569;
    n = 4*l/146097;
    l = l-(146097*n+3)/4;
    i = 4000*(l+1)/1461001;
    l = l-1461*i/4+31;
    j = 80*l/2447;
    k = l-2447*j/80;
    l = j/11;
    j = j+2-12*l;
    i = 100*(n-49)+i+l;

    year  = i;
    month = j;
    day   = k;

    if (year<1801 || year>2099) {
      cerr << "year " << year << " is out of possible range!" << endl;
      exit(EXIT_FAILURE);
    }
  }



/*!Compute for given year, month and day the Julian day
 * 
 * Example: YEAR = 1970, MONTH = 1, DAY = 1, JD = 2440588. 
 * Reference: Fliegel, H. F. and van Flandern, T. C., Communications of 
 * the ACM, Vol. 11, No. 10 (October, 1968).
 *
 * http://www.usno.navy.mil/USNO/astronomical-applications/
 *
 * \param[in]  year       year
 * \param[in]  month      month
 * \param[in]  day        day
 * \returns{Julian day}
 */
//---------------------------------------------------
  long Calendar::julday(int year, int month, int day)
//---------------------------------------------------
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



/* For given MJD, compute date
 */
//----------------------------------------------------------------------------
  void Calendar::compute_date(double MJD, int& year, int& month, int& day,
                                       int& hour, int& minute, double& second)
//----------------------------------------------------------------------------
  {
    using namespace std;
    long int     J;
    double       f,ds;

    J      = (long int)(MJD + 2400001.0);
    caldat(J,year,month,day);

    // round second to 5th decimal place
    ds = 24.0*3600.0*(MJD-floor(MJD));
    f = round_double(ds, 5);

    hour   = int(f / 3600.0);
    f     -= hour*3600.0;
    minute = int(f / 60.0);
    f     -= minute*60.0;
    second = f;
  }



/* Compute Modified Julian data
 */
//---------------------------------------------------------------------------
  double Calendar::compute_MJD(int year, int month, int day, int hour,
                                                  int minute, double second)
//---------------------------------------------------------------------------
  {
    long int  J;
    double    MJD;

    J    = julday(year,month,day);
    MJD  = double(J-2400001);  // count from midnight not 12 o'clock
                               // the 0.5 is added late with the hours
    MJD += hour/24.0 + minute/1440.0 + second/86400.0;

    return MJD;
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
