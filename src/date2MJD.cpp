/*! \file    date2MJD.cpp
 *  \author  Machiel Bos
 *
 * Simple year/month/day to MJD converter
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
//===================================================================
  #include <cstdio>
  #include <cstdlib>
  #include <cmath>
  #include <iostream>
  #include <ostream>
  #include "Calendar.h"

//===================================================================
// Main program
//===================================================================

//--------------------------------
  int main(int argc, char* argv[])
/*--------------------------------
 * Accept line arguments
 */
  {
    int       year,month,day,hour,minute;
    double    MJD,second;
    Calendar  calendar;

    if (argc!=7) {
      printf ("Correct input: date2MJD year month day hour minute second\n");
    } else {
      year  = atoi(argv[1]);
      month = atoi(argv[2]);
      day   = atoi(argv[3]);
      hour  = atoi(argv[4]);
      minute= atoi(argv[5]);
      second= atof(argv[6]);
      MJD = calendar.compute_MJD(year,month,day,hour,minute,second);
      printf ("year   : %4d\n",year);
      printf ("month  : %4d\n",month);
      printf ("day    : %4d\n",day);
      printf ("hour   : %4d\n",hour);
      printf ("minute : %4d\n",minute);
      printf ("second : %lf\n",second);
      printf ("MJD    : %lf\n",MJD);
    }
    return EXIT_SUCCESS;
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
