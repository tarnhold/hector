/*! \file    date2MJD.cpp
 *  \author  Machiel Bos
 *
 * Simple year/month/day to MJD converter
 *
 * \date 5/12/2002  Machiel Bos
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
      printf ("second : %.5lf\n",second);
      printf ("MJD    : %.10lf\n",MJD);
    }
    return EXIT_SUCCESS;
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
