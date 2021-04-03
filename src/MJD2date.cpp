/*! \file    MJD2date.cpp
 *  \author  Machiel Bos
 *  \version 1.0
 *
 * Simple MJD to date converter
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

  int main(int argc, char* argv[])
  {

    int       year,month,day,hour,minute;
    double    MJD,second;
    Calendar  calendar;

    if (argc!=2) {
      printf ("Correct input: MJD2date MJD\n");
    } else {
      MJD    = atof(argv[1]);
      calendar.compute_date(MJD,year,month,day,hour,minute,second);
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
