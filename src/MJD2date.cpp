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

//===================================================================
// Subroutines
//===================================================================

//----------------------------------------------------------
  void caldat(long int jul, int& year, int& month, int& day)
/*----------------------------------------------------------
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


//===================================================================
// Main program
//===================================================================

  int main(int argc, char* argv[])
  {

    int       year,month,day,hour,minute;
    double    f,MJD,second;
    long int  J;

    if (argc!=2) {
      printf ("Correct input: MJD2date MJD\n");
    } else {
      MJD    = atof(argv[1]);
      J      = (long int)(MJD + 2400001.0);
      caldat(J,year,month,day);
      f      = 24.0*(MJD-floor(MJD));
      hour   = int(f);
      f      = 60.0*(f - hour);
      minute = int(f);
      second = 60.0*(f - minute);
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
