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

//===================================================================
// Subroutine
//===================================================================

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
//-----------------------------------------
  long julday(int year, int month, int day)
//-----------------------------------------
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


//===================================================================
// Main program
//===================================================================

//--------------------------------
  int main(int argc, char* argv[])
/*--------------------------------
 * Accept line arguments
 */
  {
    int       year,month,day;
    double    MJD,hour,minute,second;
    long int  J;

    if (argc!=7) {
      printf ("Correct input: date2MJD year month day hour minute second\n");
    } else {
      year  = atoi(argv[1]);
      month = atoi(argv[2]);
      day   = atoi(argv[3]);
      J     = julday(year,month,day);
      MJD   = double(J-2400001);  // count from midnight not 12 o'clock
                                  // the 0.5 is added late with the hours
      hour  = atof(argv[4]);
      minute= atof(argv[5]);
      second= atof(argv[6]);
      MJD += hour/24.0 + minute/1440.0 + second/86400.0;
      printf ("year   : %4d\n",year);
      printf ("month  : %4d\n",month);
      printf ("day    : %4d\n",day);
      printf ("hour   : %4.0lf\n",hour);
      printf ("minute : %4.0lf\n",minute);
      printf ("second : %lf\n",second);
      printf ("Julian : %ld\n",J);
      printf ("MJD    : %lf\n",MJD);
    }
    return EXIT_SUCCESS;
  }
