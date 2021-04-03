/*! \file    Calendar.h
 *  \author  Machiel Bos
 *  \version 1.1
 *
 * Header file for Calendar.cpp
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
    };

  #endif
