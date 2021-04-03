/*! \file    EstimateTrend.cpp
 *  \author  Machiel Bos
 *  \version 1.0
 *
 * The main program that simply calls the minimizer class. 
 *
 * \date 16/1/2012  Santa Clara
 */
//=============================================================================
  #include "Minimizer.h"
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <ctime>

//=============================================================================
// Main program
//=============================================================================

  int main(void)
  {
    using namespace std;
    cout << endl 
         << "***************************" << endl 
         << "    Hector, version 1.0    " << endl    
         << "***************************" << endl;

    Minimizer      minimizer;
    time_t         start,end;

    time(&start);
    minimizer.solve();
    time(&end);
    cout << "Total computing time: " << difftime(end,start) << " sec" << endl;

    return EXIT_SUCCESS;
  }
