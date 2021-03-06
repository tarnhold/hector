/*! \filename RemoveOutliers.h
 *  \author   Machiel Bos
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

  #ifndef REMOVEOUTLIERS__
    #define REMOVEOUTLIERS__
    #include "Observations.h"
    #include "DesignMatrix.h"
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class RemoveOutliers
    {
      private:
        double   factor;
       
      public:
        RemoveOutliers(void);
        void     compute_LeastSquares(int& m, double **r, double **t);
        int      compare (const void * a, const void * b);
        void     data_snooping(void);
    };

  #endif 

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
