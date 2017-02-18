/*! \file   FullCov.h
 *  \author Machiel Bos
 *
 * Header file for the class that uses the full covariance matrix and
 * afterwards eliminates the rows and gaps for missing data.
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
//=============================================================================

  #ifndef __FULLCOV
    #define __FULLCOV
    #include "MLEBase.h"
    #include "NoiseModel.h"
    #include "Observations.h"
    #include "DesignMatrix.h"
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class FullCov : public MLEBase
    {
      private:
        int     N;
        double  *C,*dummyt,*dummyH,*dummyx,*gamma_x,*A,*y,*r,*dummy;

      public:
        FullCov(void);
        ~FullCov(void);
        virtual void    compute_LeastSquares(double *param);
        virtual void    compute_BIC_cs(double *BIC_c);
    };
  
  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
