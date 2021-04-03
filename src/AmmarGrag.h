/*! \file    AmmarGrag.h
 *  \author  Machiel Bos
 *
 * Header file for AmmarGrag.cpp
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
 *  along with Hector.  If not, see <http://www.gnu.org/licenses/>
 *
 */
//==============================================================================

  #ifndef __AMMARGRAG
    #define __AMMARGRAG
    #include <omp.h>
    #include <complex>
    #include <fftw3.h>
    #include "MLEBase.h"
    #include "NoiseModel.h"

    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class AmmarGrag : public MLEBase
    {
      private:
        const double            tpi,LARGE,EPS;
        int                     max_order,Nthreads,ny,*index;
        double                  *A1,*A2,*y1,*y2,*G1,*G2,*l1,*l2,*gamma_x;
        double                  *dummy,*Qy,*QA,*Qt,*M,*M2,*t1,*t2;
        fftw_complex            *F_H,*F_x,*F_l1,*F_l2,*F_dummy;
        fftw_plan               plan_backward,plan_forward;

        double    step1(double *gamma_x, double **l1, double **l2);
        void      step2(int n_columns, fftw_complex *F_B, 
						     double *B1, double *B2);

      public:
        AmmarGrag(void);
        ~AmmarGrag(void);
        virtual void    compute_LeastSquares(double *param);
        virtual void    compute_BIC_cs(double *BIC_c);
    };

  #endif 
