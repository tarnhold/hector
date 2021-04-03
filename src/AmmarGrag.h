/*! \file    AmmarGrag.h
 *  \author  Machiel Bos
 *
 * Header file for AmmarGrag.cpp
 *
 * \date 9/2/2012  CIIMAR, Porto
 */
//==============================================================================

  #ifndef __AMMARGRAG
    #define __AMMARGRAG
    #if OMP == 1
      #include <omp.h>
    #endif
    #include <complex>
    #include <fftw3.h>
    #include "Likelihood.h"
    #include "NoiseModel.h"

    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class AmmarGrag : public Likelihood
    {
      protected:
        int                     max_order,Nthreads,ny;
        double                  *A1,*A2,*y1,*y2,*G1,*G2,*l1,*l2,*gamma_x;
        double                  *dummy,*Qy,*QA,*Qt,*M;
        fftw_complex            *F_H,*F_x,*F_F,*F_l1,*F_l2,*F_dummy;
        fftw_plan               plan_backward,plan_forward;

        double    step1(double *gamma_x, double **l1, double **l2);
        void      step2(int n_columns, fftw_complex *F_B, 
						     double *B1, double *B2);

      public:
        AmmarGrag(void);
        ~AmmarGrag(void);
        virtual void    compute_LeastSquares(double *param,
                                  double& lndeterminant, double& sigma_eta);
    };

  #endif 
