/*! \file   FullCov.h
 *  \author Machiel Bos
 *
 * Header file for the class that uses the full covariance matrix and
 * afterwards eliminates the rows and gaps for missing data.
 * 
 * \date  7/10/2012  Tomar
 * \date 12/ 7/2015  Santa Clara
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
        virtual void    prepare_covariance(double *param);
        virtual void    compute_LeastSquares(double *param);
    };
  
  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
