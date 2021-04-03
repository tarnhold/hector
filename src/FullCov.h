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
    #include "Likelihood.h"
    #include "NoiseModel.h"
    #include "Observations.h"
    #include "DesignMatrix.h"
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class FullCov : public Likelihood
    {
      private:
        int     N;
        double  *C,*dummyH,*dummyx,*gamma_x,*A,*y,*r,*dummy;

      public:
        FullCov(void);
        ~FullCov(void);
        virtual void    compute_LeastSquares(double *param,
                                   double& lndeterminant, double& sigma_eta);
    };
  
  #endif
  
