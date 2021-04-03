/*! \file   LevinsonNoGap.h
 *  \author Machiel Bos
 *
 * Header file for the class that uses the Levinson method on data that
 * has no data gaps. The base class contains the mactrick routine plus
 * the necessary variables.
 * 
 * \date 18/1/2012  Santa Clara.
 */
//=============================================================================

  #ifndef __LEVINSONNOGAP
    #define __LEVINSONNOGAP
    #include "Likelihood.h"
    #include "NoiseModel.h"

    //--- ATLAS is written in C
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    }

    class LevinsonNoGap : public Likelihood
    {
      private:
        int     Nparam;
        double  *U,*V,*dummyH,*dummyx,*gamma_x,*A,*y,*r,*dummy;
        void    mactrick(double *gamma_x, double& lndeterminant_,
                                             	double *A, double *y);

      public:
        LevinsonNoGap(void);
        ~LevinsonNoGap(void);
        virtual void    compute_LeastSquares(double *param,
                                double& lndeterminant_, double& sigma_eta_);
    };
  
  #endif
  
