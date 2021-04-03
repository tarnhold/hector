/*! \file   LevinsonGap.h
 *  \author Machiel Bos
 *
 * Header file for the class that uses the Levinson method on data that
 * has data gaps. The base class contains the mactrick routine plus
 * the necessary variables.
 * 
 * \date 3/2/2012  CIIMAR, Porto
 */
//=============================================================================

  #ifndef __LEVINSONGAP
    #define __LEVINSONGAP
    #include "Likelihood.h"
    #include "NoiseModel.h"

    //--- ATLAS is written in C
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    }

    class LevinsonGap : public Likelihood
    {
      private:
        int     Nparam;
        double  *U,*V,*dummyH,*dummyx,*gamma_x,*A,*y,*r,*dummy,*G,*dummyF;
        double  *M,*t,*dummy1,*dummy2,*dummy3,*QA,*Qy,*Qt;
        void    mactrick(double *gamma_x, double& lndeterminant_,
                                             double *A, double *y, double *G);
      public:
        LevinsonGap(void);
        ~LevinsonGap(void);
        virtual void    compute_LeastSquares(double *param,
                              double& lndeterminant_, double& sigma_eta_);
    };
  
  #endif
  
