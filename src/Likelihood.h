/*! \file    Likelihood.h
 *  \author  Machiel Bos
 *
 * Header file for Likelihood.cpp. Since the Minimizer class calls this
 * class repeatedly, it makes sense to turn this class into a singleton
 * to store the design matrix and observations.
 *
 * \date 13/1/2012  Coimbra library
 */
//=============================================================================

  #ifndef __LIKELIHOOD
    #define __LIKELIHOOD
    #include "Control.h"
    #include <fftw3.h>
    #include "DesignMatrix.h"
    #include "Observations.h"
    #include "NoiseModel.h"

    class Likelihood
    {
      private:
        static bool         instanceFlag;
        static Likelihood   *singleton;

      protected:
        const double        tpi;
        int                 m,n,Nparam,Ngaps;
        double              *theta,*H,*x,*t,*C_theta,sigma_eta,lndeterminant;
        double              *F;

      public:
        Likelihood(void);
        virtual ~Likelihood(void)=0;
        static Likelihood*  getInstance(void);
        void                show_leastsquares(void);
        int                 get_Nparam(void) {return Nparam;};
        virtual void        compute_LeastSquares(double *param, 
			       double& lndeterminant_, double& sigma_eta_)=0;
        double              compute(double *param);
        void                show_matrix(const char name[], double *A,
							      int m, int n);
        void                show_Fmatrix(const char name[], fftw_complex *A,
							      int m, int n);
    };

  #endif
