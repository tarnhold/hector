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
        int                 n,m,Ngaps,Nparam;
        double              *x,*H,*F,*theta,*C_theta,sigma_eta;

      public:
        Likelihood(void);
        virtual ~Likelihood(void);
        static Likelihood*  getInstance(void);
        virtual void        compute_LeastSquares(double *param, 
			       double& lndeterminant_, double& sigma_eta_)=0;
        double              compute(double *param);
        void                show_leastsquares(void);
        void                show_matrix(const char name[], double *A,
                                                              int m, int n);
        void                set_sigma_eta(double sigma_eta_) {
						   sigma_eta = sigma_eta_;};
        double              get_sigma_eta(void) {return sigma_eta;};
    };

  #endif
