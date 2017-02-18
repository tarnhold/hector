/*! \file    Likelihood.h
 *  \author  Machiel Bos
 *
 * Header file for Likelihood.cpp. I create here a singleton which is the
 * interface to AmmarGrag/FullCov.
 *
 * \date 9/6/2016  Santa Clara
 */
//==============================================================================

  #ifndef __LIKELIHOOD
    #define __LIKELIHOOD
    #include "MLEBase.h"
    #include "AmmarGrag.h"
    #include "FullCov.h"

    class Likelihood 
    {
      private:
        std::string  methodname;
        MLEBase*     method;
        Likelihood(void);
        ~Likelihood(void);

      public:
        //--- Meyers singleton
        static Likelihood& getInstance(void) {
          static Likelihood theLikelihood;
          return theLikelihood;
        }

        void    reset_method(void);
        void    prepare_covariance(double *param);
        void    compute_LeastSquares(double *param);
        double  compute(double *param, bool quick=false);
        void    show_leastsquares(void);
        void    compute_L_and_ICs(double *param);
        double  get_sigma_eta(void) {return method->get_sigma_eta();};
        double  get_ln_L(void) {return method->get_ln_L();};
        double  get_AIC(void) {return method->get_AIC();};
        double  get_BIC(void) {return method->get_BIC();};
    };

  #endif 

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
