/*! \file    NoiseModel.h
 *  \author  Machiel Bos
 *
 * Header file for NoiseModel.cpp
 *
 * \date 13/1/2011  Santa Clara
 */
//=============================================================================

  #ifndef __NOISEMODEL
    #define __NOISEMODEL
    #include "Control.h"
    #include "NoiseModelBaseClass.h"
    #include <fftw3.h>
    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };


    class NoiseModel 
    {
      private:
        static bool          instanceFlag;
        static NoiseModel    *singleton;
        int                  Nparam,Nmodels,*NparamIndv;
        double               *phi,*fraction_fixed,*h,*w,sigma_fixed;
        char                 **noisemodel;
        NoiseModelBaseClass  **modelIndv;
        fftw_complex         *F_h,*F_w;
        fftw_plan            plan_backward,plan_forward;
        const gsl_rng_type   *T_random;
        gsl_rng              *r_random;

      public:
        static NoiseModel*   getInstance(void);

        NoiseModel(void);
        ~NoiseModel(void);
        void    get_covariance(double *param,int m,double *gamma_x);
        void    show(double *param, double *error);
        int     get_Nparam(void) {return Nparam;};
        double  compute_penalty(double *param);
        void    setup_PSD(void);
        double  compute_G(double lambda);
        void    setup_MonteCarlo(int m);
        void    create_noise(int m, double *w);
    };
 
  #endif 
