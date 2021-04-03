/*! \file   Powerlaw.h
 *  \author Machiel Bos
 *
 * Header file for Powerlaw.cpp
 *
 * \date 17/1/2012   CIIMAR, Porto
 */
//==============================================================================

  #ifndef __POWERLAW
    #define __POWERLAW
    #include "Control.h"
    #include "Observations.h"
    #include "NoiseModelBaseClass.h"
    #include <cmath>
    #include <string>

    class Powerlaw : public NoiseModelBaseClass
    {
      private:
        const double   pi;
        bool           estimate_spectral_index;
        double         param_PSD[2],d_fixed;
        std::string    unit;

      public:
        Powerlaw(double d_fixed_ = sqrt(-1.0));
        void      get_covariance(double *param, int m, double *gamma_x);
        void      show(double *param, double *error, double sigma_eta);
        int       get_Nparam(void);
        double    get_d_fixed(void) {return d_fixed;};
        double    compute_penalty(double *param);
        void      set_noise_parameters(double *params_fixed);
        double    compute_G(double lambda);
        void      compute_impulse_response(int m, double* h);
    };

  #endif
 
