/*! \file   PowerlawApprox.h
 *  \author Machiel Bos
 *
 * Header file for PowerlawApprox.cpp
 *
 * \date 3/2/2012   CIIMAR, Porto
 */
//==============================================================================

  #ifndef __POWERLAWAPPROX
    #define __POWERLAWAPPROX
    #include "Control.h"
    #include "NoiseModelBaseClass.h"
    #include "Observations.h"
    #include <math.h>
    #include <fftw3.h>
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };


    class PowerlawApprox : public NoiseModelBaseClass
    {
      private:
        const double   pi;
        bool           estimate_spectral_index;
        double         d_fixed;
        int            ms;

      public:
        PowerlawApprox(double d_fixed_ = sqrt(-1.0));
        void      get_covariance(double *param, int m, double *gamma_x);
        void      show(double *param, double *error);
        int       get_Nparam(void);
        double    get_d_fixed(void) {return d_fixed;};
        double    compute_penalty(double *param);
        void      setup_PSD(void);
        double    compute_G(double lambda);
        void      compute_impulse_response(int m, double* h);
    };

  #endif
