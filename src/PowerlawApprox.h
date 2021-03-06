/*! \file   PowerlawApprox.h
 *  \author Machiel Bos
 *
 * Header file for PowerlawApprox.cpp
 *
 *  This script is part of Hector 1.7.2
 *
 *  Hector is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  Hector is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Hector. If not, see <http://www.gnu.org/licenses/>
 *
 */
//==============================================================================

  #ifndef __POWERLAWAPPROX
    #define __POWERLAWAPPROX
    #include "Control.h"
    #include "NoiseModelBaseClass.h"
    #include "Observations.h"
    #include <math.h>
    #include <string>
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
        std::string    unit;
        int            ms;

      public:
        PowerlawApprox(double d_fixed_ = sqrt(-1.0));
        void      get_covariance(double *param, int m, double *gamma_x);
        void      show(double *param, double *error, double sigma);
        int       get_Nparam(void);
        double    get_d_fixed(void) {return d_fixed;};
        double    compute_penalty(double *param);
        void      set_noise_parameters(double *params_fixed);
        double    compute_G(double lambda);
        void      compute_impulse_response(int m, double* h);
    };

  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
