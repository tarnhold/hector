/*! \file    VaryingPeriodic.cpp
 *  \author  Machiel Bos
 *
 * Header file for VaryingPeriodic.cpp
 *
 *  This script is part of Hector 1.9
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

  #ifndef __VARYINGPERIODIC
    #define __VARYINGPERIODIC
    #include "NoiseModelBaseClass.h"
    #include <complex>
    #include <cstdlib>
    #include "Control.h"
    #include "cblas.h"

    class VaryingPeriodic : public NoiseModelBaseClass
    {
      private:
        const double pi,NaN;
        std::string  unit;
        int          Nparam;
        double       omega_0,phi_fixed;

         void mactrick(int m, double *gamma_x, double *h);

      public:
        VaryingPeriodic(double T0 = 365.25);
        void    get_covariance(double *param, int m, double *gamma_x);
        void    show(double *param, double *error, double sigma);
        int     get_Nparam(void) {return Nparam;};
        double  get_d_fixed(void) {return 0.0;};
        double  compute_penalty(double *param);
        void    set_noise_parameters(double *params_fixed);
        double  compute_G(double lambda);
        void    compute_impulse_response(int m, double* h);
    };

  #endif
