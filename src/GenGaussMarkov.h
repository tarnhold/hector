/*! \file    GenDaussMarkov.cpp
 *  \author  Machiel Bos
 *
 * Header file for GenGaussMarkov.cpp
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

  #ifndef __GENGAUSSMARKOV
    #define __GENGAUSSMARKOV
    #include "NoiseModelBaseClass.h"
    #include <complex>
    #include <cstdlib>
    #include "Control.h"
    #include "Observations.h"

    class GenGaussMarkov : public NoiseModelBaseClass
    {
      private:
        int          Nparam;
        std::string  unit;
        double       d_fixed,phi_fixed;
        double       hyperg_2F1(double a, double b, double c, double z);
        double       backward(double a, double b, double c, double z,
                                                        double F, double Fp1);

      public:
        explicit GenGaussMarkov(double d_fixed_ = sqrt(-1.0));
        void    get_covariance(double *param, int m, double *gamma_x);
        void    show(double *param, double *error, double sigma);
        int     get_Nparam(void) {return Nparam;};
        double  get_d_fixed(void) {return d_fixed;};
        double  compute_penalty(double *param);
        void    set_noise_parameters(double *params_fixed);
        double  compute_G(double lambda);
        void    compute_impulse_response(int m, double* h);
    };

  #endif
