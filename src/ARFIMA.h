/*! \file   ARFIMA.h
 *  \author Machiel Bos
 *
 * Header file for ARFIMA.cpp
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
 *  along with Hector.  If not, see <http://www.gnu.org/licenses/>
 *
 */
//=============================================================================
  
  #ifndef __ARFIMA
    #define __ARFIMA
    #include <complex>
    #include <string>
    #include "Control.h"
    #include "NoiseModelBaseClass.h"
    #include "HyperGeo.h"

    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class ARFIMA : public NoiseModelBaseClass
    {
      private:
        bool                 estimate_spectral_index;
        int                  p,q,Nparam;
        double               d_PSD,*psi_PSD,d_fixed;
        std::complex<double> *rho_PSD;
        std::string          unit;
        void                 compute_psi(double *MA, double *psi);
        void                 ZindeWalsh(double *AR, double *MA,
						    int m, double *gamma_x);
        void                 DoornikOoms(double *AR, double d, double *MA,
						    int m, double *gamma_x);

      public:
        ARFIMA(double d_fixed_ = sqrt(-1.0));
        ~ARFIMA(void);
        void    find_roots(double *AR, std::complex<double> *rho);
        void    find_coefficients(std::complex<double> *rho, double *AR);
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
