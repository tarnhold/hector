/*! \file    NoiseModel.h
 *  \author  Machiel Bos
 *
 * Header file for NoiseModel.cpp
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
        const double         NaN,tpi;
        static bool          instanceFlag;
        static NoiseModel    *singleton;
        int                  Nparam,Nmodels,*NparamIndv;
        double               *phi,*fraction_fixed,*h,*w,sigma_fixed;
        std::string          noisemodel[15],unit;
        NoiseModelBaseClass  **modelIndv;
        fftw_complex         *F_h,*F_w;
        fftw_plan            plan_backward,plan_forward;
        const gsl_rng_type   *T_random;
        gsl_rng              *r_random;

        NoiseModel(void);
        ~NoiseModel(void);
        NoiseModel(const NoiseModel &) = delete;
        NoiseModel& operator=(const NoiseModel &) = delete;
        double  compute_fraction(int i, double *param);

      public:
        //--- Meyers singleton
        static NoiseModel& getInstance(void) {
          static NoiseModel theNoiseModel;
          return theNoiseModel;
        }

        void    get_covariance(double *param,int m,double *gamma_x);
        void    show(double *param, double *error, double sigma);
        int     get_Nparam(void) {return Nparam;};
        double  compute_penalty(double *param);
        void    set_noise_parameters(double *params_fixed);
        double  compute_G(double lambda);
        void    setup_MonteCarlo(int m);
        void    create_noise(int m, double *w);
    };
 
  #endif 

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
