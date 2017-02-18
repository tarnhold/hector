/*! \file   Minimizer.h
 *  \author Machiel Bos
 *
 * Header file for Minimizer.cpp
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

  #ifndef __MINIMIZER
    #define __MINIMIZER
    #include "Control.h"
    #include "Likelihood.h"
    #include "NoiseModel.h"
    #include "DesignMatrix.h"
    #include <gsl/gsl_multimin.h>
    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>

    class Minimizer
    {
      private:
        int                  Nparam;
        double               *param,*error;
        bool                 randomise_first_guess;
        const gsl_rng_type   *T_random;
        gsl_rng              *r_random;

        void   fill_X(int k0, double s0, int k1, double s1, double *X);
        void   show_L_and_ICs(void);
        void   compute_inv_Fisher(double* C);
        void   compute_confidence_intervals(double *param);

      public:
        Minimizer(void);
        ~Minimizer(void);
        Minimizer(const Minimizer &) = delete;
        Minimizer& operator=(const Minimizer &) = delete;

        int   get_Nparam(void) {return Nparam;};
        void  get_param(double *param_);
        void  solve(void);
    };
 
  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
