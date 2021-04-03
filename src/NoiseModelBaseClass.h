/*! \file    NoiseModelBaseClass.h
 *  \author  Machiel Bos
 *
 * Template of class to which the stationary noise models must conform.
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

  #ifndef __NOISEMODELBASECLASS
    #define __NOISEMODELBASECLASS

    class NoiseModelBaseClass
    {
      public:
        virtual ~NoiseModelBaseClass(void) {};
        virtual void    get_covariance(double *param,int m,double *gamma_x)=0;
        virtual void    show(double *param, double *error, double sigma)=0;
        virtual int     get_Nparam(void)=0;
        virtual double  get_d_fixed(void)=0;
        virtual double  compute_penalty(double *param)=0;
        virtual void    set_noise_parameters(double *params_fixed)=0;
        virtual double  compute_G(double lambda)=0;
        virtual void    compute_impulse_response(int m, double* h)=0;
    };
 
  #endif 
