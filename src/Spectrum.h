/*! \file    Spectrum.h
 *  \author  Machiel Bos
 *
 * Header file for ModelSpectrum.cpp.
 *
 *  This script is part of Hector 1.8
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

  #ifndef __SPECTRUM
    #define __SPECTRUM
    #include <complex>
    #include <fftw3.h>

    class Spectrum
    {
      private:
        int             n,segments,N,L,K,m;
        const double    pi;
        fftw_plan       plan_forward;
        fftw_complex    *Y;
        double          *dummy,fraction,dt,*f,fs,scale;
        double          (Spectrum::*windowfunction)(double);

      public:
        Spectrum(int n_, int segments_, double fs_);
        ~Spectrum(void);
        double    Parzen(double s);
        double    Hann(double s);
        void      remove_mean(double *x, int n);
        double    welch(double *y, double *G);
        void      get_freq(double *f_);
    };

  #endif
