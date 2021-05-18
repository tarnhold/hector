/*! \file    ModelSpectrum.h
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

  #ifndef __MODELSPECTRUM
    #define __MODELSPECTRUM
    #include <complex>
    #include "NoiseModel.h"
    #include "Spectrum.h"

    class ModelSpectrum
    {
      private:
        int            n_simulations,m,ms,segments;
        bool           montecarlo;
        const double   pi,tpi;
        double         fs,sigma,scale;
        void           compute_confidence_levels(void);

      public:
        ModelSpectrum(void);
        void           compute(void);
    };

  #endif
