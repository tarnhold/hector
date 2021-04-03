/*! \file    ModelSpectrum.h
 *  \author  Machiel Bos
 *
 * Header file for ModelSpectrum.cpp. This class provides a subroutine
 * for creating the spectrum of an ARFIMA noise model.
 *
 * \date 22/7/2011  CIIMAR,  Porto
 * \date 25/7/2011  Coimbra
 */
//==============================================================================

  #ifndef __MODELSPECTRUM
    #define __MODELSPECTRUM
    #include <complex>
    #include "NoiseModel.h"

    class ModelSpectrum
    {
      private:
        const double   pi,tpi;
        double         fs,sigma,scale;

      public:
        ModelSpectrum(void);
        void           compute(void);
    };

  #endif
