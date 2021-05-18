/*! \file    Likelihood.h
 *  \author  Machiel Bos
 *
 * Header file for Likelihood.cpp. I create here a singleton which is the
 * interface to AmmarGrag/FullCov.
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

  #ifndef __LIKELIHOOD
    #define __LIKELIHOOD
    #include "MLEBase.h"
    #include "AmmarGrag.h"
    #include "FullCov.h"

    class Likelihood 
    {
      private:
        std::string  methodname;
        MLEBase*     method;
        Likelihood(void);
        ~Likelihood(void);

      public:
        //--- Meyers singleton
        static Likelihood& getInstance(void) {
          static Likelihood theLikelihood;
          return theLikelihood;
        }

        void    reset_method(void);
        void    compute_LeastSquares(double *param);
        double  compute(double *param);
        void    show_leastsquares(void); 
        void    compute_BIC_cs(double *BIC_c) {method->compute_BIC_cs(BIC_c);};
        void    compute_L_and_ICs(double *param);
        double  get_sigma_eta(void) {return method->get_sigma_eta();};
        double  get_ln_det_I(void) {return method->get_ln_det_I();};
        double  get_ln_L(void) {return method->get_ln_L();};
        double  get_AIC(void) {return method->get_AIC();};
        double  get_BIC(void) {return method->get_BIC();};
        double  get_BIC_tp(void) {return method->get_BIC_tp();};
        double  get_BIC_c(void) {return method->get_BIC_c();};
    };

  #endif 
