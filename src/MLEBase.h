/*! \file    MLEBase.h
 *  \author  Machiel Bos
 *
 * Header file for MLEBase.cpp. 
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

  #ifndef __MLEBASE
    #define __MLEBASE
    #include "Control.h"
    #include <vector>
    #include "DesignMatrix.h"
    #include "Observations.h"
    #include "NoiseModel.h"

    class MLEBase
    {
      private:
        const double        NaN;

      protected:
        const double        tpi;
        bool                Kashyap;
        int                 n,m,Ngaps,Nparam;
        double              *t,*x,*H,*F,*theta,*C_theta,*C_thetaInv,sigma_eta;
        double              ln_det_C,ln_det_I,AIC,BIC,BIC_c,ln_L;
        double              beta_size,beta_spacing,extra_penalty,BIC_tp;

      public:
        MLEBase(void);
        virtual ~MLEBase(void);
        virtual void  compute_LeastSquares(double *param)=0;
        virtual void  compute_BIC_cs(double *BIC_c)=0;
        void          show_matrix(const char name[], double *A, int m, int n);
        double        compute(double *param);
        void          show_leastsquares(void);
        void          compute_L_and_ICs(double *param);
        double        get_sigma_eta(void) {return sigma_eta;};
        double        get_ln_det_I(void) {return ln_det_I;};
        double        get_ln_L(void) {return ln_L;};
        double        get_AIC(void) {return AIC;};
        double        get_BIC(void) {return BIC;};
        double        get_BIC_tp(void) {return BIC_tp;};
        double        get_BIC_c(void) {return BIC_c;};
    };

  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
