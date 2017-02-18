/*! \file   Minimizer.h
 *  \author Machiel Bos
 *
 * Header file for Minimizer.cpp
 *
 * \date 13/1/2012  Santa Clara
 */
//=============================================================================

  #ifndef __MINIMIZER
    #define __MINIMIZER
    #include "Control.h"
    #include "Likelihood.h"
    #include "NoiseModel.h"
    #include "DesignMatrix.h"
    #include <gsl/gsl_multimin.h>

    class Minimizer
    {
      private:
        int            Nparam;
        double         *param,*error;

        void   fill_X(int k0, double s0, int k1, double s1, double *X);
        void   show_L_and_ICs(void);
        void   compute_inv_Fisher(double* C);
        void   compute_confidence_intervals(double *param);

      public:
        Minimizer(void);
        ~Minimizer(void);
        int   get_Nparam(void) {return Nparam;};
        void  get_param(double *param_);
        void  solve(void);
    };
 
  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
