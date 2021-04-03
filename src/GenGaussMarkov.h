/*! \file    GenDaussMarkov.cpp
 *  \author  Machiel Bos
 *
 * Header file for GenGaussMarkov.cpp
 *
 * \date 5/6/2012  Santa Clara
 */
//==============================================================================

  #ifndef __GENGAUSSMARKOV
    #define __GENGAUSSMARKOV
    #include "NoiseModelBaseClass.h"
    #include <complex>
    #include <cstdlib>
    #include "Control.h"

    class GenGaussMarkov : public NoiseModelBaseClass
    {
      private:
        int     Nparam;
        double  d_fixed,phi_fixed;
        double  hyperg_2F1(double a, double b, double c, double z);
        double  backward(double a, double b, double c, double z,
                                                        double F, double Fp1);

      public:
        GenGaussMarkov(double d_fixed_ = sqrt(-1.0));
        void    get_covariance(double *param, int m, double *gamma_x);
        void    show(double *param, double *error);
        int     get_Nparam(void) {return Nparam;};
        double  get_d_fixed(void) {return d_fixed;};
        double  compute_penalty(double *param);
        void    setup_PSD(void);
        double  compute_G(double lambda);
        void    compute_impulse_response(int m, double* h);
    };

  #endif
