/*! \file    White.h
 *  \author  Machiel Bos
 *
 * Header file for White.cpp
 *
 * \date 19/ 9/2010  Coimbra
 * \date  7/10/2012  Santa Clara
 */
//=============================================================================
  
  #ifndef __WHITE
    #define __WHITE
    #include <string>
    #include "NoiseModelBaseClass.h"

    class White : public NoiseModelBaseClass
    {
      private:
        std::string   unit;
        int           Nparam;

      public:
        White(void);
        ~White(void) {};
        void    get_covariance(double *param, int m, double *gamma_x);
        void    show(double *param, double *error, double sigma);
        int     get_Nparam(void) {return Nparam;};
        double  get_d_fixed(void) {return 0.0;};
        double  compute_penalty(double *param);
        void    set_noise_parameters(double *params_fixed);
        double  compute_G(double lambda);
        void    compute_impulse_response(int m, double* h);
    };

  #endif
