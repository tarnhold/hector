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
    #include "NoiseModelBaseClass.h"

    class White : public NoiseModelBaseClass
    {
      private:
        int  Nparam;

      public:
        White(void);
        ~White(void) {};
        void    get_covariance(double *param, int m, double *gamma_x);
        void    show(double *param, double *error);
        int     get_Nparam(void) {return Nparam;};
        double  compute_penalty(double *param);
        void    setup_PSD(void);
        double  compute_G(double lambda);
    };

  #endif
