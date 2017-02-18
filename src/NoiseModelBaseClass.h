/*! \file    NoiseModelBaseClass.h
 *  \author  Machiel Bos
 *
 * Template of class to which the stationary noise models must conform.
 *
 * \date  7/10/2012  Santa Clara
 * \date 15/ 7/2015  Santa Clara
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

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
