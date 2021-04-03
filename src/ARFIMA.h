/*! \file   ARFIMA.h
 *  \author Machiel Bos
 *
 * Header file for ARFIMA.cpp
 *
 *  \date 19/9/2010  Coimbra
 */
//=============================================================================
  
  #ifndef __ARFIMA
    #define __ARFIMA
    #include <complex>
    #include "Control.h"
    #include "NoiseModelBaseClass.h"
    #include "HyperGeo.h"

    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class ARFIMA : public NoiseModelBaseClass
    {
      private:
        bool                 estimate_spectral_index;
        int                  p,q,Nparam;
        double               d_PSD,*psi_PSD,d_fixed;
        std::complex<double> *rho_PSD;
        void                 compute_psi(double *MA, double *psi);
        void                 ZindeWalsh(double *AR, double *MA,
						    int m, double *gamma_x);
        void                 DoornikOoms(double *AR, double d, double *MA,
						    int m, double *gamma_x);

      public:
        ARFIMA(double d_fixed_ = sqrt(-1.0));
        ~ARFIMA(void);
        void    find_roots(double *AR, std::complex<double> *rho);
        void    find_coefficients(std::complex<double> *rho, double *AR);
        void    get_covariance(double *param, int m, double *gamma_x);
        void    show(double *param, double *error);
        int     get_Nparam(void) {return Nparam;};
        double  compute_penalty(double *param);
        void    setup_PSD(void);
        double  compute_G(double lambda);
        void    compute_impulse_response(int m, double* h);
    };

  #endif
