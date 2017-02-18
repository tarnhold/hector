/*! \file    MLEBase.h
 *  \author  Machiel Bos
 *
 * Header file for MLEBase.cpp. 
 *
 * For offset detection I need more flexibility. Therefore I have split
 * compute_LeastSquares into prepare_covariance and compute_LeastSquares.
 * The first subroutine performs Cholesky decomposition or computes l1 and l2
 * whichever are the preperatory tasks needed to speed up Least-Squares.
 * If matrix C remains constant then I only have to do this once and afterwards
 * I can make changes in the design matrix H and quickly recompute the 
 * likelihood.
 *
 * \date 13/1/2012  Coimbra library
 * \date  7/6/2016  Santa Clara
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
        int                 n,m,Ngaps,Nparam;
        double              *t,*x,*H,*F,*theta,*C_theta,sigma_eta;
        double              ln_det_C,AIC,BIC,ln_L;

      public:
        MLEBase(void);
        virtual ~MLEBase(void);
        virtual void  prepare_covariance(double *param)=0;
        virtual void  compute_LeastSquares(double *param)=0;
        void          show_matrix(const char name[], double *A, int m, int n);
        double        compute(double *param, bool quick=false);
        void          show_leastsquares(void);
        void          compute_L_and_ICs(double *param);
        double        get_sigma_eta(void) {return sigma_eta;};
        double        get_ln_L(void) {return ln_L;};
        double        get_AIC(void) {return AIC;};
        double        get_BIC(void) {return BIC;};
    };

  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
