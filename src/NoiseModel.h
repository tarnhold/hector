/*! \file    NoiseModel.h
 *  \author  Machiel Bos
 *
 * Header file for NoiseModel.cpp
 *
 * \date 13/1/2011  Santa Clara
 */
//=============================================================================

  #ifndef __NOISEMODEL
    #define __NOISEMODEL
    #include "Control.h"
    #include "NoiseModelBaseClass.h"

    class NoiseModel 
    {
      private:
        static bool          instanceFlag;
        static NoiseModel    *singleton;
        int                  Nparam,Nmodels,*NparamIndv;
        double               *phi,*fraction_PSD;
        char                 **noisemodel;
        NoiseModelBaseClass  **modelIndv;

      public:
        static NoiseModel*   getInstance(void);

        NoiseModel(void);
        ~NoiseModel(void);
        void    get_covariance(double *param,int m,double *gamma_x);
        void    show(double *param, double *error);
        int     get_Nparam(void) {return Nparam;};
        double  compute_penalty(double *param);
        void    setup_PSD(void);
        double  compute_G(double lambda);

    };
 
  #endif 
