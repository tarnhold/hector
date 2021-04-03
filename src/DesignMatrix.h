/*! \file   DesignMatrix.h
 *  \author Machiel Bos
 *
 * Header file for DesignMatrix.cpp
 *
 * \date 13/1/2012   Coimbra library
 */
//=============================================================================

  #ifndef __DESIGNMATRIX
    #define __DESIGNMATRIX
    #include "Control.h"
    #include "Observations.h"

    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class DesignMatrix
    {
      private:
        static bool            instanceFlag;
        static DesignMatrix    *singleton;
        const double           tpi;
        std::string            unit;
        int                    n,m,Ngaps,n_periodic_signals,n_offsets;
        int                    n_postseismiclog,n_postseismicexp;
        double                 dt,*H,*F,*periods,th;
        std::vector<double>    offsets;
        std::vector<LogEntry>  postseismiclog;
        std::vector<ExpEntry>  postseismicexp;
        bool                   quadratic_term;
        bool                   seasonal_signal,halfseasonal_signal;

      public:
        DesignMatrix(void);
        ~DesignMatrix(void);
        static DesignMatrix* getInstance(void);
        static void          destroyInstance(void);
        void    get_H(int& n_, double **H_); 
        void    get_F(double **F_); 
        void    show_results(double *theta, double *error);
        void    compute_xhat(const double *theta);
    };

  #endif
