/*! \file   DesignMatrix.h
 *  \author Machiel Bos
 *
 * Header file for DesignMatrix.cpp
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

  #ifndef __DESIGNMATRIX
    #define __DESIGNMATRIX
    #include "Control.h"
    #include "Observations.h"
    #include <gsl/gsl_sf.h>

    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class DesignMatrix
    {
      private:
        static DesignMatrix    *singleton;
        const double           tpi,pi;
        std::string            unit;
        bool                   estimate_multitrend,varying_seasonal;
        int                    n,m,Ngaps,n_periodic_signals,n_offsets,n_breaks;
        int                    n_postseismiclog,n_postseismicexp,n_ssetanh;
        int                    degree_polynomial,n_channels,index_offset;
        int                    varyingseasonal_N,index_seasonal;
        double                 dt,*H,*F,*periods,th,t0,*t,median;
        std::vector<double>    offsets,breaks;
        std::vector<LogEntry>  postseismiclog;
        std::vector<ExpEntry>  postseismicexp;
        std::vector<TanhEntry> ssetanh;
        bool                   seasonal_signal,halfseasonal_signal;
        bool                   estimate_multivariate;
        DesignMatrix(void);
        ~DesignMatrix(void);

        DesignMatrix(const DesignMatrix &) = delete;
        DesignMatrix& operator=(const DesignMatrix &) = delete;
 
        void compute_Amp(double Ac, double As,
                          double sigma_in, double& Amp, double& sigma_out);
        void compute_Pha(double Ac, double As, double sigma_Ac,
                          double sigma_As, double& Pha, double& sigma_out);

      public:
        //--- Normal pointer singleton
        static DesignMatrix* getInstance(void) {
          if (singleton==NULL) {
            singleton = new DesignMatrix();
          }
          return singleton;
        }
         
        //--- Destroy singleton
        static void resetInstance() {
          delete singleton;
          singleton = NULL; // so GetInstance will still work.
        }

        void    get_H(int& n_, double **H_); 
        void    get_F(double **F_); 
        int     get_offset_index(void) {return index_offset;};
        void    show_results(double *theta, double *error);
        void    compute_xhat(const double *theta);
    };

  #endif

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
