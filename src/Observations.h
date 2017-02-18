/*! \file    Observations.h
 *  \author  Machiel Bos
 *
 * Header file for Observations.cpp
 *
 * \date 27/7/2011  CIIMAR, Porto
 */
//==============================================================================

  #ifndef __OBSERVATIONS
    #define __OBSERVATIONS
    #include <vector>
    #include <cstring>
    #include "Control.h"

    typedef struct {
      double  MJD,T;
    } LogEntry;

    typedef struct {
      double  MJD,T;
    } ExpEntry;

    class Observations
    {
      private:
        const double          NaN;
        int                   Ngaps;
        double                fs,scale_factor;
        bool                  PSMSL_monthly,interpolate_data;
        void                  (Observations::*read_observations)(std::string 
								     filename);
        std::vector<double>   t,x,xhat,offsets,breaks;
        std::vector<LogEntry> postseismiclog;
        std::vector<ExpEntry> postseismicexp;
        std::string 	      extension;
 
        Observations(void);
        void   read_header(std::fstream& fp, int component);
        void   read_external_header(std::fstream& fp);
        void   clean_offsets(void);
        void   write_header(std::fstream& fp);
        void   read_PSMSL_monthly(std::string filename);
        void   read_mom(std::string filename);
        void   read_enu(std::string filename);
        void   read_neu(std::string filename);
        void   determine_fs(void);

      public:
        //--- Meyers singleton
        static Observations& getInstance(void) {
          static Observations theObservations;
          return theObservations;
        }

        void   save_mom(bool header=false);
        void   get_values(int& m, double **t_, double **x_);
        double get_fs(void) {return fs;};
        double estimate_lag1(void);
        void   make_continuous(bool interpolate);
        void   remove_gaps(void);
        int    number_of_gaps(void) {return Ngaps;};
        int    number_of_offsets(void) {return offsets.size();};
        int    number_of_observations(void) {return t.size();};
        void   set_xhat(const double *xhat_);
        void   set_one_x(int index, double value);
        void   get_offsets(std::vector<double>& offsets_);
        void   get_breaks(std::vector<double>& breaks_);
        void   get_postseismiclog(std::vector<LogEntry>& postseismiclog_);
        void   get_postseismicexp(std::vector<ExpEntry>& postseismicexp_);
        void   add_offset(double MJD);
        void   change_offset(size_t column, double MJD);
    };

  #endif
