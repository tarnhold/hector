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

    class Observations
    {
      private:
        static bool         instanceFlag;
        static Observations *singleton;

        const double        NaN;
        int                 Ngaps;
        double              fs,scale_factor;
        bool                PSMSL_monthly,interpolate_data,first_difference;
        void                (Observations::*read_observations)(char filename[]);
        std::vector<double> t,x,xhat,offsets;
        char                name[80],extension[20];
 
        void   read_header(std::fstream& fp, int component);
        void   read_external_header(std::fstream& fp);
        void   clean_offsets(void);
        void   write_header(std::fstream& fp);
        void   read_PSMSL_monthly(char filename[]);
        void   read_mom(char filename[]);
        void   read_enu(char filename[]);
        void   read_neu(char filename[]);
        void   read_pos(char filename[]);
        void   determine_fs(void);

      public:
        Observations(void);
        static Observations* getInstance(void);
        void   save_mom(bool header=false);
        void   get_values(int& m, double **t_, double **x_);
        void   get_name(char name_[]) {strcpy(name_,name);};
        double get_fs(void) {return fs;};
        double estimate_lag1(void);
        void   make_continuous(bool interpolate);
        void   remove_gaps(void);
        void   take_first_difference(void);
        int    number_of_gaps(void) {return Ngaps;};
        int    number_of_observations(void) {return t.size();};
        void   set_xhat(const double *xhat_);
        void   set_one_x(int index, double value);
        void   get_offsets(std::vector<double>& offsets_);
    };

  #endif
