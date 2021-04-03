/*! \file    Control.h
 *  \author  Machiel Bos
 *
 * A small class for reading the parameters that control the behaviour of
 * EstimateTrend. 
 *
 * \date 12/1/2012   CIIMAR, Porto
 */
//============================================================================

//============================================================================
// Class definition
//============================================================================

#ifndef __CONTROL
  #define __CONTROL
  #include <cstdio>
  #include <string>

  /*! Following the Namelist tradition of FORTRAN, this class lets you define
   *  parameter values in a ASCII file which is read at run-time. 
   */
  class Control 
  {
    private:
      FILE     *fp_in;
      char     ctl_file[100];
      void     find_label(const char label[]);
      void     open_file(void);
      void     close_file(void);

    public:
      static Control* getInstance(std::string filename = "estimatetrend.ctl");
      Control(std::string file);
      ~Control(void);
      void     get_name_list(const char label[], char **value, int& n);
      void     get_string(const char label[], char value[]);
      bool     get_bool(const char label[]);
      int      get_int(const char label[]);
      double   get_double(const char label[]);
  };

#endif
