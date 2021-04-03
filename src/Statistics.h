/*! \filename Statistics.h
 *  \author   Machiel Bos
 *  \version 1.0
 *
 * Header file for Statistics.cpp
 *
 * \date 18/3/2013  Santa Clara
 */
//==============================================================================

  #ifndef __STATISTICS
    #define __STATISTICS
  
    class Statistics
    {
      private:

      public:
        void    DurbinWatson(void);
        void    QQplot(void);
        void    show(void);
    };

  #endif
