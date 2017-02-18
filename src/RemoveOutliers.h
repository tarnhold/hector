/*! \filename RemoveOutliers.h
 *  \author   Machiel Bos
 *
 * \date  6/7/2012  Santa Clara
 * \date 30/9/2015  Santa Clara
 */
//==============================================================================

  #ifndef REMOVEOUTLIERS__
    #define REMOVEOUTLIERS__
    #include "Observations.h"
    #include "DesignMatrix.h"
    extern "C" {
      #include "cblas.h"
      #include "clapack.h"
    };

    class RemoveOutliers
    {
      private:
        int      SegmentLength;
        double   factor,StudentT;
       
      public:
        RemoveOutliers(void);
        void     compute_LeastSquares(int& m, double **r);
        int      compare (const void * a, const void * b);
        void     data_snooping(void);
    };

  #endif 

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
