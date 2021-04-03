/*! \file     Statistics.cpp
 *  \author   Machiel Bos
 *  \version  1.1
 *
 * Apply some simple statistical tests to the observations/residuals.
 *
 * \date  29/1/2013  Santa Clara
 */
//==============================================================================
  #include "Statistics.h"
  #include "Observations.h"
  #include <iostream>
  #include <ostream>
  #include <cmath>
  #include <cstdlib>
  #include <cstdio>
  #include <gsl/gsl_statistics.h>
  #include <gsl/gsl_sort.h>
  #include <gsl/gsl_cdf.h>

//==============================================================================
// Subroutines
//==============================================================================


/* Durbin-Watson test for serial correlation
 */
//-----------------------------------
  void Statistics::DurbinWatson(void)
//-----------------------------------
  {
    int            i,m;
    double         numerator,denominator,*x,*t;
    Observations   *observations=Observations::getInstance();

    using namespace std;
    //--- If there is estimated signal, compute test on raw observations
    denominator = numerator = 0.0;
    observations->get_values(m,&t,&x);
    for (i=1;i<m;i++) {
      if (!(isnan(x[i]) || isnan(x[i-1]))) {
        numerator   += pow(x[i]-x[i-1],2.0);
        denominator += pow(x[i],2.0);
      }
    }

    //--- Free memory
    delete[] t;
    delete[] x;

    cout << "Durbin-Watson statistics is:" << numerator/denominator << endl;
  }



/* Create file for Q-Q plot
 */
//-----------------------------
  void Statistics::QQplot(void)
//-----------------------------
  {
    int            i,n;
    double         *prob,*x,*y,*t,mean,sigma;
    Observations   *observations=Observations::getInstance();
    FILE           *fp;

    using namespace std;
    //--- Read observations, after removing missing data and sort it
    observations->remove_gaps();
    observations->get_values(n,&t,&x);
    mean  = gsl_stats_mean(x,1,n);
    sigma = gsl_stats_sd_m(x,1,n,mean);
    for (i=0;i<n;i++) x[i] -= mean;
    gsl_sort(x,1,n);

    //--- Allocate memory space
    prob = new double[n];
    y    = new double[n];

    //--- Compute inverse of normal distribution.
    for (i=0;i<n;i++) {
      prob[i] = (static_cast<double>(i)+1.0)/static_cast<double>(n+1);
      y[i] = gsl_cdf_gaussian_Pinv(prob[i],sigma);
    }

    //--- Write x and y to file
    fp = fopen("QQplot.dat","w");
    for (i=0;i<n;i++) {
      fprintf(fp,"%12.5le  %12.5le\n",x[i],y[i]);
    }
    fclose(fp);

    //--- Free memory
    delete[] t;
    delete[] x;
    delete[] y;
    delete[] prob;
  }



/*! call various statistical tests
 */
//---------------------------
  void Statistics::show(void)
//---------------------------
  {
    DurbinWatson();
    QQplot();
  }
