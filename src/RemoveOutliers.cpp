/*! \file    RemoveOutliers.cpp
 *  \author  Machiel Bos
 *  \version 1.0
 *
 * A C++ version of my Gnu-Octave script to remove outliers using data
 * snooping and throwing everything away that falls outside 3 times the
 * quartile range.
 *
 * \date 6/7/2012  Santa Clara
 */
//==============================================================================
  #include <cstdlib>
  #include <cmath>
  #include <iostream>
  #include <ostream>
  #include <algorithm>
  #include "RemoveOutliers.h"

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!-------------------------------
  RemoveOutliers::RemoveOutliers(void)
//---!!-------------------------------
  {
    Control   &control=Control::getInstance();
   
    using namespace std;
    try { 
      factor = control.get_double("IQ_factor");
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
  }


 
/* Compute least-squares
 */
//-------------------------------------------------------------
  void RemoveOutliers::compute_LeastSquares(int& m, double **r)
//-------------------------------------------------------------
  {
    Observations   &observations=Observations::getInstance();
    DesignMatrix   *designmatrix=NULL;
    int            n;
    double         *t,*x,*H,*theta,*Ctheta;
#ifdef DEBUG
    int            i,j;
    FILE           *fp;
#endif

    using namespace std;
    //--- get observations and design matrix H
    observations.remove_gaps();
    observations.get_values(m,&t,&x);
    cout << "after removal gaps, m=" << m << endl;
    designmatrix = DesignMatrix::getInstance();
    designmatrix->get_H(n,&H);

#ifdef DEBUG
    fp = fopen("x.dat","w");
    for (i=0;i<m;i++) {
      fprintf(fp,"%e\n",x[i]);
    }
    fclose(fp);
    fp = fopen("H.dat","w");
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) fprintf(fp,"%e ",H[i+j*m]);
      fprintf(fp,"\n");
    }
    fclose(fp);
#endif

    theta = new double[n];

    //--- Copy original x vector to r
    *r = new double[m];
    cblas_dcopy(m,x,1,*r,1);

    //--- Simulate dgels using ATLAS subroutines (ordinary least-squares)
    Ctheta = new double[n*n];
    memset(Ctheta,0,(n*n)*sizeof(double)); // put Ctheta to zero
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,n,m,1.0,H,m,0.0,Ctheta,n);
    clapack_dpotrf(CblasColMajor,CblasUpper,n,Ctheta,n);
    cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,H,m,x,1,0.0,theta,1); 
    clapack_dpotrs(CblasColMajor,CblasUpper,n,1,Ctheta,n,theta,n); 

    //--- Compute r = xhat - x (matrix dummy contains original A matrix)
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
					m,1,n,1.0,H,m,theta,n,-1.0,*r,m);

#ifdef DEBUG
    fp = fopen("theta.dat","w");
    for (i=0;i<n;i++) {
      fprintf(fp,"%e\n",theta[i]);
    }
    fclose(fp);
    fp = fopen("r.dat","w");
    for (i=0;i<m;i++) {
      fprintf(fp,"%e\n",(*r)[i]);
    }
    fclose(fp); 
#endif

    //--- Make sure a new design matrix will be constucted
    designmatrix->resetInstance();
       
    //--- free memory
    delete[] theta;
    delete[] x;
    delete[] t;
    delete[] H;
    delete[] Ctheta;
  }



/*! Needed for qsort
 */
//------------------------------------------------------------
  int RemoveOutliers::compare (const void * a, const void * b)
//------------------------------------------------------------
  {
    using namespace std;
    if (isnan(*(double*)a) || isnan(*(double*)b)) {
      cerr << "Something wrong here!" << endl;
      exit(EXIT_FAILURE);
    } else { 
      return static_cast<int>( *(double*)a - *(double*)b );
    }
  }



/*! Data snooping
 */
//----------------------------------------
  void RemoveOutliers::data_snooping(void)
//----------------------------------------
  {
    using namespace std;
    const double   NaN=sqrt(-1.0);
    Observations   &observations=Observations::getInstance();
    int            m,bad_points,i;
    vector<double> q;
    double         *r,interquartile,median,lower_boundary,upper_boundary;

    do {
      q.clear();
      compute_LeastSquares(m,&r);
      for (i=0;i<m;i++) q.push_back(r[i]);
      sort(q.begin(),q.end());
#ifdef DEBUG
      for (i=0;i<m;i++) cout << "i=" << i << ", r=" << r[i] 
						<< ", q=" << q[i] << endl;
#endif
      interquartile  = q[int(0.75*m)] - q[int(0.25*m)];
      median         = q[int(0.50*m)];
      lower_boundary = median - factor*interquartile;
      upper_boundary = median + factor*interquartile;
#ifdef DEBUG
      cout << "median = " << median << ", lower_boundary=" << lower_boundary
           << ", upper_boundary=" << upper_boundary << endl;
#endif

      bad_points=0;
      for (i=0;i<m;i++) {
        if (r[i]<lower_boundary || r[i]>upper_boundary) {
          //cout << "bad point index = " << i << ", value=" << r[i] << endl;
          observations.set_one_x(i,NaN);
          bad_points++;
        }
      }
      cout << "Found " << bad_points << " bad points." << endl;
      
      //--- free memory for next run
      delete[] r;
      q.clear();
    
    } while (bad_points>0);
  }



//==============================================================================
// Main Program
//==============================================================================

  int main(int argc, char *argv[])
  {
    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      Control &control = Control::getInstance("removeoutliers.ctl");
      (void) control;
    } else if (argc==2) {
      Control &control = Control::getInstance(argv[1]);
      (void) control;
    } else {
      cerr << "correct usage: removeoutliers [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    cout << endl
         << "************************************" << endl
         << "    removeoutliers, version " << VERSION << endl
         << "************************************" << endl;

    //--- Now it's safe to call the other classes
    RemoveOutliers  outliers;
    Observations    &observations=Observations::getInstance();

    //--- Start removeoutliers
    outliers.data_snooping();
    observations.save_mom(true);

    return EXIT_SUCCESS;
  }
