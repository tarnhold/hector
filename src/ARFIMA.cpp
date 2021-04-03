/*! \file   ARFIMA.cpp
 *  \author Machiel Bos
 *
 *  Implementation of the ARFIMA noise model using the tricks of Sowell (1992),
 *  Doornik and Ooms (2003) and Zinde-Wash (1988).
 *
 *  \date 10/10/2012   Santa Clara
 */
//=============================================================================
  #include "ARFIMA.h"
  #include <gsl/gsl_poly.h>
  #include <cmath>
  #include <math.h>
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstring>

  #define TINY 1.0e-7
//  #define DEBUG

//=============================================================================
// Subroutines
//=============================================================================


/*! Read the p and q values from the control file. 
 *
 * \param[in] d_fixed_  fractional noise parameter
 */
//---!!--------------------------
  ARFIMA::ARFIMA(double d_fixed_)
//---!!--------------------------
  {
    Control   *control = Control::getInstance();

    using namespace std;
    //--- See how many AR and MA coefficients we require
    try {
      p = control->get_int("AR_p");
      q = control->get_int("MA_q");
    }
    catch (const char* str) {
      cerr << str << endl;
      exit(EXIT_FAILURE);
    }

    //--- See if we need to set d to zero, creating an ARMA model. In theory
    //    we can set d to any value between -0.5 and 0.5 but only 0 makes
    //    practical sense.
    if (isnan(d_fixed_)==false) {
      estimate_spectral_index = false;
      d_fixed = d_fixed_;
      Nparam  = p+q; // include 'd' in the counting of parameters
    } else {
      estimate_spectral_index = true;
      Nparam  = p+q+1; // include 'd' in the counting of parameters
    }

#ifdef DEBUG
    cout << "p=" << p << ", q=" << q << ", Nparam=" << Nparam << endl;
#endif

    //--- Set following arrays to NULL to avoid problems later on
    psi_PSD = NULL;
    rho_PSD = NULL;
  }



//---!!----------------
  ARFIMA::~ARFIMA(void)
//---!!----------------
  {
    if (psi_PSD!=NULL)   delete[] psi_PSD;
    if (rho_PSD!=NULL)   delete[] rho_PSD;
  }



/*! Use the GSL to find roots of the charactistic polynomial
 *  I follow Sowell's definition: Phi(L) = 1 + phi_1L + phi_2L^2 + ..
 *
 * \param[in]   AR    array of p autoregressive coefficients
 * \param[out]  rho   array of p distinct roots (complex numbers)
 */
//--------------------------------------------------------------
  void ARFIMA::find_roots(double *AR, std::complex<double> *rho)
//--------------------------------------------------------------
  {
     int                         i;
     double                      roots[2*p],coeff[p+1];
     gsl_poly_complex_workspace  *w;
    
     using namespace std;
     if (p==0) {
       rho = NULL;
     } else if (p==1) {
       rho[0] = complex<double>(AR[0],0.0);
     } else {
       //--- find roots  
       coeff[0] = 1.0;
       for (i=0;i<p;i++) coeff[1+i] = -AR[i];
       if (fabs(coeff[p])<TINY) {
         w = gsl_poly_complex_workspace_alloc (p);
         gsl_poly_complex_solve (coeff, p, w, roots);
       } else {
         w = gsl_poly_complex_workspace_alloc (p+1);
         gsl_poly_complex_solve (coeff, p+1, w, roots);
       }
     
       for (i=0;i<p;i++) {
         if (i==p-1 && fabs(coeff[p])<TINY) {
           rho[i] = 0.0;
         } else {
           rho[i] = 1.0/complex<double>(roots[2*i],roots[2*i+1]);
         }
         //--- Check if stationarity condition has been met
         if (abs(rho[i])>1.0) {
           cerr << "root " << i << " is larger than 1 : " << rho[i] << endl;
       //    exit(EXIT_FAILURE);
         }
       }
       //--- free memory
       gsl_poly_complex_workspace_free (w);
     }
  }



/*! Find coefficients out of given roots. I played with the idea to use
 *  the roots in the minimization process instead of the AR coefficients but
 *  rejected it because I cannot handle complex parameters. Therefore, I need
 *  an easy way to convert back and forth between roots and AR coefficients.
 *
 * \param[in]  rho   array of p distinct roots (complex numbers)
 * \param[out] AR    array of p autoregressive coefficients
 */
//---------------------------------------------------------------------
  void ARFIMA::find_coefficients(std::complex<double> *rho, double *AR)
//---------------------------------------------------------------------
  {
    using namespace std;
    int               i,j,k,n;
    complex<double>   coeff,*dummyC;

    //--- complex numbers need to be stored
    dummyC = new complex<double>[p];
    for (i=0;i<p;i++) dummyC[i]=0.0;

    //--- Compute the number of possible subsets
    n = 1;
    n <<= p;
    for (i=0;i<n;i++) {
      coeff = 1.0;
      k     = 0;
      for (j=0;j<p;j++) {
        if ((i & 1<<j)>>j) {
          coeff *= rho[j];
          k++;
        }
      }
      //--- Avoid the empty set solution 
      if (k>0) {
        dummyC[k-1] += coeff;
      }
    }
    for (i=0;i<p;i++) {
      AR[i] = pow(-1.0,static_cast<double>(i))*real(dummyC[i]);
    }
    delete[] dummyC;
  }



/*! Compute pure ARMA covariance matrix. NOTE: change of sign of MA 
 *  coefficients needed!
 */
//-----------------------------------------------------------------------
  void ARFIMA::ZindeWalsh(double *AR, double *MA, int m, double *gamma_x)
//-----------------------------------------------------------------------
  {
    using namespace std;
    int              i,j,k;
    double           *alpha;
    complex<double>  *rho,*zeta,*xi,g;

#ifdef DEBUG
    cout << "ZindeWalsh" << endl;
    for (i=0;i<p;i++) cout << "AR i:" << i << " ,  " << AR[i] << endl;
    for (i=0;i<q;i++) cout << "MA i:" << i << " ,  " << MA[i] << endl;
#endif
    //--- Compute alpha 
    alpha = new double[q+1];
    alpha[0] = 1.0;
    for (i=0;i<q;i++) alpha[0] += pow(-MA[i],2.0);
    alpha[0] *= 0.5;
    for (i=1;i<=q;i++) {
      alpha[i] = MA[i-1];
      for (j=0;j<q-i;j++) {
        alpha[i] += MA[j]*MA[j+i];
      }
    }

    //--- Compute zeta
    zeta = new complex<double>[p];
    rho  = new complex<double>[p];
 
    find_roots(AR,rho);
    for (i=0;i<p;i++) {
      zeta[i] = complex<double>(1.0,0.0);
      for (j=0;j<p;j++) {
        zeta[i] *= (1.0 - rho[i]*rho[j]);
        if (j!=i) zeta[i] *= (rho[i]-rho[j]);
      }
      //--- Already do inversion to get zeta[j]
      zeta[i] = pow(rho[i],static_cast<double>(p-1))/zeta[i];
    }
#ifdef DEBUG
    for (i=0;i<p;i++) cout << "rho i:" << i << " ,  " << rho[i] << endl;
    for (i=0;i<p;i++) cout << "zeta i:" << i << " ,  " << zeta[i] << endl;
#endif

    //--- Compute xi
    xi = new complex<double>[p];
    for (i=0;i<p;i++) {
      xi[i] = 0.0;
      for (j=0;j<=q;j++) {
        xi[i] += alpha[j]*(pow(rho[i],static_cast<double>(j)) +
				(pow(rho[i],static_cast<double>(-j))));
      }
      xi[i] *= zeta[i];
    }
#ifdef DEBUG
    for (i=0;i<p;i++) cout << "xi i:" << i << " ,  " << xi[i] << endl;
#endif

    //*** Finally compute gamma_x ******************
    if (p==0 && q>0) {  //--- pure MA
      gamma_x[0] = 2.0*alpha[0];
      for (i=1;i<=q;i++) gamma_x[i] = alpha[i];
      memset(&gamma_x[q+1],0.0,(m-q-1)*sizeof(double));
    } else if (q==0 && p>0) { //--- pure AR
      for (i=0;i<m;i++) {
        g = 0.0;
        for (j=0;j<p;j++) 
	  g += zeta[j]*pow(rho[j],static_cast<double>(i));
        gamma_x[i] = real(g);
      }
    } else if (q>0 && p>0) {
      for (i=0;i<m;i++) {
        g = 0.0;
        for (j=0;j<p;j++) 
          g += xi[j]*pow(rho[j],static_cast<double>(i));
        for (j=i+1;j<=q;j++) {
          for (k=0;k<p;k++) {
            g += alpha[j]*zeta[k]*(pow(rho[k],static_cast<double>(j-i))
                                         -pow(rho[k],static_cast<double>(i-j)));
          }
        }
        gamma_x[i] = real(g);
      }
    } else {
      cerr << "Both q and p are zero!" << endl;
      exit(EXIT_FAILURE);
    }
#ifdef DEBUG
//    for (i=0;i<m;i++) cout << "gamma i:" << i << " ,  " << gamma_x[i] << endl;
#endif
       
    //--- Free memory
    delete[] alpha;
    delete[] zeta;
    delete[] xi;
    delete[] rho;
  }  



/*! Compute psi of Doornik and Ooms out of the MA coefficients
 */
//-------------------------------------------------
  void ARFIMA::compute_psi(double *MA, double *psi)
//-------------------------------------------------
  {
    int     i,j,k;
    double  *theta;

    theta = new double[q+1];
    theta[0] = 1.0;
    for (i=0;i<q;i++) theta[i+1]=MA[i];
    for (k=-q;k<=q;k++) {
      psi[k+q] = 0.0;
      for (j=abs(k);j<=q;j++) {
        psi[k+q] += theta[j]*theta[j-abs(k)];
      }
    }
    delete[] theta;
  }



/*! Compute the covariance matrix for an ARFIMA model
 *
 * \param[in]  AR      : array containing p phi coefficients.
 * \param[in]  p       : number of phi coefficients.
 * \param[in]  d       : the value of the Fractional Integrated noise (-1<d<0.5)
 * \param[in]  MA      : array containing q theta coefficients.
 * \param[in]  m       : number of observations and thus covariance coefficients
 * \param[out] gamma_x : the computed covariance column
 */
//--------------------------------------------------------------------------
  void ARFIMA::DoornikOoms(double *AR, double d, double *MA,
						     int m, double *gamma_x)
//--------------------------------------------------------------------------
  {
    using namespace std;
    int                i,j,k,N,h_max;
    complex<double>    scale,*rho=NULL,*zeta=NULL,*G=NULL;
    double             offset,dummy,*psi=NULL,*FI_fraction=NULL;
    double	       **C=NULL,a_max,c_max;
    HyperGeo           hypergeo;

    //--- Compute FI_fraction  - FI = fractionally integrated
    h_max = -p+q+m-1;
    if (h_max<0) {
      cerr << "Very weird, -p+q+m-1 < 0!" << endl;
      exit(EXIT_FAILURE);
    }
    N = (p+q) + h_max + 1;
    FI_fraction = new double[N];

    //--- Compute fraction with Pochhammer symbols. If d==0 then the
    //    Gauss hypergeometric function becomes infinite for negative
    //    integer values.
    if (fabs(d)<TINY) { 
      cerr << "Oops! This should not happen... aborting!" << endl;
      exit(EXIT_FAILURE);
    } else {
      FI_fraction[h_max] = exp(gamma(1.0-2.0*d))/
					pow(exp(gamma(1.0-d)),2.0); //--- i=0
      //--- forward recurrence
      for (i=1;i<=(p+q);i++) {
        FI_fraction[h_max+i]   = FI_fraction[h_max+i-1]*
		(d+static_cast<double>(i-1))/(1.0-d+static_cast<double>(i-1));
      }
      //--- backward recurrence
      for (i=0;i>-h_max;i--) {
        FI_fraction[h_max+i-1] = FI_fraction[h_max+i]*
		(d-static_cast<double>(i))/(1.0-d-static_cast<double>(i));
      }
    }
#ifdef DEBUG
    for (i=0;i<m;i++) cout << "FI_fraction=" << FI_fraction[i] << endl;
#endif

    //--- Compute coefficients given by Sowell and Doornik & Ooms for AR part
    if (p>0) {
      //--- We also need the roots of the AR polynomial
      rho  = new complex<double>[p];
      find_roots(AR,rho);
      zeta = new complex<double>[p];

      //--- Compute zeta coefficients (ignore single rho_j in definition) which
      //    is absorbed by the 'G' function.
      for (j=0;j<p;j++) {
        zeta[j] = complex<double>(1.0,0.0);
        for (i=0;i<p;i++) zeta[j] *= (1.0 - rho[i]*rho[j]);
        for (k=0;k<p;k++) {
          if (k!=j) zeta[j] *= (rho[j]-rho[k]);
        }
        //--- Already do inversion to get zeta[j]
        zeta[j] = 1.0/zeta[j];
      }
#ifdef DEBUG
      for (i=0;i<p;i++) cout << "rho=" << rho[i] << ", zeta=" << zeta[i]<< endl;
#endif

      //--- Matrix C is NOT the covariance matrix but some useful matrix
      //    defined by Sowell. It runs over all p roots (distict) and 2*q+m
      //    observations (including MA filter overlap at both ends)
      C = new double* [p];
      G = new complex<double>[2*h_max+1];

      for (i=0;i<p;i++) {
        C[i] = new double[N];

        //--- for each i between [1,p], G runs from -h to h.
        a_max =  d + static_cast<double>(h_max); 
        c_max = -d + static_cast<double>(h_max) + 1.0; 

        //--- First, assume G=1/rho*(2F1(d+h_max1,1,-d+h_max,rho)-1) = 1.0
        G[2*h_max] = (hypergeo.compute(a_max,1.0,c_max,rho[i])-1.0)/rho[i];

        //--- Use stable backward recursion to compute all 2F1's
        for (j=2*h_max,offset=1.0;j>0;j--,offset+=1.0) {
          G[j-1] = (a_max-offset)/(c_max-offset)*(1.0+rho[i]*G[j]);
        }

#ifdef DEBUG
        cout << "2F1(" << d << ",1," << -d+1.0 << ";" << rho[i] << ")=" <<
			hypergeo.compute(d,1.0,-d+1.0,rho[i]) << endl;
        cout << "rho=" << rho[i] << ", d=" << d << ", scale="<< scale << endl;
#endif
        for (j=-h_max;j<=p+q;j++) {
          //--- Note that I here already include the zeta parameter in order
          //    to have a real (not complex) C.
          C[i][h_max+j] = FI_fraction[h_max+j]*real(zeta[i]*
 			         (pow(rho[i],2.0*p)*G[h_max+j] +
					 pow(rho[i],2.0*p-1) + G[h_max-j]));
        }
      }
    }

    //--- Compute psi coefficients
    psi   = new double[2*q+1];  //--- theta[0] = 1!
    compute_psi(MA,psi);

    //--- We have enough information: construct covariance matrix!
    if (p==0) { // q cannot be zero because otherwise we never reachs this
      for (i=0;i<m;i++) {
        gamma_x[i] = 0.0;
        for (k=-q;k<=q;k++) {
          gamma_x[i] += psi[k+q]*FI_fraction[h_max+k-i];
        }
      }
    } else {
      for (i=0;i<m;i++) {
        gamma_x[i] = 0.0;
        for (j=0;j<p;j++) {
          for (k=-q;k<=q;k++) {
            gamma_x[i] += psi[q+k]*C[j][h_max+p+k-i];
          }
        }
      }
    }
#ifdef DEBUG
    for (i=0;i<2*h_max+1;i++) {
      cout << "G[" << i << "]=" << G[i] << endl;
    }
    for (i=0;i<m;i++) {
      cout << "C[" << i << "]=" << C[0][i] << endl;
    }
#endif

    //--- Clean up mess
    if (rho!=NULL)     delete[] rho;
    if (zeta!=NULL)    delete[] zeta;
    if (psi!=NULL)     delete[] psi;
    if (p>0) {
      for (i=0;i<p;i++)  delete[] C[i];
      delete[] C;
      delete[] G;
    } 
    delete[] FI_fraction;
  }



/*! Interface with NoiseModel.cpp
 */
//------------------------------------------------------------------
  void ARFIMA::get_covariance(double *param, int m, double *gamma_x)
//------------------------------------------------------------------
  {
    double  d,*AR,*MA;
    int     i;

    //--- To make things more readable, rewrite param as AR,d,MA
    AR = &param[0];
    MA = &param[p];
    if (estimate_spectral_index==true) {
      d = param[p+q];
    } else {
      d = d_fixed;
    }

    //--- Construct ARFIMA covariance array
    if (fabs(d)<TINY) ZindeWalsh(AR,MA,m,gamma_x);
    else              DoornikOoms(AR,d,MA,m,gamma_x);
  }



/*! Parse parameter list and show results on screen
 *
 * \param[in] param : list of AR,d,MA
 */
//-----------------------------------------------
  void ARFIMA::show(double *param, double *error)
//-----------------------------------------------
  {
    int  i;
    Control  *control=Control::getInstance();

    using namespace std;
    for (i=0;i<Nparam;i++) {
      if (i<p) {
        printf("AR[%1d] = %8.4lf +/- %6.4lf\n",i+1,param[i],error[i]);
      } else if (i>=p && i<(p+q)) {
        printf("MA[%1d] = %8.4lf +/- %6.4lf\n",i+1-p,param[i],error[i]);
      } else if (i==(p+q)) {
        if (control->get_bool("firstdifference")==true) {
          printf("d     = %8.4lf +/- %6.4lf\n",param[i]+1.0,error[i]);
        } else {
          printf("d     = %8.4lf +/- %6.4lf\n",param[i],error[i]);
        }
      }
    }
  }



/*! Limit the range of values of valid ARFIMA parameters
 */
//---------------------------------------------
  double ARFIMA::compute_penalty(double *param)
//---------------------------------------------
  {
    using namespace std;
    int               i;
    double            penalty=0.0,LARGE=1.0e8;
    complex<double>   rho[p];

    //--- AR : I think create_ARFIMA does not choke on rho's<1
    if (p>0) {
#ifdef DEBUG
      for (i=0;i<p;i++) cout << "before i=" << i << ", " << param[i] << endl;
#endif
      find_roots(param,rho);
      for (i=0;i<p;i++) {
        if (abs(rho[i])>0.99) {
          penalty += (abs(rho[i])-0.99)*LARGE;
          rho[i] *= 0.99/abs(rho[i]);
        }
      }
      find_coefficients(rho,param);
#ifdef DEBUG
      for (i=0;i<p;i++) cout << "rho: i=" << i << ", " << rho[i] << endl;
      for (i=0;i<p;i++) cout << "after i=" << i << ", " << param[i] << endl;
#endif
    } 

    //--- d
    if (estimate_spectral_index==true) {
      if (param[p+q]>0.499) {
        penalty   += (param[p+q]-0.499)*LARGE;
        param[p+q] = 0.499;
      } else if (param[p+q]<-0.999) {
        penalty   += (-0.999-param[p+q])*LARGE;
        param[p+q] = -0.999;
      }
    }

    //--- MA does not have a contribution to the penalty value

#ifdef DEBUG
    cout << "Penalty: phi=" << param[0] << endl;
    if (estimate_spectral_index==true) cout<<"Penalty: d="<< param[p+q] << endl;
#endif

    return penalty; 
  }



/*! Store ARFIMA parameters
 */
//----------------------------
  void ARFIMA::setup_PSD(void)
//----------------------------
  {
    int     i;
    double  *AR,*MA;

    using namespace std;
    //--- allocate memory
    AR      = new double[p];
    rho_PSD = new complex<double>[p];
    MA      = new double[q];
    psi_PSD = new double[2*q+1];

    for (i=0;i<p;i++) {
      cout << "Enter parameter value of AR[" << i+1 << "]: ";
      cin >> AR[i];
    }
    for (i=0;i<q;i++) {
      cout << "Enter parameter value of MA[" << i+1 << "]: ";
      cin >> MA[i];
    }
    if (estimate_spectral_index==true) {
      cout << "Enter fractional difference value of d: ";
      cin >> d_PSD;
    } else {
      d_PSD = d_fixed;
    }

    //--- Now that parameters are known, compute roots of AR and psi
    find_roots(AR,rho_PSD);
    compute_psi(MA,psi_PSD);

    //--- clear up mess
    delete[] AR;
    delete[] MA;
  }



/*! Compute PSD for given frequency
 */
//---------------------------------------
  double ARFIMA::compute_G(double lambda)
//---------------------------------------
  {
    using namespace std;
    int              i;
    complex<double>  dummyC;
    double           G,dummy;

    //--- Compute MA part
    dummy = psi_PSD[0];
    for (i=1;i<=q;i++) {
      dummy += 2.0*psi_PSD[i]*cos(static_cast<double>(i)*lambda);
    }
    G = dummy;

    //--- Compute AR part
    dummyC = complex<double>(1.0,0.0);
    for (i=0;i<p;i++) {
      dummyC *= (1.0 - 2.0*rho_PSD[i]*cos(lambda) + pow(rho_PSD[i],2.0));
    }
    //--- If all goes well, this product produces a pure real number
    G /= real(dummyC);
    if (imag(dummyC)>1.0e-5) {
      cerr << "Somthing is wrong here: imag(dummyC)="<<imag(dummyC)<<endl;
      exit(EXIT_FAILURE);
    }

    //--- Add FI part
    G /= pow(2.0*sin(0.5*lambda),2.0*d_PSD);

    return G;
  }



/*! Compute impulse response: h
 */
//-------------------------------------------------------
  void ARFIMA::compute_impulse_response(int m, double* h)
//-------------------------------------------------------
  {
    int     i,j,k;
    double  *AR,*MA,*psi,*b,I;

    using namespace std;
    //--- Allocate memory space for psi and b
    psi = new double[m];
    b   = new double[m];
    AR  = new double[p];
    MA  = new double[q];

    //--- Ask for noise parameter values
    for (i=0;i<p;i++) {
      cout << "Enter parameter value of AR[" << i+1 << "]: ";
      cin >> AR[i];
    }
    for (i=0;i<q;i++) {
      cout << "Enter parameter value of MA[" << i+1 << "]: ";
      cin >> MA[i];
    }
    if (estimate_spectral_index==true) {
      cout << "Enter value of fractional difference d:";
      cin >> d_fixed;
    }

    //--- b is the standard power-law impulse response
    b[0] = 1.0;
    for (i=1,I=1.0;i<m;i++,I+=1.0) b[i] = (d_fixed+I-1.0)/I*b[i-1];

    //-- Write ARMA as an infinite MA model with psi coefficients
    psi[0] = 1.0;
    for (i=1;i<m;i++) {
      if (i<=q) psi[i] = MA[i-1];
      else      psi[i] = 0.0;
      if (i<=p) k=i;
      else      k=p;
      for (j=0;j<k;j++) {
        psi[i] += psi[i-j-1]*AR[j];
      }
    }

    //--- Apply Hassler and Kokoszka formula
    for (i=0;i<m;i++) {
      h[i] = 0.0;
      for (j=0;j<=i;j++) {
        h[i] += psi[j]*b[i-j];
      }
    }

#ifdef DEBUG
    for (i=0;i<m;i++) {
      cout << i << ", psi="<<psi[i]<<", b="<<b[i] << ", h=" << h[i] << endl;
    }
#endif

    //--- free memory space
    delete[] b;
    delete[] psi;
    delete[] AR;
    delete[] MA;
  }
