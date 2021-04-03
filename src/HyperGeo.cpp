/*! \file   HyperGeo.cpp
 *  \author Machiel Bos
 *
 * The Gnu Scientific Library does not provide the Gauss Hypergeometric
 * function with complex arguments. This is a tricky function but fortunately
 * I found the MSc dissertation by John Pearson, University of Oxford, on
 * this topic which comes with good Matlab scripts. I have converted these
 * into C++ and applied his selection scheme using transformations to 
 * make sure the right method is chosen.
 *
 * \date 10/10/2012  Santa Clara
 *
 * References:
 * John Pearson (2009): "Computation of Hypergeometric Functions, Master 
          Thesis, University of Oxford. 
 */
//==============================================================================
  #include "HyperGeo.h"
  #include <cmath>
  #include <cstdlib>
  #include <iostream>
  #include <iomanip>
  #include <ostream>

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!---------------------------------------
  HyperGeo::HyperGeo(void) : tol(1.0e-10),
                             TINY(1.0e-8),
                             rho(0.9),
                             pi(4.0*atan(1.0))
//---!!---------------------------------------
  {
  }



/*! Check if complex number z is an integer along the real-axis
 */
//-------------------------------------------------
  bool HyperGeo::is_integer(std::complex<double> z)
//-------------------------------------------------
  {
    using namespace std;
    if (fabs(imag(z))<TINY && fabs(fmod(real(z),1.0))<TINY) return true;
    else                                                    return false;
  }



/*! A little wrapper to  gsl_sf_lngamma_complex_e
 */
//--------------------------------------------------------------
  std::complex<double> HyperGeo::lngamma(std::complex<double> z)
//--------------------------------------------------------------
  {
    using namespace std;
    gsl_sf_result    lnr,arg;
    complex<double>  w;

    gsl_sf_lngamma_complex_e(real(z),imag(z),&lnr,&arg);
    w = complex<double>(lnr.val*cos(arg.val),lnr.val*sin(arg.val));  

    return w;
  }



/*! Computes Gauss Hypergeometric function 2F1 using the single fraction
 *  method. Adapted from John Pearson's Masters dissertation Matlab script.
 *  He says it works well for:
 *   |a|,|b|<50  |c|>1  |z|<0.9
 *   |a|,|b|<30  |c|<1  |z|<0.9
 */
//----------------------------------------------------------------------------
  std::complex<double> HyperGeo::singlefraction2F1(std::complex<double> a,
       std::complex<double> b, std::complex<double> c, std::complex<double> z)
//----------------------------------------------------------------------------
  {
    using namespace std;
    int               i;
    bool              convergence_reached;
    const int         max_iter = 500;
    complex<double>   a1[max_iter],b1[max_iter],c1[max_iter],d1[max_iter],I;

    a1[0] = 0.0;    a1[1] = c;
    b1[0] = 1.0;    b1[1] = a*b*z;
    c1[0] = 1.0;    c1[1] = c;
    d1[0] = 1.0;    d1[1] = (c+a*b*z)/c;

    i = 2;
    I = 3.0;
    do {
      a1[i] = (a1[i-1]+b1[i-1])*(I-1.0)*(c+I-2.0);
      b1[i] = b1[i-1]*(a+I-2.0)*(b+I-2.0)*z;
      c1[i] = c1[i-1]*(I-1.0)*(c+I-2.0);

      if (isinf(abs(a1[i])) || isinf(abs(b1[i])) || isinf(abs(c1[i]))) {
        throw "HyperGeo::singlefraction2F1 - a1,b1,c1 become infinite";
        return NAN;
      } 

      d1[i]=(a1[i]+b1[i])/c1[i];

      convergence_reached = abs(d1[i]-d1[i-1])/abs(d1[i-1])<tol &&
                                abs(d1[i-1]-d1[i-2])/abs(d1[i-2])<tol;
      i++;
      I+=1.0;
    } while (i<max_iter && !convergence_reached);

    if (i==max_iter) {
      throw "HyperGeo::singlefraction2F1 - too many iterations";
      return NAN;
    }

    return d1[i-1];
  }



/*! A Runka-Kutta fourth order integrator
 */
//---------------------------------------------------------------------------
  std::complex<double> HyperGeo::deivprk2f1(std::complex<double> a,
       std::complex<double> b, std::complex<double> c, 
					       std::complex<double> z, int n)
//---------------------------------------------------------------------------
  {
    using namespace std;
    //--- Simulate nested fuctions by creating a local class
    class Local {
      public:
        complex<double>f1(complex<double>u,complex<double>v,complex<double>w) {
          return w;
        }
        complex<double>f2(complex<double> u,complex<double>v,complex<double>w,
                          complex<double> a,complex<double>b,complex<double>c){
          complex<double>   A1,A2,A3;

          A1 = u*(1.0-u);
          A2 = c-(a+b+1.0)*u;
          A3 = -a*b;
          return (-1.0/A1*(A2*w+A3*v));
        };
    };

    Local              local;
    int                i;
    complex<double>    x[n+1],y1[n+1],z1[n+1];
    complex<double>    h,k1y,k1z,k2y,k2z,k3y,k3z,k4y,k4z,dydx;

    //--- compute stepsize
    h = z/static_cast< complex<double> >(n);

    //--- derivative at x=0.0
    dydx = a*b/c;

    //--- Initial values
    y1[0] = 1.0;
    z1[0] = dydx;
    y1[1] = singlefraction2F1(a,b,c,h);
    z1[1] = dydx*singlefraction2F1(a+1.0,b+1.0,c+1.0,h);

    x[0] = 0.0;
    x[1] = h;
    for (i=1;i<n;i++) {
      k1y=local.f1(x[i],y1[i],z1[i]);
      k1z=local.f2(x[i],y1[i],z1[i],a,b,c);
      k2y=local.f1(x[i]+0.5*h,y1[i]+0.5*k1y*h,z1[i]+0.5*k1z*h);
      k2z=local.f2(x[i]+0.5*h,y1[i]+0.5*k1y*h,z1[i]+0.5*k1z*h,a,b,c);
      k3y=local.f1(x[i]+0.5*h,y1[i]+0.5*k2y*h,z1[i]+0.5*k2z*h);
      k3z=local.f2(x[i]+0.5*h,y1[i]+0.5*k2y*h,z1[i]+0.5*k2z*h,a,b,c);
      k4y=local.f1(x[i]+h,y1[i]+k3y*h,z1[i]+k3z*h);
      k4z=local.f2(x[i]+h,y1[i]+k3y*h,z1[i]+k3z*h,a,b,c);
      x[i+1] =x[i]+h;
      y1[i+1]=y1[i]+h/6.0*(k1y+2.0*k2y+2.0*k3y+k4y);
      z1[i+1]=z1[i]+h/6.0*(k1z+2.0*k2z+2.0*k3z+k4z);
    }
#ifdef DEBUG
    cout << scientific << setprecision(15) << "y1=" << y1[i] << endl;
#endif

    return y1[n];
  }



/*! RK4 method which checks convergence
 */
//----------------------------------------------------------------------------
  std::complex<double> HyperGeo::RK4(std::complex<double> a,
       std::complex<double> b, std::complex<double> c, std::complex<double> z)
//----------------------------------------------------------------------------
  {
    using namespace std;
    const int         max_iter=13;
    int               i,n=250;
    complex<double>   h=NAN,h_old;

    h = deivprk2f1(a,b,c,z,n);
    i = 0;
    do {
      n *= 2; //--- double number of grid points : half stepsize
      i++;
      h_old = h;
      h = deivprk2f1(a,b,c,z,n);
    } while (i<max_iter && abs((h_old-h)/h)>tol);

    if (i==max_iter) {
      cerr << "HyperGeo::compute - too many iterations in RK4" << endl;
      exit(EXIT_FAILURE);
    } 

    return h;
  }



/*! Taylor method B
 */
//----------------------------------------------------------------------------
  std::complex<double> HyperGeo::TaylorB(std::complex<double> a,
       std::complex<double> b, std::complex<double> c, std::complex<double> z)
//----------------------------------------------------------------------------
  {
    using namespace std;
    const int          Nmax=250000;
    int                j;
    complex<double>    J,r[Nmax],A[Nmax];

    //--- Initialise r(j) as detailed in Section 4.2
    r[0]=a*b/c;
    r[1]=(a+1.0)*(b+1.0)/2.0/(c+1.0);

    //--- Initialise A(j) as detailed in Section 4.2
    A[0]=1.0+z*r[0];
    A[1]=A[0]+z*z*a*b/c*r[1];

    j=2;
    J=complex<double>(2.0,0.0);
    while (j<Nmax) {
      //--- Update r(j) and A(j) in terms of previous values
      r[j]=(a+J)*(b+J)/(J+1.0)/(c+J);
      A[j]=A[j-1]+(A[j-1]-A[j-2])*r[j]*z;
      //--- If stopping criterion is satisfied, terminate computation
      if (abs(A[j]-A[j-1])/abs(A[j-1])<tol && 
				abs(A[j-1]-A[j-2])/abs(A[j-2])<tol) {
        break;
      }
      j++;
      J += 1.0;
    }

    //-- If 500 terms have been computed before stopping criterion has been
    //   satisfied, state this
    if (j==Nmax) {
      cerr << "a=" << a << ", b=" << b << ", c=" << c << ", z=" << z << endl;
      throw "HyperGeo::TaylorB - too many iterations!";
      return NAN;
    } else {
      return A[j];
    }
  }
    
    

/*! Try to compute 2F1 with singlefraction first, if that fails then try RK4.
 */
//-----------------------------------------------------------------------------
  std::complex<double> HyperGeo::_2F1(std::complex<double> a,
        std::complex<double> b, std::complex<double> c, std::complex<double> z)
//-----------------------------------------------------------------------------
  {
    using namespace std;
    complex<double>   h=NAN;

    if (abs(z)<rho) {
      if ((abs(c)<1.0      && abs(a)<50.0 && abs(b)<50.0) ||
          (abs(c)>1.0-TINY && abs(a)<30.0 && abs(b)<30.0 && abs(c)<50.0)) {
        try {
          h = singlefraction2F1(a,b,c,z);
        } catch (const char * str) {
#ifdef DEBUG
          cerr << str << endl;
          cerr << "HyperGeo::_2F1 - Now trying RK4 method..." << endl;
#endif
          try {
            h = RK4(a,b,c,z);
          } catch (const char * str) {
            cerr << str << endl;
            cerr << "HyperGeo:: 1) Am giving up..." << endl;
            exit(EXIT_FAILURE);
          }
        }
      } else {
        cerr << "Am giving up..." << endl;
        exit(EXIT_FAILURE);
      }
    } else {
      cerr << "HyperGeo::_2F1 - abs(z)>" << rho << "!" << endl;
      exit(EXIT_FAILURE);
    }

    return h;
  }



/*! Compute 2F1 but apply first transformation if z is to far away from
 *  zero.
 */
//----------------------------------------------------------------------------
  std::complex<double> HyperGeo::compute(std::complex<double> a,
       std::complex<double> b, std::complex<double> c, std::complex<double> z)
//----------------------------------------------------------------------------
  {
    using namespace std;
    complex<double>   h,fraction1,fraction2;

    //--- Check for special cases
    if (abs(z)<TINY) {
      return 1.0;
    }

    //--- For large values of a,b or c it is better to use the Taylor
    //    method which is simply computing the summation. However, it 
    //    does not like all transformation. Probably because c gets smaller
    //    than a or b. Numerical Recipes get away with Runge Kutta while I
    //    don't! Must find out why. Is their Bulirsch–Stoer algorithm is that 
    //    good?
    if (abs(a)>30.0-TINY || abs(b)>30.0-TINY || abs(c)>30.0-TINY) {
      try {
        h = TaylorB(a,b,c,z);
      } catch (const char * str) {
        cerr << str << endl;
        exit(EXIT_FAILURE);
      }

    //--- Try single fraction and RK4
    } else {

      //--- No transformation needed
      if (abs(z)<rho) {
        h = _2F1(a,b,c,z);

      //--- 1-z closer to zero
      } else if (abs(1.0-z)<rho && is_integer(c-a-b)==false) {
        fraction1 = exp(lngamma(c)+lngamma(c-a-b)-lngamma(c-a)-lngamma(c-b));
        fraction2 = exp(lngamma(c)+lngamma(a+b-c)-lngamma(a)-lngamma(b));
        h = fraction1*_2F1(a,b,a+b-c+1.0,(1.0-z)) +
                pow(1.0-z,c-a-b)*fraction2*_2F1(c-a,c-b,c-a-b+1.0,1.0-z);
        //--- Check if argument of z is within permissable range
        if (fabs(arg(1.0-z))>pi-TINY) h = conj(h); 

      //--- z/(z-1) is closer to zero
      } else if (abs(z/(z-1.0))<rho) {
        h = pow(1.0-z,-a)*_2F1(a,c-b,c,z/(z-1.0));

      //--- 1/(1-z) is closer to zero
      } else if (abs(1.0/(1.0-z))<rho && is_integer(a-b)==false) {
        fraction1 = exp(lngamma(c)+lngamma(b-a)-lngamma(b)-lngamma(c-a));
        fraction2 = exp(lngamma(c)+lngamma(a-b)-lngamma(a)-lngamma(c-b));
        h = pow(1.0-z,-a)*fraction1*_2F1(a,c-b,a-b+1.0,1.0/(1.0-z)) +
                   pow(1.0-z,-b)*fraction2*_2F1(b,c-a,b-a+1.0,1.0/(1.0-z)); 
        //--- Check if argument of z is within permissable range
        if (fabs(arg(1.0-z))>pi-TINY) h = conj(h);

      //--- 1-1/z closer to zero
      } else if (abs(1.0-1.0/z)<rho && is_integer(c-a-b)==false) {
        fraction1 = exp(lngamma(c)+lngamma(c-a-b)-lngamma(c-a)-lngamma(c-b));
        fraction2 = exp(lngamma(c)+lngamma(a+b-c)-lngamma(a)-lngamma(b));
        h = pow(z,-a)*fraction1*_2F1(a,a-c+1.0,a+b-1.0,1.0-1.0/z) +
                   pow(z,a-c)*fraction2*_2F1(c-a,1.0-a,c-a-b+1.0,1.0-1.0/z);
        //--- Check if argument of z is within permissable range
        if (fabs(arg(z))>pi-TINY || fabs(arg(1.0-z))>pi-TINY) h = conj(h);

      //--- We're unlucky
      } else {
        cerr << "HyperGeo::compute - Another method is needed" << endl;
        exit(EXIT_FAILURE);
      }
    }

    return h;
  }
