/*! \file   HyperGeo.h
 *  \author Machiel Bos
 *
 * Header file for HyperGeo.cpp
 *
 * \date 10/10/2012  Santa Clara
 */
//==============================================================================

  #ifndef __HYPERGEO
    #define __HYPERGEO
    #include <complex>
    #include <gsl/gsl_sf.h>

    class HyperGeo 
    {
      private:
         std::complex<double>  a,b,c,z;
         const double	       tol,TINY,rho,pi;

         bool                  is_integer(std::complex<double> z);
         std::complex<double>  singlefraction2F1(
				  std::complex<double> a,
        			  std::complex<double> b, 
				  std::complex<double> c, 
				  std::complex<double> z);
         std::complex<double>  deivprk2f1(
				  std::complex<double> a,
        			  std::complex<double> b, 
				  std::complex<double> c, 
				  std::complex<double> z, int n);
         std::complex<double>  RK4(
				  std::complex<double> a,
        			  std::complex<double> b, 
				  std::complex<double> c, 
				  std::complex<double> z);
         std::complex<double>  TaylorB(
				  std::complex<double> a,
        			  std::complex<double> b, 
				  std::complex<double> c, 
				  std::complex<double> z);
         std::complex<double>  _2F1(
				  std::complex<double> a,
        			  std::complex<double> b, 
				  std::complex<double> c, 
				  std::complex<double> z);
      public:
         HyperGeo(void);
         std::complex<double>  lngamma(std::complex<double> z);
         std::complex<double>  compute(
				  std::complex<double> a,
        			  std::complex<double> b, 
				  std::complex<double> c, 
				  std::complex<double> z);
    };

  #endif
