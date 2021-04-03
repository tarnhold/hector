/*! \file   HyperGeo.h
 *  \author Machiel Bos
 *
 * Header file for HyperGeo.cpp
 *
 *  This script is part of Hector 1.7.2
 *
 *  Hector is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  Hector is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Hector. If not, see <http://www.gnu.org/licenses/>
 *
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
