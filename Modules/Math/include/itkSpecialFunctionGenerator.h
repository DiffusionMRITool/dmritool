/**
 *       @file  itkSpecialFunctionGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-25-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSpecialFunctionGenerator_h
#define __itkSpecialFunctionGenerator_h

#include "utlCore.h"
#include "utlVector.h"
#include "utlMatrix.h"


namespace utl
{


/** lagurre polynomial \f$ L_n^a(x) \f$ */
template < class T >
T
Lagurre ( const int n, const double a, const T x );


/** \brief  gamma function.
 *
 *  \param  x real value.
 *        If x is positive integer, the result is factorial(x). 
 *        If x is a half of a positive integer, the result is analytical with gamma(0.5)=\sqrt(\pi). 
 *        If x is other real number, it returns the result of gsl_sf_gamma(). 
 **/
inline double 
Gamma(const double x);

inline double 
GammaLower (const double s, const double x);

/**
 * \brief bessel_Ja :  Regular Cylindrical Bessel Function  \f$ J_a(x) \f$
 *
 * \param a : should be an integer or half of an integer
 * \param x : float number
 * \note gsl_sf_bessel_Jn(a,x) in gsl only works when a is integer
 *
 * \return \f$ J_a(x) \f$
 */
inline double 
BesselJa(const double a, const double x);


/**
* @brief GetRotatedSHCoefficient 
*        rotate SH coefficient vector by a given rotation matrix
*
* @param shInput input SH coefficients
* @param rotationMatrix rotation matrix
* @author Jian Cheng
*
*  Reference: "Efficient and accurate rotation of finite spherical harmonics expansions", 
*  Journal of Computational Physics 231 (2012) 243–250 
*/
template<class T> 
NDArray<T,1>
GetRotatedSHCoefficient(const NDArray<T,1>& shInput, const NDArray<T,2>& rotationMatrix);

/**
 * \brief  get the SH coefficients from the symmetric tensor with eigenvalues (e1,e2,e2), e1>e2, and (theta,phi) is the angular direction of the e1 axis.
 */
template < class T >
std::vector<T>
GetSymmetricTensorSHCoef(const T b, const T e1, const T e2, const int lMax, const T theta=0, const T phi=0 );

/**
 * \brief  get the derivatives of the SH coefficients with respect to (e1, e2) 
 * in the symmetric tensor with eigenvalues (e1,e2,e2), e1>e2, and (theta,phi) is the angular direction of the e1 axis.
 * The returned matrix is a CImg<T>(utl::rank2dimSH(lMax),2)
 */
template < class T >
std::vector< std::vector<T> >
GetSymmetricTensorSHCoefDerivative(const T b, const T e1, const T e2, const int lMax, const T theta=0, const T phi=0 );


/** get the legendre coefficient vector of \f$ \exp(-a\times x^2) \f$, 
 * i.e.\f$ \exp(-a\times x^2) = \sum_{l=0}^{lMax} A_l(a) P_l(x) \f$  */
inline double
GetExpLegendreCoef(const double a, const int l );

/** get the derivative \f$\frac{\partial A_l(a)}{\partial a}\f$ of the legendre coefficient vector of \f$ \exp(-a\times x^2) \f$, 
 * i.e.\f$ \exp(-a\times x^2) = \sum_{l=0}^{lMax} A_l(a) P_l(x) \f$  */
inline double
GetExpLegendreCoefDerivative(const double a, const int l );

}


#if !defined(__itkSpecialFunctionGenerator_hxx)
#include "itkSpecialFunctionGenerator.hxx"
#endif 

#endif 
