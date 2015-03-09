/**
 *       @file  itkSpecialFunctionGenerator.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "12-29-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSpecialFunctionGenerator_hxx
#define __itkSpecialFunctionGenerator_hxx


#include "itkSpecialFunctionGenerator.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include "utl.h"
#include <tr1/cmath>

#include "itkSpecialFunctionGenerator.h"

namespace utl
{

double
Hyperg1F1(double a, double b, double x) 
{
  // Handle possible underflows with a Kummer transformation
  // if(x>=-500.0) 
  //   return gsl_sf_hyperg_1F1(a,b,x);
  // else 
  //   return std::exp(x)*gsl_sf_hyperg_1F1(b-a,b,-x);
  return std::tr1::conf_hyperg(a,b,x);
}

template < class T >
T
Lagurre (const int n, const double a, const T x )
{
  if (n==1)
    return gsl_sf_laguerre_1 (a,x); 
  else if (n==2)
    return gsl_sf_laguerre_2 (a,x); 
  else if (n==3)
    return gsl_sf_laguerre_3 (a,x); 
  else
    return gsl_sf_laguerre_n(n,a,x);
}    

double 
Gamma(const double x)
{
  utlException(std::fabs(x)<1e-8, "wrong input x, x=0");

  if (utl::IsInt(2*x) && x>0)
    return GammaHalfInteger(x);
  else
    return gsl_sf_gamma(x);
}

double 
GammaLower ( const double s, const double x )
{
  double gamma_whole = Gamma(s);
  double gamma_upper = gsl_sf_gamma_inc(s,x);
  return gamma_whole - gamma_upper;
}

double
BesselJa(const double a, const double x)
{
  int int_a = int(a);
  int int_2a = int(2.0*a);
  if (std::abs(int_a-a)<1e-8) // integer
    {
    if (a<0)
      return (int_a%2?(-1):1)*BesselJa(-a,x);
    if (int_a==0)
      return gsl_sf_bessel_J0(x);
    else if (int_a==1)
      return gsl_sf_bessel_J1(x);
    else 
      return gsl_sf_bessel_Jn(int_a, x);
    }
  else if (std::abs(int_2a-2.0*a)<1e-8) // half of an intger
    {
    utlException (a<0, "the parameter a should be positive in this routine. \
      Note: a can be negative in theory... \n\
      When a is an integer, J_a(x)=(-1)^a J_{-a}(x); When a is not an integer, J_a(x) and J_{-a}(x) are linearly independent.");
    // j_l(x) = \sqrt(\pi / (2x)) J_{l+0.5}(x)
    return std::sqrt(2.0*x/M_PI) * gsl_sf_bessel_jl(int(a-0.5),x);
    }
  else  // use std::tr1::cyl_bessel_j  for non-negative order.
    utlException(true, "a should be an integer or a half of an integer!");
  return 0;
}


template<class T> 
NDArray<T,1>
GetRotatedSHCoefficients(const NDArray<T,1>& shInput, const NDArray<T,2>& rotationMatrix)
{
  typedef NDArray<double,2> MatrixDouble;
  typedef NDArray<double,1> VectorDouble;
  static utl_shared_ptr<MatrixDouble > gradt3 = utl::ReadGrad<double>(3, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  // static int rank = utl::DimToRankSH(shInput.size());
  static int rank = 10;
  static utl_shared_ptr< MatrixDouble> shMatrix = utl::ComputeSHMatrix(rank, *gradt3, CARTESIAN_TO_SPHERICAL);
  static std::vector< MatrixDouble > shMatrixInv(rank/2);
  static bool isFirstTime = true; 
  if (isFirstTime)
    {
    for ( int l = 2; l <= rank; l += 2 ) 
      {
      MatrixDouble shMatrixL;
      shMatrix->GetNColumns(utl::RankToDimSH(l-2), 2*l+1, shMatrixL);
      utl::PInverseMatrix(shMatrixL, shMatrixInv[l/2-1]);
      }
    isFirstTime = false;
    }

  int rankReal = utl::DimToRankSH(shInput.Size());
  if (rankReal>rank)
    {
    rank = rankReal;
    shMatrix = utl::ComputeSHMatrix(rank, *gradt3, CARTESIAN_TO_SPHERICAL);
    shMatrixInv = std::vector<MatrixDouble >(rank/2);
    for ( int l = 2; l <= rank; l += 2 )
      {
      MatrixDouble shMatrixL;
      shMatrix->GetNColumns(utl::RankToDimSH(l-2), 2*l+1, shMatrixL);
      utl::PInverseMatrix(shMatrixL, shMatrixInv[l/2-1]);
      }
    }

  // MatrixDouble gradt3Rotated = (rotationMatrix.GetTranspose() * gradt3->GetTranspose()).GetTranspose();
  MatrixDouble gradt3Rotated =  (*gradt3) * rotationMatrix;
  utl_shared_ptr<MatrixDouble > shMatrixRotated = utl::ComputeSHMatrix(rankReal, gradt3Rotated,CARTESIAN_TO_SPHERICAL);
  VectorDouble shRotated(shInput);
  MatrixDouble tmp;
  for ( int l = 2; l <= rankReal; l += 2 ) 
    {
    int colstart = utl::RankToDimSH(l-2);
    shMatrixRotated->GetNColumns(colstart, 2*l+1, tmp);
    VectorDouble sfValueRotatedL =  tmp* shInput.GetSubVector(utl::GetRange(colstart, colstart+2*l+1));
    VectorDouble shRotatedL = shMatrixInv[l/2-1]*sfValueRotatedL;
    shRotated.SetSubVector(utl::GetRange(colstart, colstart+shRotatedL.Size()), shRotatedL);
    }
  return shRotated;
}

template < class T >
std::vector<T>
GetSymmetricTensorSHCoef ( const T b, const T e1, const T e2, const int lMax, const T theta, const T phi )
{
  utlException(!utl::IsInt(0.5*lMax), "lMax should be even, lMax="<<lMax);
  utlException(e1<e2-1e-10, "e1 should be more than e2, e1="<<e1 << ", e2="<<e2);

  std::vector<T> coef_vec(utl::RankToDimSH(lMax),0.0);
  double a = (e1 - e2)*b; 
  double expbe2 = std::exp(-b*e2);

  for ( int l = 0; l <= lMax; l += 2 ) 
    {
    double A = GetExpLegendreCoef(a,l);
    for ( int m = -l; m <= l; m += 1 ) 
      {
      int jj = utl::GetIndexSHj(l,m);
      coef_vec[jj] = 4.0*M_PI/(2.0*l+1.0) * expbe2 * A * itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,theta, phi);
      }
    }
  return coef_vec;
}   // -----  end of method SpecialFunctions<T>::<method>  -----

template < class T >
std::vector< std::vector<T> >
GetSymmetricTensorSHCoefDerivative ( const T b, const T e1, const T e2, const int lMax, const T theta, const T phi )
{
  utlException(!utl::IsInt(0.5*lMax), "lMax should be even, lMax="<<lMax);
  utlException(e1<e2-1e-10, "e1 should be more than e2, e1="<<e1 << ", e2="<<e2);

  std::vector<std::vector<T> > coef(2);
  coef[0] = std::vector<T>(utl::RankToDimSH(lMax),0.0);
  coef[1] = std::vector<T>(utl::RankToDimSH(lMax),0.0);

  double a = (e1 - e2)*b; 
  double expbe2 = std::exp(-b*e2);

  for ( int l = 0; l <= lMax; l += 2 ) 
    {
    double A = GetExpLegendreCoef(a,l);
    double dA = GetExpLegendreCoefDerivative(a,l);
    for ( int m = -l; m <= l; m += 1 ) 
      {
      int jj = utl::GetIndexSHj(l,m);
      coef[0][jj] = 4.0*M_PI/(2.0*l+1.0) * expbe2 * b * dA * itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,theta, phi);
      coef[1][jj] = 4.0*M_PI/(2.0*l+1.0) * expbe2 * (-b) * (A + dA);
      }
    }
  return coef;
}   // -----  end of method SpecialFunctions<T>::<method>  -----



}


#endif 
