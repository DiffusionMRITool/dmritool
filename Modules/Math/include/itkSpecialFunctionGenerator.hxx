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
#include "itkUnaryFunctorLookUpTable.h"

#include "utlRotationMatrixFromVectors.h"
#include "itkSHCoefficientsRotation.h"

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
    return BesselJInteger(int_a, x);
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

double
BesselJInteger(const int n, const double x)
{
  if (n<0)
    return (n%2?(-1):1)*BesselJInteger(-n,x);
  if (n==0)
    return gsl_sf_bessel_J0(x);
  else if (n==1)
    return gsl_sf_bessel_J1(x);
  else 
    return gsl_sf_bessel_Jn(n, x);
}

double
BesselJIntegerPrime(const int n, const double x)
{
  if (n==0)
    return -BesselJInteger(1,x);
  else
    return 0.5*(BesselJInteger(n-1,x) - BesselJInteger(n+1,x));
}

// template<class T> 
// NDArray<T,1>
// GetRotatedSHCoefficients(const NDArray<T,1>& shInput, const NDArray<T,2>& rotationMatrix)
// {
//   typedef NDArray<double,2> MatrixDouble;
//   typedef NDArray<double,1> VectorDouble;
//   static utl_shared_ptr<MatrixDouble > gradt3 = utl::ReadGrad<double>(3, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
//   // static int rank = utl::DimToRankSH(shInput.size());
//   static int rank = 10;
//   static utl_shared_ptr< MatrixDouble> shMatrix = utl::ComputeSHMatrix(rank, *gradt3, CARTESIAN_TO_SPHERICAL);
//   static std::vector< MatrixDouble > shMatrixInv(rank/2);
//   static bool isFirstTime = true; 
//   if (isFirstTime)
//     {
//     for ( int l = 2; l <= rank; l += 2 ) 
//       {
//       MatrixDouble shMatrixL = shMatrix->GetNColumns(utl::RankToDimSH(l-2), 2*l+1);
//       utl::PInverseMatrix(shMatrixL, shMatrixInv[l/2-1]);
//       }
//     isFirstTime = false;
//     }

//   int rankReal = utl::DimToRankSH(shInput.Size());
//   if (rankReal>rank)
//     {
//     rank = rankReal;
//     shMatrix = utl::ComputeSHMatrix(rank, *gradt3, CARTESIAN_TO_SPHERICAL);
//     shMatrixInv = std::vector<MatrixDouble >(rank/2);
//     for ( int l = 2; l <= rank; l += 2 )
//       {
//       MatrixDouble shMatrixL = shMatrix->GetNColumns(utl::RankToDimSH(l-2), 2*l+1);
//       utl::PInverseMatrix(shMatrixL, shMatrixInv[l/2-1]);
//       }
//     }

//   // MatrixDouble gradt3Rotated = (rotationMatrix.GetTranspose() * gradt3->GetTranspose()).GetTranspose();
//   MatrixDouble gradt3Rotated =  (*gradt3) * rotationMatrix;
//   utl_shared_ptr<MatrixDouble > shMatrixRotated = utl::ComputeSHMatrix(rankReal, gradt3Rotated,CARTESIAN_TO_SPHERICAL);
//   VectorDouble shRotated(shInput);
//   MatrixDouble tmp;
//   for ( int l = 2; l <= rankReal; l += 2 ) 
//     {
//     int colstart = utl::RankToDimSH(l-2);
//     tmp = shMatrixRotated->GetNColumns(colstart, 2*l+1);
//     VectorDouble sfValueRotatedL =  tmp* shInput.GetSubVector(utl::GetRange(colstart, colstart+2*l+1));
//     VectorDouble shRotatedL = shMatrixInv[l/2-1]*sfValueRotatedL;
//     shRotated.SetSubVector(utl::GetRange(colstart, colstart+shRotatedL.Size()), shRotatedL);
//     }
//   return shRotated;
// }

template < class T >
utl_shared_ptr<utl::NDArray<T,1> >
GetSymmetricTensorSHCoef ( const T b, const T e1, const T e2, const int lMax, const T theta, const T phi )
{
  utlException(!utl::IsInt(0.5*lMax), "lMax should be even, lMax="<<lMax);
  utlException(e1<e2-1e-10, "e1 should be more than e2, e1="<<e1 << ", e2="<<e2);
  
  typedef utl::NDArray<T, 1> VectorType;
  typedef utl_shared_ptr<VectorType>  VectorPointer;

  // utlPrintVar6(true, b, e1, e2, lMax, theta, phi);
  VectorPointer coef_vec (new VectorType(utl::RankToDimSH(lMax),T(0.0)));
  double a = (e1 - e2)*b; 
  double expbe2 = itk::lutExpValue(-b*e2);

  for ( int l = 0; l <= lMax; l += 2 ) 
    {
    double A = utl::GetExpLegendreCoef(a,l, &itk::lutExpValue);
    if (std::fabs(theta)<1e-10)
      {
      int jj = utl::GetIndexSHj(l,0);
      (*coef_vec)[jj] = 4.0*M_PI/(2.0*l+1.0) * expbe2 * A * std::sqrt(1+2*l)/(2*utl::SQRTPI);
      }
    else
      {
      for ( int m = -l; m <= l; m += 1 ) 
        {
        int jj = utl::GetIndexSHj(l,m);
        (*coef_vec)[jj] = 4.0*M_PI/(2.0*l+1.0) * expbe2 * A * itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,theta, phi);
        }
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
  double expbe2 = itk::lutExpValue(-b*e2);

  for ( int l = 0; l <= lMax; l += 2 ) 
    {
    double A = utl::GetExpLegendreCoef(a,l, &itk::lutExpValue);
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


inline double
GetExpProductLegendreCoef(const double a, const double b, const int l )
{
  if (!utl::IsInt(0.5*l))
    return 0;

  if (a==b)
    return l==0? itk::lutExpValue(-b) : 0;
  else
    return itk::lutExpValue(-b)* utl::GetExpLegendreCoef(a-b,l, &itk::lutExpValue);
}

template < class T >
utl_shared_ptr<utl::NDArray<T,1> >
ComputeDWISHCoefficientsForGPDCylinder(const T radius, const T diffusivity, const T deltaBig, const T deltaSmall, const T qq, const int lMax, const T theta, const T phi)
{
  double tau = deltaBig - deltaSmall/3.0; 
  double q_norm = qq*2.0*M_PI;
  double a = tau*q_norm*q_norm*diffusivity;

  typedef utl::NDArray<double, 1> VectorType;
  typedef utl_shared_ptr<VectorType>  VectorPointer;
  const int m_max = 60;
  const VectorType besselJZeros(utl::BesselJPrimeZerosOrder1,m_max);
  VectorType beta_m = besselJZeros / radius; 
  // utl::PrintUtlVector(besselJZeros, "zeros");
  // utl::PrintUtlVector(beta_m, "beta_m");

  double sum_orth=0.0;
  // utlPrintVar4(true, diffusivity, deltaBig, deltaSmall, radius);
  for ( int m = 0; m < m_max; ++m ) 
    {
    double dam = diffusivity*beta_m[m]*beta_m[m];
    double numerator = 2.0*dam*deltaSmall -2 
      + 2.0*itk::lutExpValue(-dam*deltaSmall)
      + 2.0*itk::lutExpValue(-dam*deltaBig)
      - itk::lutExpValue(-dam*(deltaBig - deltaSmall))
      - itk::lutExpValue(-dam*(deltaBig + deltaSmall));
    double beta_m2 = beta_m[m]*beta_m[m];
    double denominator = dam*dam*beta_m2 * (radius*radius*beta_m2 - 1);
    // utlPrintVar2(true, numerator, denominator);
    sum_orth += numerator/denominator; 
    }

  // utlPrintVar3(true, q_norm, deltaSmall, sum_orth);
  double b = 2.0* q_norm*q_norm/(deltaSmall*deltaSmall) * sum_orth;

  VectorPointer coef_vec(new VectorType(utl::RankToDimSH(lMax),0.0));
  for ( int l = 0; l <= lMax ; l+=2 ) 
    {
    double A = utl::GetExpProductLegendreCoef(a, b, l);
    // utlPrintVar4(true, a, b, l, A);
    if (std::fabs(theta)<1e-10)
      {
      int jj = utl::GetIndexSHj(l,0);
      (*coef_vec)[jj] = 4.0*M_PI/(2.0*l+1.0) * A * std::sqrt(1+2*l)/(2*utl::SQRTPI);
      }
    else
      {
      for ( int m = -l; m <= l; m += 1 ) 
        {
        int jj = utl::GetIndexSHj(l,m);
        (*coef_vec)[jj] = 4.0*M_PI/(2.0*l+1.0) * A * itk::SphericalHarmonicsGenerator<double>::RealSH(l,m,theta, phi);
        }
      }
    }

  // utl::PrintUtlVector(*coef_vec,"coef_vec");
  return coef_vec;
}


template < class T >
double
ComputeOrientationalOrderFromSHCoefficients(const utl::NDArray<T,1>& shCoef, const utl::NDArray<T,1>& axis)
{
  utlException(axis.Size()!=3, "wrong size of axis");
  utlException(std::fabs(axis.GetSquaredTwoNorm()-1)>1e-10, "axis should have unit norm");

  const double normZ = 2.0*std::sqrt(M_PI/5.0);
  
  // only l=2 is used. 
  int N = utl::RankToDimSH(2);
  utl::NDArray<T,1> shCoef2(N);
  for ( int i = 0; i < N; ++i ) 
    shCoef2[i] = shCoef[i];

  typedef itk::SHCoefficientsRotation<double> SHRotationFilterType;
  SHRotationFilterType::Pointer shRotateFilter = SHRotationFilterType::GetInstance();  

  // assume rank and Initialize have been set to the global instance.
  utlSAException(shRotateFilter->GetMaxRank()<2)
    (shRotateFilter->GetMaxRank()).msg("need to set large rank");

  utl::NDArray<T,1> zAxis(3,0.0), shCoefRot;
  zAxis[2]=1.0;
  utl::NDArray<T,2> matRot(3,3);
  utl::RotationMatrixFromVectors(axis, zAxis, matRot);

  shCoefRot = shRotateFilter->GetRotatedSHCoefficients(shCoef2, matRot);
  double sh20 = shCoefRot[utl::GetIndexSHj(2,0)];

  return normZ * sh20;
}

double
ComputeOrientationalOrderFromSymmetricTensor(const double e1, const double e2, const double phi)
{
  double diff = e1-e2;
  if (utl::IsSame(diff,0.0, 1e-10))
    return 0;

  double sqrtdiff = std::sqrt(diff);
  double sqrte2 = std::sqrt(e2);
  double oo = ( sqrtdiff*(2.*e1+e2) - 3.*e1*sqrte2*std::atan(sqrtdiff/sqrte2) ) / (2.0* diff*sqrtdiff);

  if (utl::IsSame(phi, 0.0, 1e-10))
    return oo;
  else
    return oo*(1.0+3.*std::cos(2.*phi))/4.0;
}


}


#endif 
