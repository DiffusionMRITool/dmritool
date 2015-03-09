/**
 *       @file  itkSphericalHarmonicsGenerator.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  10-29-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */  
#ifndef __itkSphericalHarmonicsGenerator_hxx
#define __itkSphericalHarmonicsGenerator_hxx

#include "itkSphericalHarmonicsGenerator.h"
#include <complex>
#include <gsl/gsl_sf_legendre.h>
#include "itkImage.h"
#include <gsl/gsl_sf_coupling.h>

// #include "utl.h"
#include "DMRITOOLConfigure.h"
#include "utlCore.h"
#include "utlDMRI.h"
#include "utlITK.h"

namespace utl 
{
double Gamma(const double);
}

namespace itk 
{

template<class PreciseType>
SphericalHarmonicsGenerator<PreciseType>
::SphericalHarmonicsGenerator() 
{ 
}

template <class PreciseType>
SphericalHarmonicsGenerator<PreciseType>
::~SphericalHarmonicsGenerator() 
{
}

template < class PreciseType >
void
SphericalHarmonicsGenerator<PreciseType>
::PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

template<class PreciseType> 
std::complex<PreciseType> 
SphericalHarmonicsGenerator<PreciseType>
::ComplexSH(const int l, const int m, const PreciseType theta, const PreciseType phi) 
{
  int absm = std::abs(int(m));
  PreciseType sign = utl::IsEven(absm)? 1.0 : -1.0;

  std::complex<PreciseType> retval(0.0,(PreciseType)(absm*phi));
  retval = std::exp(retval) * gsl_sf_legendre_sphPlm(l,absm,cos(theta));  

  if (m<0) 
    retval = sign * std::conj(retval); 

  return retval;
}



template < class PreciseType >
PreciseType
SphericalHarmonicsGenerator<PreciseType>
::RealSH (const int l, const int m, const PreciseType theta, const PreciseType phi )
{
  std::complex<PreciseType>  cplx;
  if(m > 0) 
    {
    cplx = ComplexSH(l, m, theta, phi);
    return std::sqrt(2.0)*std::imag(cplx);
    }
  else if( m == 0 )
    {
    cplx = ComplexSH(l, m, theta, phi);
    return std::real(cplx);
    }
  else  
    {
    // NOTE: use m%2==1 is wrong, because it can be -1.
    PreciseType sign = utl::IsEven(m)? 1.0 : -1.0;
    // utlPrintVar5(true, l, m, utl::IsEven(m), m%2, -m%2);
    cplx = sign * ComplexSH(l, m, theta, phi);
    return std::sqrt(2.0)*std::real(cplx);
    }
}

template < class PreciseType >
PreciseType
SphericalHarmonicsGenerator<PreciseType>
::ComplexTripleIntegration (const int l1, const int m1, const int l2, const int m2, const int l3, const int m3)
{
  utlException(l1<std::abs(m1) || l2<std::abs(m2) || l3<std::abs(m3), "wrong input! m should be less than l");
  PreciseType result = 0.0;
  if (m1+m2+m3==0 && l1+l2>=l3 && l1+l3>=l2 && l2+l3>=l1)
    {
    result = std::sqrt( (2*l1+1)*(2*l2+1)*(2*l3+1)/(4*M_PI) ) * gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3) * gsl_sf_coupling_3j(2*l1,2*l2,2*l3,0,0,0);
    }
  else
    result = 0;
  return result;
}

template < class PreciseType >
PreciseType
SphericalHarmonicsGenerator<PreciseType>
::RealTripleIntegration (const int l1, const int m1, const int l2, const int m2, const int l3, const int m3, const bool is_precalculated)
{
  utlException(l1<std::abs(m1) || l2<std::abs(m2) || l3<std::abs(m3), "wrong input! m should be less than l");
  utlException((l1%2!=0) || (l2%2!=0) || (l3%2!=0), "l1, l2, l3 need to be even integers" );
  PreciseType result = 0;
  // Sqrt[2]*Re[SphericalHarmonicY[l, Abs[m], \[Theta], \[Phi]]] = 1/Sqrt[2]*(SphericalHarmonicY[l,-m,Theta,Phi]+SphericalHarmonicY[l,m,Theta,Phi])
  if (is_precalculated)
    {
    utlException(!utl::IsFileExist(utl::SH3Itegralhdr), "no SH3Itegralhdr.");
    typedef itk::Image<double,3> ImageType;
    static bool isFirstTime = true;
    static ImageType::Pointer image = ImageType::New();
    static int rank = -1;
    ImageType::IndexType pixelIndex;
    if (isFirstTime)
      {
      itk::ReadImage<ImageType>(utl::SH3Itegralhdr, image);
      ImageType::RegionType region = image->GetLargestPossibleRegion();
      ImageType::SizeType size = region.GetSize();
      rank = utl::DimToRankSH(size[0]);
      isFirstTime = false;
      }
    if (l1<=rank && l2<=rank && l3<=rank)
      {
      pixelIndex[0] = utl::GetIndexSHj(l1,m1);
      pixelIndex[1] = utl::GetIndexSHj(l2,m2);
      pixelIndex[2] = utl::GetIndexSHj(l3,m3);
      result = image->GetPixel(pixelIndex);
      }
    else
      result = RealTripleIntegration(l1,m1,l2,m2,l3,m3,false);
    }
  else
    {
    // NOTE: order 0 -2 -1 0 1 2 -4 -3 -2 -1 0 1 2 3 4 ...
    // if m<0, Y[l,m] = Sqrt[2] * Re[SHY[l,Abs[m]]]
    // if m=0, Y[l,m] = Re[SHY[l,m]], (imaginary part==0)
    // if m>0, Y[l,m] = Sqrt[2] * Im[SHY[l,m]]
    if (m1<0 && m2<0 && m3<0) // result = Sqrt[2]/4* Integrate[ ((-1)^(-m1)*Y[l1,m1]+Y[l1,-m1]) * ((-1)^(-m2)*Y[l2,m2]+Y[l2,-m2]) * ((-1)^(-m3)*Y[l3,m3]+Y[l3,-m3]) ]
      {
      if (m1+m2==m3)
        result = ( m3%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,m2,l3,-m3) / std::sqrt(2); // (-1)^(m1)*(-1)^(m2)*Integrate[Y[l1,m1]*Y[l2,m2]*Y[l3,-m3]]==(-1)^m3*Integrate[Y[l1,-m1]*Y[l2,-m2]*Y[l3,m3]]
      else if (m1+m3==m2)
        result = ( m2%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,-m2,l3,m3) / std::sqrt(2); 
      else if (m2+m3==m1)
        {
        // utlPrintVar3(true,l1,l2,l3);
        // utlPrintVar3(true,m1,m2,m3);
        // utlPrintVar3(true,m1%2, ( m1%2==1?-1:1 ),ComplexTripleIntegration(l1,-m1,l2,m2,l3,m3));
        // utlPrintVar1(true, ComplexTripleIntegration(l1,-m1,l2,m2,l3,m3));
        result = ( m1%2==0?1:-1 ) * ComplexTripleIntegration(l1,-m1,l2,m2,l3,m3) / std::sqrt(2); 
        }
      }
    else if ( (m1>0 && m2<0 && m3<0) || (m1<0 && m2>0 && m3<0) || (m1<0 && m2<0 && m3>0) )  // Integrate[Sin[m1*Phi]*Cos[m2*Phi]*Cos[m3*Phi]]=0
      {
      result = 0;
      }
    else if (m1>0 && m2>0 && m3<0) // result = -Sqrt[2]/4* Integrate[ (Y[l1,m1]-(-1)^(m1)*Y[l1,-m1]) * (Y[l2,m2]-(-1)^(-m2)*Y[l2,-m2]) * ((-1)^(-m3)*Y[l3,m3]+Y[l3,-m3]) ]
      {
      if (m1+m2+m3==0)
        result = (-1.0) * ( m3%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,m2,l3,m3) / std::sqrt(2); // (-1)^(m1)*(-1)^(m2)*Integrate[Y[l1,m1]*Y[l2,m2]*Y[l3,-m3]]==(-1)^m3*Integrate[Y[l1,-m1]*Y[l2,-m2]*Y[l3,m3]]
      else if (m1==m2+m3)
        result = ( m2%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,-m2,l3,-m3) / std::sqrt(2); 
      else if (m2==m1+m3)
        result = ( m1%2==0?1:-1 ) * ComplexTripleIntegration(l1,-m1,l2,m2,l3,-m3) / std::sqrt(2); 
      }
    else if (m1>0 && m2<0 && m3>0) 
      {
      if (m1+m2+m3==0)
        result = (-1.0) * ( m2%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,m2,l3,m3) / std::sqrt(2); // (-1)^(m1)*(-1)^(m2)*Integrate[Y[l1,m1]*Y[l2,m2]*Y[l3,-m3]]==(-1)^m3*Integrate[Y[l1,-m1]*Y[l2,-m2]*Y[l3,m3]]
      else if (m1==m2+m3)
        result = ( m3%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,-m2,l3,-m3) / std::sqrt(2); 
      else if (m3==m1+m2)
        result = ( m1%2==0?1:-1 ) * ComplexTripleIntegration(l1,-m1,l2,-m2,l3,m3) / std::sqrt(2); 
      }
    else if (m1<0 && m2>0 && m3>0) 
      {
      if (m1+m2+m3==0)
        result = (-1.0) * ( m1%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,m2,l3,m3) / std::sqrt(2); // (-1)^(m1)*(-1)^(m2)*Integrate[Y[l1,m1]*Y[l2,m2]*Y[l3,-m3]]==(-1)^m3*Integrate[Y[l1,-m1]*Y[l2,-m2]*Y[l3,m3]]
      else if (m3==m2+m1)
        result = ( m2%2==0?1:-1 ) * ComplexTripleIntegration(l1,-m1,l2,-m2,l3,m3) / std::sqrt(2); 
      else if (m2==m1+m3)
        result = ( m3%2==0?1:-1 ) * ComplexTripleIntegration(l1,-m1,l2,m2,l3,-m3) / std::sqrt(2); 
      }
    else if (m1>0 && m2>0 && m3>0) // Integrate[Sin[m1*Phi]*Sin[m2*Phi]*Sin[m3*Phi]]=0
      result = 0;
    else if (m1==0 && m2==0 && m3==0)
      result = ComplexTripleIntegration(l1,0,l2,0,l3,0); 
    else if (m1==0 && m2==m3)
      result = ( m3%2==0?1:-1 ) * ComplexTripleIntegration(l1,0,l2,m2,l3,-m3);
    else if (m2==0 && m3==m1)
      result = ( m1%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,0,l3,-m3);
    else if (m3==0 && m1==m2)
      result = ( m2%2==0?1:-1 ) * ComplexTripleIntegration(l1,m1,l2,-m2,l3,0);
    }
  return result;
}

template < class PreciseType >
std::complex<double>
SphericalHarmonicsGenerator<PreciseType>
::ComplexDerivativeOfTheta ( const int l, const int m, const double theta, const double phi)
{
  std::complex<double> result(0,0), ee(0,-phi);
  std::complex<double> sh_value = ComplexSH(l, m, theta, phi);
  double theta_used = std::abs(theta)>1e-6? theta : (theta>=0?1e-6:-1e-6);
  if (m==l)
    return l * 1.0/std::tan(theta_used) * sh_value;
  std::complex<double> sh_value_2 = ComplexSH(l, m+1, theta_used, phi);
  result = m * 1.0/std::tan(theta_used) * sh_value + std::sqrt( utl::Gamma(1+l-m) * utl::Gamma(2+l+m) / (utl::Gamma(l-m)*utl::Gamma(l+m+1)) ) * sh_value_2 * ee;
  return result;
}

template < class PreciseType >
std::complex<double>
SphericalHarmonicsGenerator<PreciseType>
::ComplexDerivativeOfPhi ( const int l, const int m, const double theta, const double phi)
{
  std::complex<double> sh_value = ComplexSH(l, m, theta, phi);
  std::complex<double> r1 (0,m);
  return r1*sh_value;
}

template < class PreciseType >
double
SphericalHarmonicsGenerator<PreciseType>
::RealDerivativeOfTheta ( const int l, const int m, const double theta, const double phi)
{
  if (m==0)
    return std::real(ComplexDerivativeOfTheta(l,m,theta,phi));
  else if (m>0)
    return std::sqrt(2.0) * std::imag(ComplexDerivativeOfTheta(l,m,theta,phi));
  else 
    return std::sqrt(2.0) * std::real(ComplexDerivativeOfTheta(l,-m,theta,phi));
}

template < class PreciseType >
double 
SphericalHarmonicsGenerator<PreciseType>
::RealDerivativeOfPhi ( const int l, const int m, const double theta, const double phi)
{
  if (m==0)
    return std::real(ComplexDerivativeOfPhi(l,m,theta,phi));
  else if (m>0)
    return std::sqrt(2.0) * std::imag(ComplexDerivativeOfPhi(l,m,theta,phi));
  else 
    return std::sqrt(2.0) * std::real(ComplexDerivativeOfPhi(l,-m,theta,phi));
}

}

#endif 
