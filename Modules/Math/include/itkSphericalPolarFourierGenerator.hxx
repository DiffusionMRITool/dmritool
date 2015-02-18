/**
 *       @file  itkSphericalPolarFourierGenerator.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-09-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSphericalPolarFourierGenerator_hxx
#define __itkSphericalPolarFourierGenerator_hxx


#include "itkSphericalPolarFourierGenerator.h"

#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_hyperg.h>
#include "utl.h"
#include "itkSphericalHarmonicsGenerator.h"
#include "itkSpecialFunctionGenerator.h"


namespace itk 
{

template < class PreciseType >
SphericalPolarFourierRadialGenerator<PreciseType>
::SphericalPolarFourierRadialGenerator ()
{
  m_N = 0; m_L=0;
  m_SPFType = SPF;
  m_Scale = -1; 
}  // -----  end of constructor of template class SphericalPolarFourierRadialGenerator  -----


template < class PreciseType >
PreciseType
SphericalPolarFourierRadialGenerator<PreciseType>
::GetNormalizeFacotr ( const bool isFourier) const
{
  utlShowPosition(this->GetDebug());
  utlException(m_Scale<=0, "need to set a scale");

  if (m_SPFType==SPF)
    {
    if (!isFourier)
      {
      PreciseType first_temp = 2.0 / utl::PowHalfInteger((double)this->m_Scale,(double)1.5) * PreciseType(utl::Factorial(m_N)) / utl::GammaHalfInteger(m_N+1.5);
      first_temp = std::sqrt(first_temp);
      return first_temp;
      }
    else
      {
      utlException(true, "no normalizeFacotr in dSPF basis");
      }
    }
  else if (m_SPFType==DSPF)
    {
    if (!isFourier)
      utlException(true, "no normalizeFacotr in dSPF basis");
    else
      {
      Self radial_dual = SphericalPolarFourierRadialGenerator<PreciseType>::New(); 
      radial_dual->SetNLSPF(m_N,m_L,SPF);
      radial_dual->SetScale(m_Scale);
      return radial_dual->GetNormalizeFacotr(false);
      }
    }
  else if (m_SPFType==SHORE)
    {
    if (!isFourier)
      {
      PreciseType first_temp = 2.0 / utl::PowHalfInteger((double)this->m_Scale,1.5) * utl::Factorial(m_N-m_L/2) / utl::GammaHalfInteger(m_N+m_L/2+1.5);
      first_temp = std::sqrt(first_temp);
      return first_temp;
      }
    else
      {
      Self radial_dual = SphericalPolarFourierRadialGenerator<PreciseType>::New(); 
      radial_dual->SetNLSPF(m_N,m_L,SHORE);
      radial_dual->SetScale( 1.0/(4*M_PI*M_PI*this->GetScale()) );
      PreciseType sign = (m_N)%2==0 ? 1 : -1;
      return sign*radial_dual->GetNormalizeFacotr(false);
      }
    }
  else
    utlException(true, "wrong type");
}   // -----  end of method SphericalPolarFourierGenerator<PreciseType>::Evaluate  -----

template < class PreciseType >
PreciseType
SphericalPolarFourierRadialGenerator<PreciseType>
::Evaluate ( const PreciseType qrValue, const bool isFourier) const
{
  utlShowPosition(this->GetDebug());
  utlException(m_Scale<=0, "need to set a scale");

  if (m_SPFType==SPF)
    {
    if (!isFourier)
      {
      PreciseType first_temp = 2.0 / utl::PowHalfInteger((double)this->m_Scale,(double)1.5) * PreciseType(utl::Factorial(m_N)) / utl::GammaHalfInteger(m_N+1.5);
      first_temp = std::sqrt(first_temp);
      if (this->GetDebug())
        {
        std::cout << "radial order = " << m_N << std::endl;
        std::cout << "utl::Factorial("<<m_N<<") = " << utl::Factorial(m_N) << std::endl;
        std::cout << "std::pow(m_Scale,1.5) = " << std::pow((double)this->m_Scale,(double)1.5) << std::endl;
        std::cout << "utl::Gamma(m_N+1.5) = " << utl::GammaHalfInteger(m_N+1.5) << std::endl;
        std::cout << "first_temp = " << first_temp << std::endl;
        }
      PreciseType result = first_temp * std::exp( -1.0*qrValue*qrValue / (2*this->m_Scale) ) 
        * utl::Lagurre(m_N, 0.5, qrValue*qrValue/this->m_Scale);
      return result;
      }
    else
      {
      Pointer radial_dual = SphericalPolarFourierRadialGenerator<PreciseType>::New(); 
      radial_dual->SetNLSPF(m_N,m_L,DSPF);
      radial_dual->SetScale(m_Scale);
      return radial_dual->Evaluate(qrValue, false);
      }
    }
  else if (m_SPFType==DSPF)
    {
    if (!isFourier)
      {
      PreciseType first_temp = 2.0 / utl::PowHalfInteger((double)this->m_Scale,1.5) * utl::Factorial(m_N) / utl::GammaHalfInteger(m_N+1.5);
      first_temp = std::sqrt(first_temp) * utl::PowHalfInteger((double)this->m_Scale, 1.5);
      std::vector<PreciseType> coef_laguerre = utl::GetCoefLaguerre(m_N);
      if (this->GetDebug())
        {
        std::cout << "radial order = " << m_N << std::endl;
        std::cout << "utl::Factorial("<<m_N<<") = " << utl::Factorial(m_N) << std::endl;
        std::cout << "std::pow(m_Scale,1.5) = " << std::pow((double)this->m_Scale,1.5) << std::endl;
        std::cout << "utl::Gamma(m_N+1.5) = " << utl::GammaHalfInteger(m_N+1.5) << std::endl;
        std::cout << "first_temp = " << first_temp << std::endl;
        }
      PreciseType rr = qrValue*std::sqrt(this->m_Scale);
      if (this->GetDebug())
        std::cout << "qrValue = " << qrValue << ", relative qrvalue = " << rr << std::endl;
      PreciseType sign = 4.0;
      if ( (m_L/2) % 2 == 1)
        sign = -4.0;
      sign *= utl::PowHalfInteger(M_PI, m_L+1.5)*utl::PowHalfInteger(rr,m_L) / utl::GammaHalfInteger(m_L+1.5);
      PreciseType integral = 0;
      for ( int i = 0; i <= m_N; i += 1 ) 
        {
        PreciseType first_C = utl::PowHalfInteger(2,0.5*m_L+i-0.5) *  utl::GammaHalfInteger(i+0.5*m_L+1.5);
        // std::cout << "first_C = " << first_C << std::endl;
        double a = double(i+0.5*m_L+1.5); 
        double b = double(m_L+1.5);
        double c = double(-2*M_PI*M_PI*rr*rr);
        // std::cout << "(a,b,c)=("<<a<<","<<b<<","<<c<<")" << std::endl;
        // integral += coef_laguerre[i]* first_C; 
        // NOTE: when c is small, 1F1=1 
        if (-1.0*c>1e-7)
          integral += coef_laguerre[i]* first_C *gsl_sf_hyperg_1F1(a, b, c);
        else
          integral += coef_laguerre[i]* first_C;
        // std::cout << "integral = " << integral << std::endl;
        }
      PreciseType result = sign * first_temp * integral;
      return result;
      }
    else
      {
      Pointer radial_dual = SphericalPolarFourierRadialGenerator<PreciseType>::New(); 
      radial_dual->SetNLSPF(m_N,m_L,SPF);
      radial_dual->SetScale(m_Scale);
      return radial_dual->Evaluate(qrValue, false);
      }
    }
  else if (m_SPFType==SHORE)
    {
    if (!isFourier)
      {
      PreciseType first_temp = 2.0 / utl::PowHalfInteger((double)this->m_Scale,1.5) * utl::Factorial(m_N-m_L/2) / utl::GammaHalfInteger(m_N+m_L/2+1.5);
      first_temp = std::sqrt(first_temp);
      PreciseType x2_temp = qrValue*qrValue / this->m_Scale;
      PreciseType result = first_temp * std::exp(-0.5*x2_temp)* utl::PowHalfInteger((double)x2_temp, m_L/2) * utl::Lagurre(m_N-m_L/2, m_L+0.5, x2_temp);
      if (this->GetDebug())
        {
        utlPrintVar3(true, m_N, m_L, x2_temp);
        PreciseType exp_temp = std::exp( -1.0*x2_temp / 2 );
        PreciseType lagurre_temp = utl::Lagurre(m_N-m_L/2, m_L+0.5, x2_temp);
        PreciseType pow_temp = utl::PowHalfInteger((double)x2_temp, m_L/2);
        utlPrintVar4(this->GetDebug(),m_N,m_L, std::sqrt(x2_temp), qrValue);
        utlPrintVar5(this->GetDebug(), first_temp, exp_temp, lagurre_temp, pow_temp, result);
        }
      return result;
      }
    else
      {
      Pointer radial_dual = SphericalPolarFourierRadialGenerator<PreciseType>::New(); 
      radial_dual->SetNLSPF(m_N,m_L,SHORE);
      radial_dual->SetScale( 1.0/(4*M_PI*M_PI*m_Scale) );
      PreciseType sign = (m_N)%2==0 ? 1 : -1;
      return sign*radial_dual->Evaluate(qrValue, false);
      }
    }
  else
    utlException(true, "wrong type");
  return 0;
}   // -----  end of method SphericalPolarFourierGenerator<PreciseType>::Evaluate  -----

template < class PreciseType >
PreciseType
SphericalPolarFourierGenerator<PreciseType>
::Evaluate (const PreciseType qrValue, const PreciseType theta, const PreciseType phi, const bool isFourier) const
{
  utlShowPosition(this->GetDebug());
  
  PreciseType radial_value = m_Radial->Evaluate(qrValue, isFourier);
  PreciseType spherical_value = SphericalHarmonicsGenerator<PreciseType>::RealSH(m_Radial->GetL(),m_M,theta,phi);
  utlPrintVar2(this->GetDebug(), radial_value, spherical_value);
  return radial_value*spherical_value;
}   // -----  end of method SphericalPolarFourierGenerator<PreciseType>::Evaluate  -----

template < class PreciseType >
PreciseType
SphericalPolarFourierGenerator<PreciseType>
::EvaluateRadial (const PreciseType qrValue, const bool isFourier) const
{
  return m_Radial->Evaluate(qrValue, isFourier);
}

template <class PreciseType>
void 
SphericalPolarFourierGenerator<PreciseType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  // PrintVar4(true, GetN(), GetL(), GetM(), GetScale(), os<<indent);
  typename RadialType::SPFType spfType = this->GetSPFType();
  if (spfType==SPF)
    os << indent << "SPFType: SPF basis" << std::endl;
  else if (spfType==DSPF)
    os << indent << "SPFType: Dual SPF basis" << std::endl;
  else if (spfType==SHORE)
    os << indent << "SPFType: SHORE basis" << std::endl;
  else 
    utlException(true, "wrong SPF type");
}


}



#endif 
