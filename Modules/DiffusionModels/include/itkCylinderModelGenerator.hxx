/**
 *       @file  itkCylinderModelGenerator.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-02-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkCylinderModelGenerator_hxx
#define __itkCylinderModelGenerator_hxx

#include "itkCylinderModelGenerator.h"
#include "itkSpecialFunctionGenerator.h"

namespace itk
{

template <class PreciseType>
CylinderModelGenerator<PreciseType>
::CylinderModelGenerator() : Superclass() 
{
  m_Length = 5.0; 
  m_Radius = 0.005;
  m_D0 = 2.02e-3;
  // m_DeltaBig = 0.0208;
  // m_DeltaSmall = 0.0024;

  m_CylinderAxis[0]=0;
  m_CylinderAxis[1]=0;
  m_CylinderAxis[2]=1.0;
  
  m_LUTExp = LUTExpType::New();
  m_LUTExp->SetVariableMax(0);
  m_LUTExp->SetVariableMin(-30);
  m_LUTExp->SetNumberOfBins(30*1e3);
  // m_LUTExp->BuildLUT();
}

template <class PreciseType>
typename LightObject::Pointer
CylinderModelGenerator<PreciseType>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_Length = m_Length;
  rval->m_Radius = m_Radius;
  rval->m_D0 = m_D0;
  rval->m_CylinderAxis = m_CylinderAxis;
  rval->m_LUTExp = m_LUTExp->Clone();

  return loPtr;
}

template <class PreciseType>
void
CylinderModelGenerator<PreciseType>
::BuildTable()
{
  m_LUTExp->BuildTable();
}

template <class PreciseType>
void
CylinderModelGenerator<PreciseType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar3(true, m_Length, m_Radius, m_D0, os<<indent);
  os << indent << "m_CylinderAxis = " << m_CylinderAxis << std::endl;
}

// template <class PreciseType>
// void
// CylinderModelGenerator<PreciseType>
// ::VerifyInputParameters() const
// {
// }

template <class PreciseType>
void
CylinderModelGenerator<PreciseType>
::Rotate(const MatrixType& mat)
{
  VectorType axis = m_CylinderAxis.GetVnlVector();
  axis = mat*axis;
  for ( int i = 0; i < 3; i += 1 ) 
    m_CylinderAxis[i] = axis[i];
}

template <class PreciseType>
void
CylinderModelGenerator<PreciseType>
::ComputeDWISamples ()
{
  utlShowPosition(this->GetDebug());
  this->VerifyInputParameters();
  utlGlobalException(this->m_SamplingSchemeQSpace->GetNumberOfSamples()==0, "need to set m_SamplingSchemeQSpace");
  utlGlobalException(!m_LUTExp->IsTableBuilt(), "need to build m_LUTExp first");

  const int nMax = 1000, mMax = 10, kMax = 10;
  int N = this->m_SamplingSchemeQSpace->GetNumberOfSamples();
  STDVectorPointer qVector = this->m_SamplingSchemeQSpace->GetRadiusVector();

  double deltaBig = this->m_SamplingSchemeQSpace->GetDeltaBig();
  double d0DeltaBig = m_D0*deltaBig;
  double radius2 = m_Radius*m_Radius;

  std::vector<double> nPiLength2(nMax+1);
  std::vector<double> nPiLength2D0Delta(nMax+1);
  for ( int n = 0; n <= nMax; n += 1 ) 
    {
    nPiLength2[n] = std::pow(n*M_PI/m_Length,2);
    nPiLength2D0Delta[n] = d0DeltaBig *nPiLength2[n];
    }

  this->m_DWISamples = VectorPointer(new VectorType(N));

  // m_LUTExp->Print(std::cout<<"lut = ");


  int i=0;
#pragma omp parallel for private (i) 
  for(i=0; i < N; i++) 
    {

    PointType sample = (*this->m_SamplingSchemeQSpace)[i];
    double dott = sample[0]*m_CylinderAxis[0] + sample[1]*m_CylinderAxis[1] + sample[2]*m_CylinderAxis[2];
    double angle = std::acos( dott>=1?1:(dott<=-1?-1:dott) );
    if (angle>M_PI/2) // only in [0,90]
      angle = M_PI - angle;
    if (M_PI/2-angle<1.0e-5) // NOTE: 0 and 90 is a sigular point
      angle = M_PI/2 - 1.0e-5;
    if (angle<1e-5)
      angle = 1e-5;
    double angle2sin = std::sin(2*angle);

    double signal_temp=0;
    double qq = (*qVector)[i];
    double pi2qcos = 2.0*M_PI*qq*std::cos(angle);
    double pi2qsin = 2.0*M_PI*qq*std::sin(angle);
    double pi2qsinRadius = pi2qsin*m_Radius;
    double pi2qcosLengthcos = std::cos(pi2qcos*m_Length);
    double first_temp = 2.0*std::pow(m_Radius,2) *std::pow(2*M_PI*qq,4)*angle2sin*angle2sin / (m_Length*m_Length);
    MatrixType third_temp(mMax+1, nMax+1);
    for ( int n = 0; n <= nMax; n += 1 ) 
      {
      double num2 = (n%2 ? (1+pi2qcosLengthcos) : (1-pi2qcosLengthcos));
      double tmp1 = std::pow( nPiLength2[n] - std::pow(pi2qcos,2), 2 );
      for ( int m = 0; m <= mMax; m += 1 ) 
        {
        // utlPrintVar2(true, m, n);
        double K_mn = m>0? (n>0?4:2) : (n>0?2:1);
        third_temp(m,n) = num2 *K_mn / tmp1;
        }
      }

    std::vector<double> besselJnVec(mMax+2);
    for ( int m = 0; m <= mMax; m += 1 ) 
      besselJnVec[m] = gsl_sf_bessel_Jn(m,pi2qsinRadius);

    for ( int m = 0; m <=mMax; m += 1 ) 
      {
      double J_m_d = m/(pi2qsinRadius)*besselJnVec[m] - besselJnVec[m+1];
      double J_m_d2 = J_m_d*J_m_d;
      for ( int k = 1; k <=kMax; k += 1 ) 
        {
        double gamma_km = utl::BesselJPrimeZerosTable[kMax*m+k-1];
        double gamma_km2 = gamma_km * gamma_km;
        // utlPrintVar4(true, m,k, gamma_km, J_m_d2);
        double gamma_km_Radius2_D0Delta = gamma_km2/radius2 * d0DeltaBig;
        double second_temp = J_m_d2 / std::pow(gamma_km2-std::pow(pi2qsinRadius,2) ,2) ;
        if (m!=0)
          second_temp *= gamma_km2 / (gamma_km2-m*m);
        for ( int n = 0; n <=nMax; n += 1 ) 
          {
          // double fourth_temp = std::exp(- (gamma_km_Radius2_D0Delta+nPiLength2D0Delta[n]) ); 
          // double fourth_temp = utl::ExpNegtiveLUT(gamma_km_Radius2_D0Delta + nPiLength2D0Delta[n], 30.0, 1e3); 
          double fourth_temp = m_LUTExp->GetFunctionValue( - (gamma_km_Radius2_D0Delta+nPiLength2D0Delta[n])); 
          // double val = (gamma_km_Radius2_D0Delta+nPiLength2D0Delta[n]);
          // utlPrintVar3(true, val, utl::ExpNegtiveLUT(val, 30.0, 1e3), m_LUTExp->GetFunctionValue(-val));
          signal_temp += first_temp*second_temp*third_temp(m,n)*fourth_temp;
          }
        }
      }
    (*this->m_DWISamples)[i] = signal_temp;
    if (this->GetDebug())
      {
      std::cout << "i = " << i << std::endl << std::flush;
      std::cout << "sample = " << sample << std::endl << std::flush;
      std::cout << "qq = " << qq << std::endl << std::flush;
      std::cout << "signal_temp = " << signal_temp << std::endl << std::flush;
      }
    }
}

template <class PreciseType>
void
CylinderModelGenerator<PreciseType>
::ComputeEAPSamples ()
{
  this->VerifyInputParameters();
  utlGlobalException(this->m_SamplingSchemeRSpace->GetNumberOfSamples()==0, "need to set m_SamplingSchemeRSpace");

}

template <class PreciseType>
void
CylinderModelGenerator<PreciseType>
::ComputeODFSamples ()
{
  this->VerifyInputParameters();
  utlGlobalException(this->m_SamplingSchemeRSpace->GetNumberOfSamples()==0, "need to set m_SamplingSchemeRSpace");
}

}

#endif 



