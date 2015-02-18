/**
 *       @file  itkDiffusionModelEstimationInSphericalCoordinateImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "03-10-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkDiffusionModelEstimationInSphericalCoordinateImageFilter_hxx
#define __itkDiffusionModelEstimationInSphericalCoordinateImageFilter_hxx

#include "itkDiffusionModelEstimationInSphericalCoordinateImageFilter.h"

namespace itk
{

template< class TInputImage, class TOutputImage >
DiffusionModelEstimationInSphericalCoordinateImageFilter< TInputImage, TOutputImage >
::DiffusionModelEstimationInSphericalCoordinateImageFilter() : Superclass(), 
  m_BasisSHMatrix(new MatrixType()), 
  m_BasisRadialMatrix(new MatrixType())
{
  m_SHRank = 4;
  m_RadialRank = 0;
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
DiffusionModelEstimationInSphericalCoordinateImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  
  rval->m_SHRank = m_SHRank;
  rval->m_RadialRank = m_RadialRank;
  rval->m_BasisSHMatrix = m_BasisSHMatrix;
  rval->m_BasisRadialMatrix = m_BasisRadialMatrix;
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
DiffusionModelEstimationInSphericalCoordinateImageFilter< TInputImage, TOutputImage >
::ComputeSHMatrix()
{
  itkShowPositionThreadedLogger(this->GetDebug());

  MatrixPointer qOrientations = this->m_SamplingSchemeQSpace->GetOrientationsSpherical();
  this->m_BasisSHMatrix= utl::ComputeSHMatrix(m_SHRank, *qOrientations, SPHERICAL_TO_SPHERICAL);

  if(this->GetDebug()) 
    {
    int n_s = qOrientations->Rows();
    int n_b = (m_SHRank + 1)*(m_SHRank + 2)/2;
    std::ostringstream msg;
    msg << this->ThreadIDToString() << "Generated the "<< n_s << "x" << n_b << " Bmatrix...\n";
    utl::PrintUtlMatrix(*this->m_BasisSHMatrix,"SHMatrix", " ", msg << this->ThreadIDToString());
    this->WriteLogger(msg.str());
    }
}

template< class TInputImage, class TOutputImage >
void
DiffusionModelEstimationInSphericalCoordinateImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  Superclass::VerifyInputParameters();

  utlGlobalException(this->m_SHRank<0 || this->m_RadialRank<0, "negative rank");
  utlGlobalException(!utl::IsEven(this->m_SHRank), "sh rank should be even");
}

template< class TInputImage, class TOutputImage >
void
DiffusionModelEstimationInSphericalCoordinateImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  PrintVar2(true, m_SHRank, m_RadialRank, os<<indent);
  utl::PrintUtlMatrix(*m_BasisSHMatrix, "m_BasisSHMatrix", " ", os<<indent);
  utl::PrintUtlMatrix(*m_BasisRadialMatrix, "m_BasisRadialMatrix", " ", os<<indent);
}

}


#endif 



