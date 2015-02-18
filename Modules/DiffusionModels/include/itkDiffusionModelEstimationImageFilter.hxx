/**
 *       @file  itkDiffusionModelEstimationImageFilter.hxx
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


#ifndef __itkDiffusionModelEstimationImageFilter_hxx
#define __itkDiffusionModelEstimationImageFilter_hxx

#include "itkDiffusionModelEstimationImageFilter.h"

namespace itk
{

template< class TInputImage, class TOutputImage >
DiffusionModelEstimationImageFilter< TInputImage, TOutputImage >
::DiffusionModelEstimationImageFilter() : Superclass(), 
  m_BasisMatrix(new MatrixType()), 
  m_RegularizationWeight(new VectorType())
{
  m_MD0 = 0.7e-3;
  m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
DiffusionModelEstimationImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_SamplingSchemeQSpace = m_SamplingSchemeQSpace->Clone();
  
  // NOTE: shared_ptr is thread safe, if the data is read only, thus do not need to copy the data block. 
  rval->m_RegularizationWeight = m_RegularizationWeight;
  rval->m_BasisMatrix = m_BasisMatrix;

  rval->m_MD0 = m_MD0;
  
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
DiffusionModelEstimationImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  if (!IsImageEmpty(this->m_MaskImage))
    this->VerifyMaskInformation();
  
  utlGlobalException(m_SamplingSchemeQSpace->GetBVector()->size()==0 && m_SamplingSchemeQSpace->GetRadiusVector()->size()==0, "no b values nor q values");
  utlGlobalException(m_SamplingSchemeQSpace->size()==0, "no gradients in q-space");
  InputImageConstPointer inputPtr = this->GetInput();
  utlGlobalException(!inputPtr, "no input DWIs");
  int dwiNumber = inputPtr->GetNumberOfComponentsPerPixel();
  utlGlobalException(dwiNumber!=m_SamplingSchemeQSpace->GetNumberOfSamples(), "the size of gradients and the size of DWIs are not the same");
  utlGlobalException(m_SamplingSchemeQSpace->GetRadiusVector()->size()!=m_SamplingSchemeQSpace->GetNumberOfSamples(), "the size of gradients and the size of b values are not the same");
}

template< class TInputImage, class TOutputImage >
void
DiffusionModelEstimationImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar1(true, m_MD0, os<<indent);
  
  os << indent << "m_SamplingSchemeQSpace = " << m_SamplingSchemeQSpace << std::endl;
  utl::PrintUtlMatrix(*m_BasisMatrix, "m_BasisMatrix", " ", os<<indent);
  utl::PrintUtlVector(*m_RegularizationWeight, "m_RegularizationWeight", " ", os<<indent);
}

}

#endif 


