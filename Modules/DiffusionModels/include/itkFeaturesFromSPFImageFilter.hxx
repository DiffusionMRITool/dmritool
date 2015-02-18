/**
 *       @file  itkFeaturesFromSPFImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-08-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkFeaturesFromSPFImageFilter_hxx
#define __itkFeaturesFromSPFImageFilter_hxx

#include "itkFeaturesFromSPFImageFilter.h"
#include "itkSphericalPolarFourierImageFilter.h"
#include "utl.h"

namespace itk
{  

template< class TInputImage, class TOutputImage >
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::FeaturesFromSPFImageFilter ( ) : 
  m_SPFToFeatureTransform(new MatrixType()), 
  m_Orientations(new MatrixType())
{
  m_Tau = ONE_OVER_4_PI_2;
  m_MD0 = 0.7e-3;
  m_SHRank=-1;
  m_RadialRank=-1;
  m_IsFourier=true;
  // m_IsOriginalBasis=true;
  m_ScaleImage=ScalarImageType::New();
  m_BasisType=SPF;

  m_IsInQSpace=true;

  m_BasisScale=-1;
  m_SPFIEstimator=SPFIFilterBaseType::New();
}

template< class TInputImage, class TOutputImage >
double
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::ComputeScale(const bool setScale)
{
  double scale=-1;
  if (m_BasisType==SPF)
    scale = 1.0 / (8*M_PI*M_PI*this->m_Tau*this->m_MD0);  // 714.29 (700) scale for SPF basis, dual scale is 1/(4*pi^2*scale)
  else if (m_BasisType==DSPF)
    scale = 2*this->m_Tau*this->m_MD0;  // 3.5462e-5 scale for SPF basis, dual scale is 1/(4*pi^2*scale)
  if (setScale)
    m_BasisScale = scale;
  return scale;
}

template< class TInputImage, class TOutputImage >
void
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::SetBasisScale(const double scale)
{
  double scale_old = m_BasisScale;  
  if (scale>0)
    m_BasisScale = scale;
  else
    this->ComputeScale(true);
  itkDebugMacro("setting m_BasisScale to " << m_BasisScale);

  if (scale>0 && std::fabs((scale_old-m_BasisScale)/m_BasisScale)>1e-8)
    {
    this->Modified();
    m_SPFToFeatureTransform=MatrixPointer(new MatrixType());
    }
}

  
template< class TInputImage, class TOutputImage >
void 
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::SetScaleImage(const ScalarImagePointer& scaleImage)
{
  if( scaleImage != m_ScaleImage )
    {
    std::cout << "Use adaptive scale" << std::endl << std::flush;
    m_ScaleImage = scaleImage;
    m_SPFToFeatureTransform=MatrixPointer(new MatrixType());
    this->Modified();
    }
}

template< class TInputImage, class TOutputImage >
void
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters ( )
{
  utlGlobalException(m_SHRank<0 || m_RadialRank<0, "negative rank");
  typename TInputImage::ConstPointer  inputPtr = this->GetInput();
  utlSAGlobalException(this->m_SPFIEstimator->RankToDim(false, this->m_RadialRank, this->m_SHRank)!=inputPtr->GetNumberOfComponentsPerPixel())
    (this->m_RadialRank)(this->m_SHRank)(this->m_SPFIEstimator->RankToDim(false, this->m_RadialRank, this->m_SHRank))(inputPtr->GetNumberOfComponentsPerPixel()).msg("the size of input image is not consistent with the given shRank and radialRank");
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_BasisScale = m_BasisScale;
  rval->m_MD0 = m_MD0;
  rval->m_Tau = m_Tau;
  rval->m_SHRank = m_SHRank;
  rval->m_RadialRank = m_RadialRank;
  rval->m_ScaleImage = m_ScaleImage;

  rval->m_IsFourier = m_IsFourier;
  rval->m_IsInQSpace = m_IsInQSpace;
  rval->m_Orientations = m_Orientations;

  // NOTE: shared_ptr is thread safe if different threads read the same data block.  
  // However if the data is modified in different threads, then it needs to copy the data block. 
  rval->m_SPFToFeatureTransform = m_SPFToFeatureTransform;
  rval->m_BasisType = m_BasisType;

  rval->m_SPFIEstimator = m_SPFIEstimator->Clone();
  rval->SetDebug(this->GetDebug());
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData ( )
{
  std::cout << "Use " << this->GetNumberOfThreads() << " threads!" << std::endl << std::flush;
  VerifyInputParameters();
}

template< class TInputImage, class TOutputImage >
void
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::SetSPFIEstimator ( )
{
  if (this->m_BasisType==SPF || this->m_BasisType==DSPF)
    {
    typedef SphericalPolarFourierImageFilter<TInputImage, TOutputImage> SPFIFilterType;
    this->m_SPFIEstimator = SPFIFilterType::New();

    if (this->m_BasisType==SPF && !this->m_IsFourier)
      this->m_SPFIEstimator->SetIsOriginalBasis(true);
    else
      this->m_SPFIEstimator->SetIsOriginalBasis(false);
    }
}

template< class TInputImage, class TOutputImage >
void 
FeaturesFromSPFImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent
    << "m_Tau = " << m_Tau
    << ", m_BasisScalar = " << m_BasisScale
    << ", m_MD0 = " << m_MD0
    << ", m_SHRank = " << m_SHRank
    << ", m_RadialRank = " << m_RadialRank
    << ", m_IsInQSpace = " << m_IsInQSpace
    << ", m_IsFourier = " << this->m_IsFourier
    // << ", m_IsOriginalBasis = " << this->GetIsOriginalBasis()
    // << ", m_ScaleImage = " << m_ScaleImage
    << std::endl;
  utl::PrintUtlMatrix(*m_Orientations, "m_Orientations", " ",os << indent);
  utl::PrintUtlMatrix(*m_SPFToFeatureTransform, "m_SPFToFeatureTransform", " ",os << indent);
  if (!itk::IsImageEmpty(m_ScaleImage))
    std::cout << "m_ScaleImage = " << m_ScaleImage << std::endl << std::flush;
  if (m_BasisType==SPF)
    os << indent << "Use SPF basis" << std::endl << std::flush;
  else if (m_BasisType==DSPF)
    os << indent << "Use DSPF basis" << std::endl << std::flush;
  os << indent << "m_SPFIEstimator = " << m_SPFIEstimator << std::endl << std::flush;
}

}

#endif 

