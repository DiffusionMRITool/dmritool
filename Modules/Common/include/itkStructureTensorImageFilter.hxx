/**
 *       @file  itkStructureTensorImageFilter.hxx
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkStructureTensorImageFilter_hxx
#define __itkStructureTensorImageFilter_hxx

#include "itkStructureTensorImageFilter.h"
#include "itkOuterProductVectorImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkFunctors.h"

namespace itk 
{

template< class TInputImage, class TOutputImage >
StructureTensorImageFilter < TInputImage, TOutputImage >
::StructureTensorImageFilter() : Superclass() 
{
  m_IntensityScale=1.0;
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
StructureTensorImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_IntensityScale = m_IntensityScale;
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
StructureTensorImageFilter < TInputImage, TOutputImage >
::GenerateData()
{
  typename TInputImage::ConstPointer inputPtr = this->GetInput();
  typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
  
  if (m_IntensityScale!=1.0)
    {
    typedef itk::Functor::TIMES<typename TInputImage::PixelType, double> FunctorType;
    FunctorType timeFunctor;
    typedef UnaryFunctorImageFilter<TInputImage, TInputImage, FunctorType> ScaleFilterType;
    typename ScaleFilterType::Pointer scaleFilter = ScaleFilterType::New();
    timeFunctor.SetArgument(m_IntensityScale);
    scaleFilter->SetFunctor(timeFunctor);
    scaleFilter->SetInput(inputPtr);
    scaleFilter->Update();
    gradientFilter->SetInput(scaleFilter->GetOutput());
    }
  else
    gradientFilter->SetInput(inputPtr);

  gradientFilter->Update();

  typename TOutputImage::Pointer gradient = gradientFilter->GetOutput();

  typedef itk::OuterProductVectorImageFilter<TOutputImage, TOutputImage> OuterProductFilterType;
  typename OuterProductFilterType::Pointer outproductFilter = OuterProductFilterType::New();
  outproductFilter->SetInput(gradient);
  outproductFilter->Update();

  typename TOutputImage::Pointer out = outproductFilter->GetOutput();
  this->GraftOutput(out);
}


}


#endif 

