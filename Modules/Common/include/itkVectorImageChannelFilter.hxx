/**
 *       @file  itkVectorImageChannelFilter.hxx
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkVectorImageChannelFilter_hxx
#define __itkVectorImageChannelFilter_hxx

#include "itkVectorImageChannelFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkComposeVectorImageFilter.h"

namespace itk
{

template< class TInputImage, class TOutputImage, class TFilter >
VectorImageChannelFilter< TInputImage, TOutputImage, TFilter >
::VectorImageChannelFilter() : Superclass() 
{
  m_Filter=NULL;
}

template< class TInputImage, class TOutputImage, class TFilter >
typename LightObject::Pointer
VectorImageChannelFilter< TInputImage, TOutputImage, TFilter >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_Filter = m_Filter->Clone();
  return loPtr;
}

template< class TInputImage, class TOutputImage, class TFilter >
void
VectorImageChannelFilter< TInputImage, TOutputImage, TFilter >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  if (m_Filter)
    m_Filter->Print(os<<indent << "m_Filter =\n");
}

template< class TInputImage, class TOutputImage, class TFilter >
void
VectorImageChannelFilter< TInputImage, TOutputImage, TFilter >
::GenerateData()
{
  utlException(!m_Filter, "need to set m_Filter");
  typename TInputImage::ConstPointer inputImage = this->GetInput();
  int N = inputImage->GetNumberOfComponentsPerPixel();

  typedef typename TFilter::OutputImageType FilterOutputImageType;
  typedef itk::ComposeVectorImageFilter<FilterOutputImageType, OutputImageType> ComposeImageFilterType;
  typename ComposeImageFilterType::Pointer composeImageFilter = ComposeImageFilterType::New();

  typedef itk::VectorIndexSelectionCastImageFilter<TInputImage, typename TFilter::InputImageType> IndexSelectionType;

  for ( int i = 0; i < N; ++i ) 
    {
    typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(inputImage);
    // indexSelectionFilter->Update();
    
    typename TFilter::Pointer filterCopy = m_Filter->Clone();
    filterCopy->SetInput(indexSelectionFilter->GetOutput());
    filterCopy->Update();
    typename FilterOutputImageType::Pointer out = filterCopy->GetOutput();

    composeImageFilter->SetInput(i, filterCopy->GetOutput());
    }

  composeImageFilter->Update();
  typename OutputImageType::Pointer outputImage = composeImageFilter->GetOutput();

  this->GraftOutput(outputImage);
}


}

#endif 

