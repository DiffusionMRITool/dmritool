/**
 *       @file  itkAddNoiseToDWIImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-15-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkAddNoiseToDWIImageFilter_hxx
#define __itkAddNoiseToDWIImageFilter_hxx


#include "itkAddNoiseToDWIImageFilter.h"
#include "utl.h"


namespace itk
{

template <class TInputImage, class TB0Image, class TMaskImage>
AddNoiseToDWIImageFilter<TInputImage, TB0Image, TMaskImage>
::AddNoiseToDWIImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  m_Sigma = -1.0;
  m_Noisetype = RICIAN;
}


template <class TInputImage, class TB0Image, class TMaskImage>
void 
AddNoiseToDWIImageFilter<TInputImage, TB0Image, TMaskImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar1(true, m_Sigma, os<<indent);
  PrintVar1(true, m_Noisetype, os<<indent);
  PrintVar1(this->GetInput(), this->GetInput(), os<<indent << "dwi");
  PrintVar1(this->GetB0Image(), this->GetB0Image(), os<<indent << "b0");
  PrintVar1(this->GetMaskImage(), this->GetMaskImage(), os<<indent << "mask");
}


template<class TInputImage, class TB0Image, class TMaskImage>
void AddNoiseToDWIImageFilter<TInputImage, TB0Image, TMaskImage>
::GenerateData()
{
  utlException(m_Sigma<=0, "need to set sigma to generate noise");

  typename InputImageType::ConstPointer inputPtr = this->GetInput();
  typename MaskImageType::ConstPointer maskPtr = this->GetMaskImage();
  typename MaskImageType::ConstPointer b0Ptr = this->GetB0Image();
  utlException(!inputPtr, "no dwi");
 
  typename InputImageType::Pointer outputPtr = this->GetOutput();
  outputPtr->SetRequestedRegion(inputPtr->GetRequestedRegion());
  outputPtr->SetBufferedRegion(inputPtr->GetBufferedRegion());
  outputPtr->Allocate();
 
  ImageRegionIterator<InputImageType> outputIt(outputPtr, outputPtr->GetLargestPossibleRegion());
  ImageRegionConstIterator<InputImageType> inputIt(inputPtr, inputPtr->GetLargestPossibleRegion());
  ImageRegionConstIterator<MaskImageType> maskIt;
  
  if (maskPtr)
    maskIt = ImageRegionConstIterator<MaskImageType>(maskPtr, maskPtr->GetLargestPossibleRegion());
  
  itk::ImageRegionConstIterator<B0ImageType> b0It;
  if (b0Ptr)
    b0It = ImageRegionConstIterator<B0ImageType>(b0Ptr, b0Ptr->GetLargestPossibleRegion());
 
  outputIt.GoToBegin();
  inputIt.GoToBegin();
  maskIt.GoToBegin();
  b0It.GoToBegin();
  
  PixelType pixel, zeroPixel;
  zeroPixel.SetSize(inputPtr->GetNumberOfComponentsPerPixel());
  zeroPixel.Fill(0);
  
  double realSigma = m_Sigma;
 
  while(!outputIt.IsAtEnd())
    {
    if ((!maskPtr || (maskPtr && maskIt.Get()>0) )
        && (!b0Ptr || (b0Ptr && b0It.Get()>0)))
      {
      pixel = inputIt.Get();
      if (b0Ptr)
        realSigma = m_Sigma * b0It.Get();
      if (m_Noisetype==GAUSSIAN)
        pixel = utl::AddNoise<PixelType>(pixel, pixel.GetSize(), realSigma, false);
      if (m_Noisetype==RICIAN)
        pixel = utl::AddNoise<PixelType>(pixel, pixel.GetSize(), realSigma, true);
      outputIt.Set(pixel);
      }
    else
      outputIt.Set(zeroPixel);

    ++inputIt;
    ++outputIt;
    
    if (b0Ptr)
      ++b0It;
    if (maskPtr)
      ++maskIt;
    } 
}

}

#endif 
