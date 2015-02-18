/**
 *       @file  itkMultiVolumeImageToVectorImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-15-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkMultiVolumeImageToVectorImageFilter_hxx
#define __itkMultiVolumeImageToVectorImageFilter_hxx


#include "utlITK.h"
#include "itkMultiVolumeImageToVectorImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"


namespace itk
{

template <class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension>
void MultiVolumeImageToVectorImageFilter<TInputPixelType, TOutputPixelType, VImageDimension>
::GenerateOutputInformation()
{
  InputImagePointer  inputPtr = const_cast< InputImageType * >( this->GetInput());
  OutputImagePointer outputPtr = this->GetOutput();
  itk::CopyImageInformation<InputImageType, OutputImageType>(inputPtr, outputPtr);
  typename InputImageType::RegionType inputRegion = inputPtr->GetLargestPossibleRegion();
  typename InputImageType::SizeType inputImageSize = inputRegion.GetSize();
  outputPtr->SetNumberOfComponentsPerPixel( inputImageSize[MultiVolumeImageDimension-1] );
}

template <class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension>
void MultiVolumeImageToVectorImageFilter<TInputPixelType, TOutputPixelType, VImageDimension>
::GenerateData()
{
  // get pointers to the input and output
  InputImagePointer  inputPtr = const_cast< InputImageType * >( this->GetInput());
  OutputImagePointer outputPtr = this->GetOutput();

  outputPtr->Allocate();

  ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputPtr->GetLargestPossibleRegion());
  outputIt.GoToBegin();
  OutputImagePixelType outputPixel;
  unsigned int vectorDimension = outputPtr->GetNumberOfComponentsPerPixel();
  outputPixel.SetSize( vectorDimension );
  typename InputImageType::IndexType   inputIndex;
  typename OutputImageType::IndexType  outputIndex;

  while(!outputIt.IsAtEnd())
    {
    outputIndex = outputIt.GetIndex();
    for ( int i = 0; i < VectorImageDimension; i += 1 ) 
      {
      inputIndex[i] = outputIndex[i];
      }
    for ( int i = 0; i < vectorDimension; i += 1 )
      {
      inputIndex[MultiVolumeImageDimension-1] = i;
      outputPixel[i] = (TOutputPixelType) inputPtr->GetPixel(inputIndex);
      }
    outputIt.Set(outputPixel);
    ++outputIt;
    }
}

template <class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension>
void MultiVolumeImageToVectorImageFilter<TInputPixelType, TOutputPixelType, VImageDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar1(this->GetInput(), this->GetInput(), os<<indent << "MultiVolumeImage");
  PrintVar1(this->GetOutput(), this->GetOutput(), os<<indent << "VectorImage");
}

}

#endif 
