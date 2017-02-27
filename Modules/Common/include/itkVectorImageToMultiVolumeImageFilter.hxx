/**
 *       @file  itkVectorImageToMultiVolumeImageFilter.hxx
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

#ifndef __itkVectorImageToMultiVolumeImageFilter_hxx
#define __itkVectorImageToMultiVolumeImageFilter_hxx


#include "utlITK.h"
#include "itkVectorImageToMultiVolumeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"


namespace itk
{

template <class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension>
void VectorImageToMultiVolumeImageFilter<TInputPixelType, TOutputPixelType, VImageDimension>
::GenerateOutputInformation()
{
  InputImagePointer  inputPtr = const_cast< InputImageType * >( this->GetInput());
  
  OutputImagePointer outputPtr = this->GetOutput();
  // NOTE: outputPtr->CopyInformation(inputPtr) is wrong when inputPtr is a 2D VectorImage
  itk::CopyImageInformation(inputPtr, outputPtr);
  
  unsigned int vectorDimension = inputPtr->GetNumberOfComponentsPerPixel();
  typename OutputImageType::RegionType outRegion = outputPtr->GetLargestPossibleRegion();
  typename OutputImageType::SizeType outSize = outRegion.GetSize();
  outSize[MultiVolumeImageDimension-1] = vectorDimension;
  outRegion.SetSize(outSize);
  outputPtr->SetRegions(outRegion);
}

template <class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension>
void VectorImageToMultiVolumeImageFilter<TInputPixelType, TOutputPixelType, VImageDimension>
::GenerateData()
{
  // get pointers to the input and output
  InputImagePointer  inputPtr = const_cast< InputImageType * >( this->GetInput());
  OutputImagePointer outputPtr = this->GetOutput();

  outputPtr->Allocate();
  
  ImageRegionIteratorWithIndex<InputImageType> inputIt(inputPtr, inputPtr->GetLargestPossibleRegion());
  inputIt.GoToBegin();
  InputImagePixelType inputPixel;
  unsigned int vectorDimension = inputPtr->GetNumberOfComponentsPerPixel();
  inputPixel.SetSize( vectorDimension );
  typename InputImageType::IndexType   inputIndex;
  typename OutputImageType::IndexType  outputIndex;

  while(!inputIt.IsAtEnd())
    {
    inputIndex = inputIt.GetIndex();
    inputPixel = inputIt.Get();

    for ( unsigned int i = 0; i < VectorImageDimension; i++ )
      {
      outputIndex[i] = inputIndex[i];
      }

    for ( unsigned int i = 0; i < vectorDimension; i++ )
      {
      outputIndex[MultiVolumeImageDimension-1] = i;
      outputPtr->SetPixel(outputIndex, inputPixel[i]);
      }

    ++inputIt;
    }
}

template <class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension>
void VectorImageToMultiVolumeImageFilter<TInputPixelType, TOutputPixelType, VImageDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar1(this->GetInput(), this->GetInput(), os<<indent << "VectorImage");
  PrintVar1(this->GetOutput(), this->GetOutput(), os<<indent << "MultiVolumeImage");
}

template <class ScalarType, unsigned Dim>
void 
VectorToMultiVolumeImage ( const SmartPointer<VectorImage<ScalarType,Dim> >& image1, SmartPointer<Image<ScalarType, Dim+1> >& image2 )
{
  typedef itk::VectorImageToMultiVolumeImageFilter<ScalarType, ScalarType> ConvertorType;
  typename ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(image1);
  convertor->Update();
  image2 = convertor->GetOutput();
}

}

#endif 
