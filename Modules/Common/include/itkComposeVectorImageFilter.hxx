/**
 *       @file  itkComposeVectorImageFilter.hxx
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkComposeVectorImageFilter_hxx
#define __itkComposeVectorImageFilter_hxx

#include "itkComposeVectorImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace itk
{
//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
ComposeVectorImageFilter< TInputImage, TOutputImage >
::ComposeVectorImageFilter()
{
  OutputPixelType p;
  int nbOfComponents = NumericTraits<OutputPixelType>::GetLength(p);
  nbOfComponents = std::max( 1, nbOfComponents );  // require at least one input
  this->SetNumberOfRequiredInputs( nbOfComponents );
}

//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
void
ComposeVectorImageFilter< TInputImage, TOutputImage >
::SetInput1(const InputImageType *image1)
{
  // The ProcessObject is not const-correct so the const_cast is required here
  this->SetNthInput( 0, const_cast< InputImageType * >( image1 ) );
}

//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
void
ComposeVectorImageFilter< TInputImage, TOutputImage >
::SetInput2(const InputImageType *image2)
{
  // The ProcessObject is not const-correct so the const_cast is required here
  this->SetNthInput( 1, const_cast< InputImageType * >( image2 ) );
}

//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
void
ComposeVectorImageFilter< TInputImage, TOutputImage >
::SetInput3(const InputImageType *image3)
{
  // The ProcessObject is not const-correct so the const_cast is required here
  this->SetNthInput( 2, const_cast< InputImageType * >( image3 ) );
}

//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
void
ComposeVectorImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation(void)
{
  // Override the method in itkImageSource, so we can set the vector length of
  // the output itk::VectorImage

  this->Superclass::GenerateOutputInformation();

  OutputImageType *output = this->GetOutput();
  int n=0;
  for ( int i = 0; i < this->GetNumberOfIndexedInputs(); ++i ) 
    {
    const InputImageType * inputImage = this->GetInput(i);
    n += inputImage->GetNumberOfComponentsPerPixel();
    }
  output->SetNumberOfComponentsPerPixel( n );
}

//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
void
ComposeVectorImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  // Check to verify all inputs are specified and have the same metadata,
  // spacing etc...
  const unsigned int numberOfInputs = this->GetNumberOfIndexedInputs();
  RegionType         region;

  for ( unsigned int i = 0; i < numberOfInputs; i++ )
    {
    InputImageType *input = itkDynamicCastInDebugMode< InputImageType * >
      (this->ProcessObject::GetInput(i) );
    if ( !input )
      {
      itkExceptionMacro(<< "Input " << i << " not set!");
      }
    if ( i == 0 )
      {
      region = input->GetLargestPossibleRegion();
      }
    else if ( input->GetLargestPossibleRegion() != region )
      {
      itkExceptionMacro(<< "All Inputs must have the same dimensions.");
      }
    }
}

//----------------------------------------------------------------------------
template< typename TInputImage, typename TOutputImage >
void
ComposeVectorImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  typename OutputImageType::Pointer outputImage =
    static_cast< OutputImageType * >( this->ProcessObject::GetOutput(0) );

  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  InputIteratorContainerType inputItContainer;

  for ( unsigned int i = 0; i < this->GetNumberOfIndexedInputs(); i++ )
    {
    const InputImageType * inputImage = this->GetInput(i);

    InputIteratorType iit( inputImage, outputRegionForThread);
    iit.GoToBegin();
    inputItContainer.push_back(iit);
    }

  OutputPixelType pix;
  NumericTraits<OutputPixelType>::SetLength( pix, outputImage->GetNumberOfComponentsPerPixel() );
  while ( !oit.IsAtEnd() )
    {
    ComputeOutputPixel( pix, inputItContainer );
    oit.Set(pix);
    ++oit;
    progress.CompletedPixel();
    }
}
} // end namespace itk

#endif



