/*=========================================================================

 Program:   RGB Image to VTK Image Data Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkRGBImageToVTKImageDataFilter_hxx
#define __itkRGBImageToVTKImageDataFilter_hxx

#include "itkRGBImageToVTKImageDataFilter.h"

#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPointData.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage>
RGBImageToVTKImageDataFilter<TInputImage>
::RGBImageToVTKImageDataFilter()
{
  m_Image = NULL;
}

/**
 * Get a vtkImage as output
 */
template <class TInputImage>
vtkImageData *
RGBImageToVTKImageDataFilter<TInputImage>
::GetOutput() const
{
  return m_Image;
}

/**
 * Get a vtkImage as output
 */
template <class TInputImage>
void
RGBImageToVTKImageDataFilter<TInputImage>
::SetInput( const InputImageType *image )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
    const_cast< InputImageType * >( image ) );
}

/**
 * Get a vtkImage as output
 */
template <class TInputImage>
const typename RGBImageToVTKImageDataFilter<TInputImage>::InputImageType *
RGBImageToVTKImageDataFilter<TInputImage>
::GetInput() const
{
  return itkDynamicCastInDebugMode< const TInputImage * >(this->GetPrimaryInput());
}

/**
 * Generate the data
 */
template <class TInputImage>
void
RGBImageToVTKImageDataFilter<TInputImage>
::GenerateData()
{
  // Get input image
  InputImageConstPointer inputImage = this->GetInput();

  // Image information
  InputSizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
  InputSpacingType inputSpacing = inputImage->GetSpacing();
  InputPointType inputOrigin = inputImage->GetOrigin();
  unsigned int numberOfComponentsPerPixel
      = inputImage->GetNumberOfComponentsPerPixel();
  unsigned int rgbaSize = 4;
  unsigned int numberOfComponents
      = vnl_math_min(rgbaSize, numberOfComponentsPerPixel);

  // Image data
  m_Image = OutputImageType::New();

  int outputSize[Dimension];
  double outputSpacing[Dimension];
  double outputOrigin[Dimension];
  vtkIdType numberOfVoxels = 1;

  for (unsigned int k=0; k<Dimension; k++)
    {
    outputSize[k] = inputSize[k];
    outputSpacing[k] = inputSpacing[k];
    outputOrigin[k] = inputOrigin[k];
    numberOfVoxels *= outputSize[k];
    }

  m_Image->SetDimensions(outputSize);
  m_Image->SetSpacing(outputSpacing);
  m_Image->SetOrigin(outputOrigin);

  // Data array
  vtkSmartPointer<vtkUnsignedCharArray> data =
    vtkSmartPointer<vtkUnsignedCharArray>::New();

  data->SetNumberOfComponents(numberOfComponents);
  data->SetNumberOfTuples(numberOfVoxels);
  data->SetName("RGB_Vectors");

  // Iterator
  InputImageIteratorType
    inputIt(inputImage, inputImage->GetLargestPossibleRegion());

  inputIt.GoToBegin();

  // Populate data
  InputPixelType inputPixel;
  vtkIdType i = 0;
  double dataEntry [numberOfComponents];

  while(!inputIt.IsAtEnd())
    {
    inputPixel = inputIt.Get();
    for (unsigned int k = 0; k < rgbaSize - 1; k++)
      {
      if ( k < numberOfComponents )
        {
        dataEntry[k] = static_cast<unsigned char>( inputPixel[k] );
        }
      else
        {
        dataEntry[k] = static_cast<unsigned char>( 0 );
        }
      }
    if ( numberOfComponents == rgbaSize )
      {
      dataEntry[rgbaSize - 1]
          = static_cast<unsigned char>( inputPixel[rgbaSize - 1] );
      }
    else
      {
      dataEntry[rgbaSize - 1] = 1.0;
      }
    data->SetTuple(i, dataEntry);
    ++inputIt;
    i++;
    }

  m_Image->GetPointData()->AddArray(data);
  m_Image->GetPointData()->SetActiveScalars(data->GetName());
}

} // end namespace itk

#endif
