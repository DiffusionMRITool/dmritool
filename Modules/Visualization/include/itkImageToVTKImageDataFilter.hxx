/*=========================================================================

 Program:   Image to VTK Image Data Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkImageToVTKImageDataFilter_hxx
#define __itkImageToVTKImageDataFilter_hxx

#include "itkImageToVTKImageDataFilter.h"

#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage>
ImageToVTKImageDataFilter<TInputImage>
::ImageToVTKImageDataFilter()
{
  m_Image = NULL;
}

/**
 * Get a vtkImage as output
 */
template <class TInputImage>
vtkImageData *
ImageToVTKImageDataFilter<TInputImage>
::GetOutput() const
{
  return m_Image;
}

/**
 * Get a vtkImage as output
 */
template <class TInputImage>
void
ImageToVTKImageDataFilter<TInputImage>
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
const typename ImageToVTKImageDataFilter<TInputImage>::InputImageType *
ImageToVTKImageDataFilter<TInputImage>
::GetInput() const
{
  return itkDynamicCastInDebugMode< const TInputImage * >(this->GetPrimaryInput());
}

/**
 * Generate the data
 */
template <class TInputImage>
void
ImageToVTKImageDataFilter<TInputImage>
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
  vtkSmartPointer<vtkDoubleArray> data =
    vtkSmartPointer<vtkDoubleArray>::New();

  data->SetNumberOfComponents(numberOfComponentsPerPixel);
  data->SetNumberOfTuples(numberOfVoxels);
  data->SetName("scalars");

  // Iterator
  InputImageIteratorType
    inputIt(inputImage, inputImage->GetLargestPossibleRegion());

  inputIt.GoToBegin();

  // Populate data
  InputPixelType inputPixel;
  vtkIdType i = 0;
  double dataEntry[numberOfComponentsPerPixel];

  while(!inputIt.IsAtEnd())
    {
    inputPixel = inputIt.Get();
    for (unsigned int k = 0; k < numberOfComponentsPerPixel; k++)
      {
      dataEntry[k] = static_cast<double>( inputPixel[k] );
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
