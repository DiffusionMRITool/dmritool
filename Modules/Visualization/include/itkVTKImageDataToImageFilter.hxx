/*=========================================================================

 Program:   VTK Image Data to Image Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkVTKImageDataToImageFilter_hxx
#define __itkVTKImageDataToImageFilter_hxx

#include "itkVTKImageDataToImageFilter.h"
 
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkVersion.h"

namespace itk
{

/**
 * Constructor
 */
template <class TOutputImage>
VTKImageDataToImageFilter<TOutputImage>
::VTKImageDataToImageFilter()
{
  m_Image = NULL;
}

/**
 * Set a vtkImageData as input
 */
template <class TOutputImage>
void
VTKImageDataToImageFilter<TOutputImage>
::SetInputData( vtkImageData * inputImage )
{
  if ( inputImage != m_Image )
    {
    m_Image = inputImage;
//    m_Image->Register(this);
    this->Modified();
    }
}

/**
 * Generate the data
 */
template <class TOutputImage>
void
VTKImageDataToImageFilter<TOutputImage>
::GenerateData()
{
  // Data array
  vtkSmartPointer<vtkDoubleArray> data 
    = vtkSmartPointer<vtkDoubleArray>::New();
  data = static_cast<vtkDoubleArray *>(
    m_Image->GetPointData()->GetScalars() );

  // Image information
  int inputSize[Dimension];
  double inputSpacing[Dimension];
  double inputOrigin[Dimension];
  m_Image->GetDimensions(inputSize);
  m_Image->GetSpacing(inputSpacing);
  m_Image->GetOrigin(inputOrigin);
  int numberOfComponentsPerPixel
    = m_Image->GetPointData()->GetNumberOfComponents();

  OutputImagePointer outputImage = this->GetOutput();

  // Setup output image
  OutputSizeType outputSize;
  OutputSpacingType outputSpacing;
  OutputPointType outputOrigin;
  OutputIndexType outputIndex;
  OutputPixelType outputPixel;
  outputPixel.SetSize( numberOfComponentsPerPixel );

  for (unsigned int k=0; k<Dimension; k++)
    {
    outputSize[k] = inputSize[k];
    outputSpacing[k] = inputSpacing[k];
    outputOrigin[k] = inputOrigin[k];
    outputIndex[k] = 0;
    }

  OutputRegionType outputRegion;
  outputRegion.SetSize(outputSize);
  outputRegion.SetIndex(outputIndex);

  outputImage->SetSpacing(outputSpacing);
  outputImage->SetRegions(outputRegion);
  outputImage->SetOrigin(outputOrigin);
  outputImage->SetNumberOfComponentsPerPixel(numberOfComponentsPerPixel);

  outputImage->Allocate();
  outputPixel.Fill(0);
  outputImage->FillBuffer(outputPixel);

  // Iterator
  OutputImageIteratorType
    outputIt(outputImage, outputImage->GetLargestPossibleRegion());

  outputIt.GoToBegin();

  // Populate data
  vtkIdType i = 0;
  double *dataEntry;

  while( !outputIt.IsAtEnd() )
    {
    dataEntry = data->GetTuple(i);
    for (unsigned int k = 0; k < numberOfComponentsPerPixel; k++)
      {
      outputPixel[k] = dataEntry[k];
      }
    outputIt.Set( outputPixel );
    ++outputIt;
    i++;
    }  
}

} // end namespace itk

#endif
