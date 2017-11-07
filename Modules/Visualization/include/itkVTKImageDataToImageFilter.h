/*=========================================================================

 Program:   VTK Image Data to Image Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkVTKImageDataToImageFilter_h
#define __itkVTKImageDataToImageFilter_h

#include "itkImageSource.h"
#include "itkImageRegionIterator.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"

namespace itk
{

/** \class VTKImageDataToImageFilter
 * \brief Convert VTK image data to an ITK image.
 */
template <class TOutputImage>
class ITK_EXPORT VTKImageDataToImageFilter : 
  public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef VTKImageDataToImageFilter          Self;
  typedef ImageSource<TOutputImage>         Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VTKImageDataToImageFilter, ImageSource);

  /** Image Dimension */
  itkStaticConstMacro(Dimension, unsigned int, TOutputImage::ImageDimension);

  /** Some typedefs. */
  typedef TOutputImage                              OutputImageType;
  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename OutputImageType::RegionType      OutputRegionType;
  typedef typename OutputImageType::SizeType        OutputSizeType;
  typedef typename OutputImageType::SpacingType     OutputSpacingType;
  typedef typename OutputImageType::PointType       OutputPointType;
  typedef typename OutputImageType::IndexType       OutputIndexType;
  typedef typename OutputImageType::PixelType       OutputPixelType;
  typedef ImageRegionIterator<OutputImageType>      OutputImageIteratorType;
  typedef vtkSmartPointer<vtkImageData>             InputImageType;

  /** Set the input in the form of a vtkImageData */
  void SetInputData( vtkImageData * );

protected:
  VTKImageDataToImageFilter();
  virtual ~VTKImageDataToImageFilter() {};

  virtual void GenerateData() ITK_OVERRIDE;

private:
  VTKImageDataToImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  InputImageType m_Image;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVTKImageDataToImageFilter.hxx"
#endif

#endif
