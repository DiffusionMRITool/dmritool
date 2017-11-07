/*=========================================================================

 Program:   Image to VTK Image Data Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkImageToVTKImageDataFilter_h
#define __itkImageToVTKImageDataFilter_h

#include "itkProcessObject.h"
#include "itkImageRegionConstIterator.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"

namespace itk
{

/** \class ImageToVTKImageDataFilter
 * \brief Converts an ITK image to VTK image data.
 */
template <class TInputImage >
class ITK_EXPORT ImageToVTKImageDataFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ImageToVTKImageDataFilter    Self;
  typedef ProcessObject                Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToVTKImageDataFilter, ProcessObject);

  /** Image Dimension */
  itkStaticConstMacro(Dimension, unsigned int, TInputImage::ImageDimension);

  /** Some typedefs. */
  typedef TInputImage                               InputImageType;
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename InputImageType::RegionType       InputRegionType;
  typedef typename InputImageType::SizeType         InputSizeType;
  typedef typename InputImageType::SpacingType      InputSpacingType;
  typedef typename InputImageType::PointType        InputPointType;
  typedef typename InputImageType::IndexType        InputIndexType;
  typedef typename InputImageType::PixelType        InputPixelType;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef ImageRegionConstIterator<InputImageType>  InputImageIteratorType;

  typedef vtkSmartPointer<vtkImageData>             OutputImageType;

  /** Get the output in the form of a vtkImage.
      This call is delegated to the internal vtkImageImporter filter  */
  vtkImageData *  GetOutput() const;

  /** Set the input */
  using Superclass::SetInput;
  virtual void SetInput(const InputImageType *image);
  const InputImageType * GetInput() const;

  virtual void Update() ITK_OVERRIDE
    {
    this->GenerateData();
    }

protected:
  ImageToVTKImageDataFilter();
  virtual ~ImageToVTKImageDataFilter() {};

  virtual void GenerateData() ITK_OVERRIDE;
  
private:
  ImageToVTKImageDataFilter(const Self&); //purposely not implemented
  void operator=(const Self&);        //purposely not implemented

  OutputImageType m_Image;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToVTKImageDataFilter.hxx"
#endif

#endif
