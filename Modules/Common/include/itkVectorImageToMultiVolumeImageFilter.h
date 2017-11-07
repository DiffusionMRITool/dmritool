/**
 *       @file  itkVectorImageToMultiVolumeImageFilter.h
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

#ifndef __itkVectorImageToMultiVolumeImageFilter_h
#define __itkVectorImageToMultiVolumeImageFilter_h

#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{

/**
 *   \class   VectorImageToMultiVolumeImageFilter
 *   \brief   convert VectorImage<TOutputPixelType, VImageDimension> to Image<TInputPixelType, VImageDimension+1>
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension = 3>
class ITK_EXPORT VectorImageToMultiVolumeImageFilter :
    public ImageToImageFilter< VectorImage<TOutputPixelType, VImageDimension>, Image<TInputPixelType, VImageDimension+1>  >
{

public:
  
  itkStaticConstMacro (VectorImageDimension, unsigned int, VImageDimension);
  itkStaticConstMacro (MultiVolumeImageDimension, unsigned int, VImageDimension+1);

  typedef VectorImage<TOutputPixelType,VectorImageDimension>         InputImageType;
  typedef Image<TInputPixelType,MultiVolumeImageDimension>           OutputImageType;

  /** Standard class typedefs. */
  typedef VectorImageToMultiVolumeImageFilter                    Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType>    Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(VectorImageToMultiVolumeImageFilter, ImageToImageFilter);


  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename OutputImageType::PixelType   OutputImagePixelType;

protected:
  VectorImageToMultiVolumeImageFilter() {};
  virtual ~VectorImageToMultiVolumeImageFilter() {}
  
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  
  /** Does the real work. */
  virtual void GenerateData() ITK_OVERRIDE;

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

private:
  VectorImageToMultiVolumeImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkVectorImageToMultiVolumeImageFilter_hxx)
#include "itkVectorImageToMultiVolumeImageFilter.hxx"
#endif

#endif 

