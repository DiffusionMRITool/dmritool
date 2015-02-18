/**
 *       @file  itkMultiVolumeImageToVectorImageFilter.h
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

#ifndef __itkMultiVolumeImageToVectorImageFilter_h
#define __itkMultiVolumeImageToVectorImageFilter_h

#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{

/**
 *   \class   MultiVolumeImageToVectorImageFilter
 *   \brief   convert Image<TInputPixelType, VImageDimension+1> to VectorImage<TOutputPixelType, VImageDimension> 
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputPixelType, class TOutputPixelType, unsigned int VImageDimension = 3>
class ITK_EXPORT MultiVolumeImageToVectorImageFilter :
    public ImageToImageFilter< Image<TInputPixelType,VImageDimension+1>, VectorImage<TOutputPixelType, VImageDimension> >
{

public:
  
  itkStaticConstMacro (VectorImageDimension, unsigned int, VImageDimension);
  itkStaticConstMacro (MultiVolumeImageDimension, unsigned int, VImageDimension+1);

  typedef Image<TInputPixelType,MultiVolumeImageDimension>           InputImageType;
  typedef VectorImage<TOutputPixelType,VectorImageDimension>         OutputImageType;

  /** Standard class typedefs. */
  typedef MultiVolumeImageToVectorImageFilter                    Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType>    Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MultiVolumeImageToVectorImageFilter, ImageToImageFilter);


  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;
  
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename OutputImageType::PixelType   OutputImagePixelType;

protected:
  MultiVolumeImageToVectorImageFilter() {};
  virtual ~MultiVolumeImageToVectorImageFilter() {}
  
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Does the real work. */
  virtual void GenerateData();

  virtual void GenerateOutputInformation();

private:
  MultiVolumeImageToVectorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkMultiVolumeImageToVectorImageFilter_hxx)
#include "itkMultiVolumeImageToVectorImageFilter.hxx"
#endif

#endif 

