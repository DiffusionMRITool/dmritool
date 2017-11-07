/**
 *       @file  itkVectorImageChannelFilter.h
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkVectorImageChannelFilter_h
#define __itkVectorImageChannelFilter_h


#include "itkImageToImageFilter.h"

namespace itk
{

/**
 *   \class   VectorImageChannelFilter
 *   \brief   TInputImage and TOutputImage are VectorImage, TFilter is an image filter which works for itk::Image. 
 *   This class performs TFilter on each channel of TInputImage, then compose results into output. 
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template< class TInputImage, class TOutputImage, 
    class TFilter=ImageToImageFilter< Image<typename TInputImage::InternalPixelType, TInputImage::ImageDimension>, TOutputImage> >
class ITK_EXPORT VectorImageChannelFilter
  :public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef VectorImageChannelFilter                        Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorImageChannelFilter, ImageToImageFilter);
  
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::IndexType      InputImageIndexType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SpacingType    InputImageSpacingType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  
  typedef TOutputImage                            OutputImageType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename OutputImageType::IndexType     OutputImageIndexType;
  typedef typename OutputImageType::SizeType      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType   OutputImageSpacingType;
  typedef typename OutputImageType::PixelType     OutputImagePixelType;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  
  itkStaticConstMacro(InputImageDimension, unsigned int,  TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  itkSetMacro(Filter, typename TFilter::Pointer);
  itkGetMacro(Filter, typename TFilter::Pointer);

protected:
  VectorImageChannelFilter();
  virtual ~VectorImageChannelFilter() {};

  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;

  typename TFilter::Pointer m_Filter;


private:
  VectorImageChannelFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};

}

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkVectorImageChannelFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkVectorImageChannelFilter_hxx)
#include "itkVectorImageChannelFilter.hxx"
#endif

#endif 
