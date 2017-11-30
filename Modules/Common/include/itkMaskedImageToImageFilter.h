/**
 *       @file  itkMaskedImageToImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-20-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkMaskedImageToImageFilter_h
#define __itkMaskedImageToImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkThreadLogger.h"
#include "utlSTDHeaders.h"
#include "utlITKMacro.h"


namespace itk
{

/**
 *   \class   MaskedImageToImageFilter
 *   \brief   ImageToImageFilter with mask and threaded logger support
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template< class TInputImage, class TOutputImage, class TMaskImage=Image<double, 3> >
class ITK_EXPORT MaskedImageToImageFilter
  :public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MaskedImageToImageFilter                        Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MaskedImageToImageFilter, ImageToImageFilter);
  
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
  
  typedef TMaskImage                            MaskImageType;
  typedef typename MaskImageType::Pointer       MaskImagePointer;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,  TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  typedef ThreadLogger                          LoggerType;
  typedef typename LoggerType::Pointer          LoggerPointer;
  typedef std::vector<LoggerPointer>            LoggerVectorType;
  typedef utl_shared_ptr<LoggerVectorType>      LoggerVectorPointer;
  
  /** Set/Get the mask image. */
  itkSetObjectMacro(MaskImage, MaskImageType);
  itkGetObjectMacro(MaskImage, MaskImageType);
  itkGetConstObjectMacro(MaskImage, MaskImageType);
  void SetMaskImage(const std::string& file);
  
  /** Set/Get the mask image. */
  itkSetObjectMacro(Logger, LoggerType);
  itkGetObjectMacro(Logger, LoggerType);
  itkGetConstObjectMacro(Logger, LoggerType);

  itkSetGetMacro(LogLevel, int);

  bool IsMaskUsed()
    {
    return !IsImageEmpty(this->m_MaskImage);
    }


protected:
  MaskedImageToImageFilter();
  virtual ~MaskedImageToImageFilter();

  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;

  virtual void VerifyMaskInformation() const;

  virtual void VerifyInputParameters() const
    {
    if (!IsImageEmpty(this->m_MaskImage))
      this->VerifyMaskInformation();
    }
  
  /** Use single thread for MKL or openmp to avoid confliction with itk threader. 
   * Otherwise it has worse performance due to thread conflication (but correct solution).  */
  virtual void InitializeThreadedLibraries();
  
  /** create m_LoggerVector only when m_Debug is true and multiple threads are used  */
  void CreateLoggerVector ();
  
  void WriteLogger(const std::string& str, const LoggerBase::PriorityLevelType level=LoggerBase::DEBUG) const;

  std::string ThreadIDToString() const;

  MaskImagePointer m_MaskImage;

  LoggerPointer        m_Logger;
  LoggerVectorPointer  m_LoggerVector;

  /** It is used for ThreadLogger. It is -1 as default, and it is set as thread id when *this is cloned.  */
  int m_ThreadID;
  
  /** Debug verbose log level. Default value is 1 (LOG_NORMAL), which only shows some process information. 
   * It is normally set as utl::LogLevel, but each class could have its own logLevel. 
   * If it is LOG_MUTE, mute all output messeage. 
   * If it is LOG_LARGE, print some large matrix etc.*/
  int m_LogLevel;

private:
  MaskedImageToImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

};

}

// Define instantiation macro for this template.
#define ITK_TEMPLATE_MaskedImageToImageFilter(_, EXPORT, TypeX, TypeY)                  \
  namespace itk                                                                   \
  {                                                                               \
  _( 2 ( class EXPORT MaskedImageToImageFilter< ITK_TEMPLATE_2 TypeX > ) )              \
  namespace Templates                                                             \
  {                                                                               \
  typedef MaskedImageToImageFilter< ITK_TEMPLATE_2 TypeX > MaskedImageToImageFilter##TypeY; \
  }                                                                               \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkMaskedImageToImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkMaskedImageToImageFilter_hxx)
#include "itkMaskedImageToImageFilter.hxx"
#endif


#endif 

