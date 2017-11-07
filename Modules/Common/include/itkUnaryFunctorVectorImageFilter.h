/**
 *       @file  itkUnaryFunctorVectorImageFilter.h
 *      @brief  
 *     Created  "09-11-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef itkUnaryFunctorVectorImageFilter_h
#define itkUnaryFunctorVectorImageFilter_h

#include "itkFunctorBaseVectorImageFilter.h"

namespace itk
{
/** \class UnaryFunctorVectorImageFilter
 * \brief Implements vector-valued generic operation on one image.
 *
 *  Calculation loops over image domain and extracts vectors along x/y/z/t-axis. 
 *
 *  m_Functor is a functor with utl::Vector as input and utl::Vector as output. 
 *  m_Functor needs to define GetOutputDimension to obtain the NumberOfComponentsPerPixel in output image. 
 *
 * \author Jian Cheng
 *
 * \ingroup ITKCommon
 *
 */
template< typename TInputImage, typename TOutputImage, typename TFunction, class TMaskImage=Image<double,3> >
class ITK_EXPORT UnaryFunctorVectorImageFilter
  :public FunctorBaseVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage >
{
public:
  /** Standard class typedefs. */
  typedef UnaryFunctorVectorImageFilter                         Self;
  typedef FunctorBaseVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(UnaryFunctorVectorImageFilter, FunctorBaseVectorImageFilter);

  /** Some typedefs. */
  typedef TFunction FunctorType;

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
  
  typedef typename Superclass::MaskImageType      MaskImageType;

protected:
  UnaryFunctorVectorImageFilter();
  virtual ~UnaryFunctorVectorImageFilter() {}

  void BeforeThreadedGenerateData () ITK_OVERRIDE;

  /** UnaryFunctorVectorImageFilter can produce an image which is a different
   * resolution than its input image.  As such, UnaryFunctorVectorImageFilter
   * needs to provide an implementation for
   * GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below.
   *
   * \sa ProcessObject::GenerateOutputInformaton()  */
  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  /** UnaryFunctorVectorImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData() routine
   * which is called for each processing thread. The output image data is
   * allocated automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to the
   * portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId) ITK_OVERRIDE;

protected: 
  

private:
  UnaryFunctorVectorImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};


/** 
 * By defining a functor which maps a vector to another vector, 
 * this function perform the functor on an image (VectorImage or NDImage, along axis). 
 * */
template <class ImageType, class ImageOutType, class OpFunctor, class MaskImageType=Image<double,4> >
void
UnaryVectorOPImage(const itk::SmartPointer<ImageType>& image, itk::SmartPointer<ImageOutType>& outImage, const OpFunctor& func, const itk::SmartPointer<MaskImageType>& mask=nullptr, int numberOfThreads=-1, int vectorAxis=3)
{
  typedef itk::UnaryFunctorVectorImageFilter<ImageType, ImageOutType, OpFunctor, MaskImageType> UnaryFunctorFilterType;
  typename UnaryFunctorFilterType::Pointer filter = UnaryFunctorFilterType::New();

  if (!itk::IsImageEmpty(mask))
    filter->SetMaskImage(mask);
  filter->SetInput(image);
  if (vectorAxis>=0)
    filter->SetVectorAxis(vectorAxis);
  filter->SetFunctor(func);
  filter->SetDebug(utl::IsLogDebug());
  filter->SetLogLevel(utl::LogLevel);
  if (numberOfThreads>0)
    filter->SetNumberOfThreads(numberOfThreads);

  filter->Update();

  outImage = filter->GetOutput();       
}

} // end namespace itk

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkUnaryFunctorVectorImageFilter_hxx)
#include "itkUnaryFunctorVectorImageFilter.hxx"
#endif

#endif
