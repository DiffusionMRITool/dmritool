/**
 *       @file  itkMultiVariableFunctorVectorImageFilter.h
 *      @brief  
 *     Created  "09-13-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef itkMultiVariableFunctorVectorImageFilter_h
#define itkMultiVariableFunctorVectorImageFilter_h

#include "itkFunctorBaseVectorImageFilter.h"

namespace itk
{
/** \class MultiVariableFunctorVectorImageFilter
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
class ITK_EXPORT MultiVariableFunctorVectorImageFilter
  :public FunctorBaseVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiVariableFunctorVectorImageFilter                         Self;
  typedef FunctorBaseVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiVariableFunctorVectorImageFilter, FunctorBaseVectorImageFilter);

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
  MultiVariableFunctorVectorImageFilter();
  virtual ~MultiVariableFunctorVectorImageFilter() {}

  void BeforeThreadedGenerateData () ITK_OVERRIDE;

  /** Override Superclass::PropagateRequestedRegion. Otherwise, it may have region outside error.  */
  virtual void PropagateRequestedRegion(DataObject *output) ITK_OVERRIDE
    {
    }

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId) ITK_OVERRIDE;

protected: 
  

private:
  MultiVariableFunctorVectorImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};


/** 
 * By defining a functor which maps multiple vectors to another vector, 
 * this function perform the functor on a set of images (VectorImage or NDImage, along axis). 
 * */
template <class ImageType, class ImageOutType, class OpFunctor, class MaskImageType=Image<double,4> >
void
MultiVariableVectorOPImage(const std::vector<itk::SmartPointer<ImageType> >& images, itk::SmartPointer<ImageOutType>& outImage, const OpFunctor& func, const itk::SmartPointer<MaskImageType>& mask=nullptr, int numberOfThreads=-1, int vectorAxis=3)
{
  typedef itk::MultiVariableFunctorVectorImageFilter<ImageType, ImageOutType, OpFunctor, MaskImageType> FunctorImageFilterType;
  typename FunctorImageFilterType::Pointer filter = FunctorImageFilterType::New();

  if (!itk::IsImageEmpty(mask))
    filter->SetMaskImage(mask);
  for ( int i = 0; i < images.size(); ++i ) 
    {
    filter->SetInput(i, images[i]);
    }
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

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkMultiVariableFunctorVectorImageFilter_hxx)
#include "itkMultiVariableFunctorVectorImageFilter.hxx"
#endif

#endif
