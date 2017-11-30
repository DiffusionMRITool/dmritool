/**
 *       @file  itkFunctorFromStringImageFilter.h
 *      @brief  
 *     Created  "04-06-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#ifndef itkFunctorFromStringImageFilter_h
#define itkFunctorFromStringImageFilter_h

#include "itkFunctorBaseVectorImageFilter.h"
#include "utlExprtk.h"
#include "utlITKMacro.h"
#include "utlFunctors.h"

namespace itk
{
/** \class FunctorFromStringImageFilter
 * \brief Implements vector-valued generic operation on images with the same size.
 *
 *  m_Functor is not used. 
 *  m_Expression is a math expression maps one scalar value (one input) or several sacalr values (several inputs) to another scalar (output image). 
 *  This filter performs m_Expression element-wise on the images.
 *
 * \author Jian Cheng
 *
 * \ingroup ITKCommon
 *
 */
template< typename TInputImage, typename TOutputImage, class TMaskImage=Image<double,3> >
class ITK_EXPORT FunctorFromStringImageFilter
  :public FunctorBaseVectorImageFilter< TInputImage, TOutputImage, utl::Functor::VectorMultiVariableFunctionWrapper<>, TMaskImage >
{
public:
  /** Standard class typedefs. */
  typedef FunctorFromStringImageFilter                         Self;
  typedef FunctorBaseVectorImageFilter< TInputImage, TOutputImage, utl::Functor::VectorMultiVariableFunctionWrapper<>, TMaskImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FunctorFromStringImageFilter, FunctorBaseVectorImageFilter);

  /** Some typedefs. */
  typedef utl::Functor::VectorMultiVariableFunctionWrapper<> FunctorType;

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

  itkSetGetMacro(Expression, std::string);

protected:
  FunctorFromStringImageFilter();
  virtual ~FunctorFromStringImageFilter() {}

  void BeforeThreadedGenerateData () ITK_OVERRIDE;

  /** Override Superclass::PropagateRequestedRegion. Otherwise, it may have region outside error.  */
  virtual void PropagateRequestedRegion(DataObject *output) ITK_OVERRIDE
    {
    }

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId) ITK_OVERRIDE;

protected: 

  std::string m_Expression;

private:
  FunctorFromStringImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};


/** 
 * By defining a functor which maps multiple vectors to another vector, 
 * this function perform the functor on a set of images (VectorImage or NDImage). 
 * */
template <class ImageType, class ImageOutType, class MaskImageType=Image<double,4> >
void
FunctorFromStringOPImage(const std::vector<itk::SmartPointer<ImageType> >& images, itk::SmartPointer<ImageOutType>& outImage, const std::string& funcStr, const itk::SmartPointer<MaskImageType>& mask=nullptr, int numberOfThreads=-1)
{
  typedef itk::FunctorFromStringImageFilter<ImageType, ImageOutType, MaskImageType> FunctorImageFilterType;
  typename FunctorImageFilterType::Pointer filter = FunctorImageFilterType::New();

  if (!itk::IsImageEmpty(mask))
    filter->SetMaskImage(mask);
  for ( int i = 0; i < images.size(); ++i ) 
    {
    filter->SetInput(i, images[i]);
    }
  filter->SetExpression(funcStr);
  filter->SetDebug(utl::IsLogDebug());
  filter->SetLogLevel(utl::LogLevel);
  if (numberOfThreads>0)
    filter->SetNumberOfThreads(numberOfThreads);

  filter->Update();

  outImage = filter->GetOutput();       
}


} // end namespace itk

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkFunctorFromStringImageFilter_hxx)
#include "itkFunctorFromStringImageFilter.hxx"
#endif

#endif
