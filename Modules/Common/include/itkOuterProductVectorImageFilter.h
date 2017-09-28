/**
 *       @file  itkOuterProductVectorImageFilter.h
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#ifndef __itkOuterProductVectorImageFilter_h
#define __itkOuterProductVectorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "utlNDArray.h"

namespace itk
{
/** \class OuterProductVectorImageFilter
 *
 * \brief In each vxoel, multiply the input vector by a matrix.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 */
namespace Functor
{
template <class TInput, class TOutput=TInput>
class OuterProduct
{
public:
  typedef typename TInput::ValueType   InputValueType;
  typedef typename TOutput::ValueType  OutputValueType;
  typedef utl::NDArray<double, 1>      VectorType;

  OuterProduct()
  {
  }

  ~OuterProduct()
  {
  }

  bool operator!=(const OuterProduct & other) const
  {
    return !( *this == other );
  }

  bool operator==(const OuterProduct & other) const
  {
    return true;
  }

  inline TOutput operator()(const TInput & A)
  {
  int N = A.GetSize();
  TOutput out;
  out.SetSize( N*N );
  int kk=0;
  for ( int i = 0; i < N; ++i ) 
    {
    for ( int j = 0; j < N; ++j ) 
      {
      out[kk] = A[i]*A[j];
      kk++;
      }
    }

  return out;
  }

};
}

template <class TInputImage, class TOutputImage=TInputImage>
class ITK_EXPORT OuterProductVectorImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::OuterProduct<
                            typename TInputImage::PixelType,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef OuterProductVectorImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Functor::OuterProduct<
        typename TInputImage::PixelType, typename TOutputImage::PixelType> >             Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( OuterProductVectorImageFilter, UnaryFunctorImageFilter );

protected:
  OuterProductVectorImageFilter()
  {
  }

  virtual ~OuterProductVectorImageFilter()
  {
  }

  virtual void GenerateOutputInformation()
    {
    Superclass::GenerateOutputInformation();
    typename TInputImage::ConstPointer inputPtr = this->GetInput();
    typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
    int N = inputPtr->GetNumberOfComponentsPerPixel();
    outputPtr->SetNumberOfComponentsPerPixel(N*N);
    }


private:
  OuterProductVectorImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};
} // end namespace itk


#endif 
