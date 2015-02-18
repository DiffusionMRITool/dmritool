/**
 *       @file  itkMultiplyByConstantMatrixVectorImageFilter.h
 *      @brief  
 *     Created  "02-16-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkMultiplyByConstantMatrixVectorImageFilter_h
#define __itkMultiplyByConstantMatrixVectorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "utlNDArray.h"

namespace itk
{
/** \class MultiplyByConstantMatrixVectorImageFilter
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
template <class TInput, class TMatrix, class TOutput=TInput>
class MultiplyByConstantMatrix
{
public:
  typedef typename TInput::ValueType   InputValueType;
  typedef typename TOutput::ValueType  OutputValueType;
  typedef utl::NDArray<double, 1>      VectorType;

  MultiplyByConstantMatrix()
  {
  }

  ~MultiplyByConstantMatrix()
  {
  }

  bool operator!=(const MultiplyByConstantMatrix & other) const
  {
    return !( *this == other );
  }

  bool operator==(const MultiplyByConstantMatrix & other) const
  {
    return other.m_ConstantMatrix == m_ConstantMatrix;
  }

  inline TOutput operator()(const TInput & A)
  {
    // Because the user has to specify the constant we don't
    // check if the cte is not 0;
  VectorType  vectorInput(A.GetSize());
  VectorType  vectorOutput(m_ConstantMatrix.Rows());
  TOutput out;
  out.SetSize( m_ConstantMatrix.Rows() );

  utl::VectorToVector<TInput, VectorType>(A, vectorInput, A.Size());
  utl::ProductUtlMv(m_ConstantMatrix, vectorInput, vectorOutput );
  utl::VectorToVector<VectorType, TOutput>(vectorOutput, out, vectorOutput.Size());

  return out;
  }

  void SetConstantMatrix(const TMatrix& ct)
  {
    this->m_ConstantMatrix = ct;
  }

  const TMatrix & GetConstantMatrix() const
  {
    return m_ConstantMatrix;
  }

  const unsigned GetNumberOfComponentsPerPixelInOutput() const
    {
    return m_ConstantMatrix.Rows();
    }

  TMatrix m_ConstantMatrix;
};
}

template <class TInputImage, class TConstantMatrix, class TOutputImage=TInputImage>
class ITK_EXPORT MultiplyByConstantMatrixVectorImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::MultiplyByConstantMatrix<
                            typename TInputImage::PixelType, TConstantMatrix,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef MultiplyByConstantMatrixVectorImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Functor::MultiplyByConstantMatrix<
        typename TInputImage::PixelType, TConstantMatrix,
        typename TOutputImage::PixelType> >             Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiplyByConstantMatrixVectorImageFilter, UnaryFunctorImageFilter );

  /** Set the constant that will be used to multiply all the image pixels */
  void SetConstantMatrix(const TConstantMatrix& ct )
  {
    if( ct != this->GetFunctor().GetConstantMatrix() )
      {
      this->GetFunctor().SetConstantMatrix( ct );
      this->Modified();
      }
  }

  const TConstantMatrix & GetConstantMatrix() const
  {
    return this->GetFunctor().GetConstantMatrix();
  }

protected:
  MultiplyByConstantMatrixVectorImageFilter()
  {
  }

  virtual ~MultiplyByConstantMatrixVectorImageFilter()
  {
  }

  virtual void GenerateOutputInformation()
    {
    Superclass::GenerateOutputInformation();
    typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
    unsigned dim = this->GetFunctor().GetNumberOfComponentsPerPixelInOutput();
    outputPtr->SetNumberOfComponentsPerPixel(dim);
    }

  void PrintSelf(std::ostream & os, Indent indent ) const
  {
    Superclass::PrintSelf( os, indent );
    os << indent << "Constant Matrix : " << this->GetConstantMatrix() << std::endl;
  }

private:
  MultiplyByConstantMatrixVectorImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};
} // end namespace itk


#endif 
