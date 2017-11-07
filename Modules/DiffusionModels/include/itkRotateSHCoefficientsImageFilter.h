/**
 *       @file  itkRotateSHCoefficientsImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "12-30-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkRotateSHCoefficientsImageFilter_h
#define __itkRotateSHCoefficientsImageFilter_h

#include <vnl/vnl_matrix.h>
#include "utl.h"
#include "itkSpecialFunctionGenerator.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkSHCoefficientsRotation.h"


namespace itk
{
  
/** \class RotateSHCoefficientsImageFilter
 *
 * \brief Rotate SH coefficient to all input pixels.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \ingroup DiffusionModels
 * \sa UnaryFunctorImageFilter
 *
 * \author Jian Cheng
 */
namespace Functor {  
  
template< class TInput, class TOutput>
class RotateSHCoefficients
{
public:
  typedef utl::NDArray<double,2>   MatrixType;
  typedef utl::NDArray<double,1>   VectorType;
  typedef SHCoefficientsRotation<double>  SHRotationFilter;
  typedef typename SHRotationFilter::Pointer   SHRotationPointer;

  RotateSHCoefficients() : m_RotationMatrix(MatrixType(3,3))
    {
    m_RotationMatrix.SetIdentity();
    m_SHRotate = SHRotationFilter::New();
    };
  ~RotateSHCoefficients() 
    {
    };
  bool operator!=( const RotateSHCoefficients & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const RotateSHCoefficients & other ) const
    {
    return m_RotationMatrix.IsEqual(other.m_RotationMatrix, 1e-8);
    }
  inline TOutput operator()( const TInput & A ) const
    {
    const int N = A.GetSize();
    TOutput output(A);
    VectorType vec(N);
    utl::VectorToVector<TInput, VectorType >(A, vec, N);
    VectorType vecRotated = m_SHRotate->GetRotatedSHCoefficients(vec, m_RotationMatrix);
    utl::VectorToVector<VectorType, TInput>(vecRotated, output, N);
    return output;
    }

  void SetRotationMatrix(const MatrixType& mat) 
    { this->m_RotationMatrix = mat; }

  const MatrixType & GetRotationMatrix() const 
    { return m_RotationMatrix; }
  
  void SetSHRotate(const SHRotationPointer& ptr) 
    { this->m_SHRotate = ptr; }

  const SHRotationPointer & GetSHRotate() const 
    { return m_SHRotate; }
  
  MatrixType m_RotationMatrix;
  SHRotationPointer m_SHRotate;
};
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT RotateSHCoefficientsImageFilter :
public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
    Functor::RotateSHCoefficients< 
    typename TInputImage::PixelType,
    typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef RotateSHCoefficientsImageFilter             Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Functor::RotateSHCoefficients< 
      typename TInputImage::PixelType,
      typename TOutputImage::PixelType>   >             Superclass;

  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef utl::NDArray<double,2>   MatrixType;
  typedef utl::NDArray<double,1>   VectorType;
  typedef SHCoefficientsRotation<double>  RotationFilter;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RotateSHCoefficientsImageFilter, UnaryFunctorImageFilter);

  void BeforeThreadedGenerateData() ITK_OVERRIDE
    {
    this->GetFunctor().GetSHRotate()->SetTessOrder(3);
    this->GetFunctor().GetSHRotate()->SetMaxRank(10);
    this->GetFunctor().GetSHRotate()->Initialize();
    }
  
  /** Set the sigma of Rician noise distribution */
  void SetRotationMatrix(const MatrixType& mat)
    {
    if( !mat.IsEqual(this->GetFunctor().GetRotationMatrix(), 1e-8) )
      {
      this->GetFunctor().SetRotationMatrix(mat);
      this->Modified();
      }
    }
  const MatrixType & GetRotationMatrix() const
    {
    return this->GetFunctor().GetRotationMatrix();
    }

protected:
  RotateSHCoefficientsImageFilter() { }
  virtual ~RotateSHCoefficientsImageFilter() {}
   
  void PrintSelf(std::ostream &os, Indent indent) const ITK_OVERRIDE
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "RotationMatrix: " 
       << this->GetRotationMatrix()
       << std::endl;
    }


private:
  RotateSHCoefficientsImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


} // end namespace itk


#endif 

