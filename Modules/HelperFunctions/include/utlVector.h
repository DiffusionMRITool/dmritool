/**
 *       @file  utlVector.h
 *      @brief  utl::NDArray<T,1> class which uses blas mkl
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlVector_h
#define __utlVector_h

#include "utlCore.h"
#include "utlBlas.h"
// #include "utlCore.h"
#include <vector>
#include "utlSTDHeaders.h"
#include "utlExpression.h"
#include "utlNDArray.h"

namespace utl
{
/** \ingroup utlNDArray
* @{ */

// typedef NDArray<double, 1>  VectorDouble;
// typedef NDArray<float, 1>  VectorFloat;
// template <class T> class Matrix;

template <class T, class EType>
utl_shared_ptr< NDArray<T,1> >
ToVector ( const Expr<EType>& expr )
{
  utl_shared_ptr< NDArray<T,1> > vec (new NDArray<T,1>(expr));
  return vec;
}

template <class T, class EType>
utl_shared_ptr< NDArray<T,2> >
ToMatrix ( const Expr<EType>& expr )
{
  utl_shared_ptr< NDArray<T,2> > mat(new NDArray<T,2>(expr));
  return mat;
}

template <class T>
std::vector<T> 
UtlVectorToStdVector ( const NDArray<T,1>& vec )
{
  std::vector<T> v(vec.Size());
  for ( int i = 0; i < vec.Size(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
NDArray<T,1> 
StdVectorToUtlVector ( const std::vector<T>& vec )
{
  utl::NDArray<T,1> v(vec.size());
  for ( int i = 0; i < vec.size(); ++i ) 
    v[i] = vec[i];
  return v;
}

template <class T>
NDArray<T,1>
ConnectUtlVector ( const NDArray<T,1>& m1, const NDArray<T,1>& m2 )
{
  NDArray<T,1> result(m1.Size()+m2.Size());
  int m1Size = m1.Size();
  for ( int i = 0; i < m1Size; i += 1 ) 
    result[i] = m1[i];
  for ( int i = 0; i < m2.Size(); i += 1 ) 
    result[i+m1Size] = m2[i];
  return result;
}

/** http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ Y = alpha * A * X + beta * Y \f$
*  \note: Y should be pre-allocated
*
*  define several functions. 
*
*  template \<class T\>  inline bool  
*  gemv_UtlMatrixTimesVector(const bool bATrans, const T alpha, const utl::NDArray<T,2>& A, const utl::NDArray<T,1>& X, const T beta, utl::NDArray<T,1>& Y);
*
*  \f$ Y = alpha * A * X + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlMv(const utl::NDArray<T,2>& A, const utl::NDArray<T,1>& b, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
*  \f$ Y = alpha * A^T * X + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlMtv(const utl::NDArray<T,2>& A, const utl::NDArray<T,1>& b, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
* */
__utl_gemv_MatrixTimesVector(T, gemv_UtlMatrixTimesVector, Utl, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, utl::NDArray<T UTL_COMMA 1>, Size, GetData, ReSize);


/**
*  \f$ Y = alpha X^T * A + beta * Y \f$
*  \note: Y should be pre-allocated
*
*  define several functions. 
*
*  template \<class T\>  inline bool  
*  gemm_UtlVectorTimesMatrix(const bool bATrans, const T alpha, const utl::NDArray<T,1>& X, const utl::NDArray<T,2>& A, const T beta, utl::NDArray<T,1>& Y);
*
*  \f$ Y = alpha * b * A + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlvM(const utl::NDArray<T,1>& b, const utl::NDArray<T,2>& A, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
*  \f$ Y = alpha * b * A^T + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlvMt(const utl::NDArray<T,1>& b, const utl::NDArray<T,2>& A, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
* */
__utl_gevm_MatrixTimesVector(T, gemm_UtlVectorTimesMatrix, Utl, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, utl::NDArray<T UTL_COMMA 1>, Size, GetData, ReSize);
 


/**
 *   \class   NDArray
 *   \brief   NDArray<T,1> is a vector class which uses blas mkl
 *   \author  Jian Cheng
 *   \date    08-22-2014
 *   \ingroup utlNDArray Math
 */
template < class T >
class NDArray<T, 1> : public NDArrayBase<T,1>
{

public:
  
  typedef NDArray          Self;
  typedef NDArrayBase<T,1>    Superclass;

  typedef typename Superclass::ValueType          ValueType;
  
  typedef typename Superclass::SizeType           SizeType;
  typedef typename Superclass::ShapeType          ShapeType;
  
  typedef typename Superclass::Iterator           Iterator;
  typedef typename Superclass::ConstIterator      ConstIterator;
  typedef typename Superclass::Pointer            Pointer;
  typedef typename Superclass::ConstPointer       ConstPointer;
  typedef typename Superclass::Reference          Reference;
  typedef typename Superclass::ConstReference     ConstReference;

  using Superclass::Dimension;

  // NOTE: "using" to avoid hiding member functions of superclass
  using Superclass::SetData;
  using Superclass::CopyData;
  using Superclass::ReSize;
  using Superclass::operator();
  using Superclass::operator=;
  using Superclass::operator+=;
  using Superclass::operator-=;
  using Superclass::operator%=;
  using Superclass::operator/=;

public:
  /**
   * Default constructor uses compiler's default initialization of memory.
   * For efficiency, no initialization to zero is done.
   */
  NDArray() : Superclass()
  {}

#define __NDArray_vector_constructor(intType)             \
  explicit NDArray(const intType size) : Superclass()     \
  {                                                       \
  SizeType shape[1];                                      \
  shape[0]=size;                                          \
  __utl_ndarray_alloc_blah(shape);                        \
  }

__NDArray_vector_constructor(unsigned int)
__NDArray_vector_constructor(int)
__NDArray_vector_constructor(std::size_t)
__NDArray_vector_constructor(long)

#undef __NDArray_vector_constructor

  NDArray(const NDArray<T,1>& vec) : Superclass(vec)
  {
  }
  
  template<typename EType>
  NDArray(const Expr<EType>& expr) : Superclass(expr)
  {
  }

/**
 * Constructor assumes input points to array of correct size.
 * Values are copied individually instead of with a binary copy.  This
 * allows the T's assignment operator to be executed.
 */
  NDArray(const T* vec, const SizeType size) : Superclass()
  {
  SizeType shape[1];
  shape[0]=size;
  __utl_ndarray_alloc_blah(shape);
  utl::cblas_copy<T>(this->Size(), vec, 1, this->m_Data, 1);
  }

  /**
   * Constructor to initialize entire array to one value.
   */
  NDArray(const SizeType size, const T r) : Superclass()
  {
  SizeType shape[1];
  shape[0]=size;
  __utl_ndarray_alloc_blah(shape);
  this->Fill(r);
  }

  template< typename TVectorValueType >
  NDArray(const NDArray< TVectorValueType, 1> & r) : Superclass(r)
  {
  }

  explicit NDArray(const ShapeType& shape) : Superclass(shape) { }

  NDArray(const T* vec, const ShapeType& shape) : Superclass(vec,shape)   {  }

  /**
   * Constructor to initialize entire array to one value.
   */
  NDArray(const ShapeType& shape, const T r) : Superclass(shape, r)   {  }


  /** set m_Data from data, the ownership is in data. */
  void SetData( T* const data, const SizeType size )
    {
    SizeType shape[1];
    shape[0]=size;
    Superclass::SetData(data, shape);
    }

  /** copy m_Data from data  */
  void CopyData(T* const data, const SizeType size)
    {
    SizeType shape[1];
    shape[0]=size;
    Superclass::CopyData(data, shape);
    }

  inline bool ReSize(const SizeType size)
    {
    SizeType shape[1];
    shape[0]=size;
    return Superclass::ReSize(shape);
    }

// #define __Vector_Saver_Expr(saver)                                                     \
// template<typename EType>                                                               \
// inline NDArray<T,1>& operator saver (const Expr<EType>& src){                             \
//   Superclass::operator saver (src);                                                    \
//   return *this;                                                                        \
// }                                                                                      
  
//   // __Vector_Saver_Expr(=)
//   __Vector_Saver_Expr(+=)
//   __Vector_Saver_Expr(-=)
//   __Vector_Saver_Expr(%=)
//   __Vector_Saver_Expr(/=)

// #define __Vector_Saver_Scalar(saver)                    \
// inline NDArray<T,1>& operator saver (const T val){         \
//   Superclass::operator saver(val);                      \
//   return *this;                                         \
// }                                                       

//   // __Vector_Saver_Scalar(=)
//   __Vector_Saver_Scalar(+=)
//   __Vector_Saver_Scalar(-=)
//   __Vector_Saver_Scalar(%=)
//   __Vector_Saver_Scalar(/=)
  
  /** \note The separate versions are a work-around to avoid ambiguity for operator()(const ShapeType& shape), 
   * when converting an integer number is needed.   */
  UTL_ALWAYS_INLINE Reference operator()(unsigned long index)         { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator()(unsigned long index) const { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator()(long index)         { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator()(long index) const { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator()(int index)         { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator()(int index) const { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator()(unsigned int index)         { return this->m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator()(unsigned int index) const { return this->m_Data[index]; }
  
  // UTL_ALWAYS_INLINE Reference operator()(const ShapeType& shape)         { return this->m_Data[shape[0]]; }
  // UTL_ALWAYS_INLINE ConstReference operator()(const ShapeType& shape) const { return this->m_Data[shape[0]]; }


  NDArray<T,1>& operator*=(const NDArray<T,2>& mat) 
    {
    utlSAException(mat.Rows()!=this->Size())(this->Size())(mat.Rows()).msg("wrong size");
    NDArray<T,1> tmp;
    utl::ProductUtlvM(*this, mat, tmp);
    *this = tmp;
    return *this;
    }
  // NDArray<T,1> operator*(const NDArray<T,2>& mat) const
  //   {
  //   utlSAException(mat.Rows()!=this->Size())(this->Size())(mat.Rows()).msg("wrong size");
  //   NDArray<T,1> tmp;
  //   utl::ProductUtlvM(*this, mat, tmp);
  //   return tmp;
  //   }
  
  template< typename EType >
  NDArray<T,1>& operator*=(const Expr<EType>& expr ) 
    {
    utl_shared_ptr<NDArray<T,2> > mat = ToMatrix(expr);
    operator*=(mat);
    return *this;
    }
  template< typename EType >
  NDArray<T,1> operator*(const Expr<EType>& expr ) const
    {
    utl_shared_ptr<NDArray<T,2> > mat = ToMatrix(expr);
    return operator*(mat);
    }


  /** vec = mat * (*this), considering *this is a column vector.  */
  void PreMultiply(const NDArray<T,2>& mat, NDArray<T,1>& vec)
    {
    utl::ProductUtlMv(mat, *this, vec);
    }
  
  void PostMultiply(const NDArray<T,2>& mat, NDArray<T,1>& vec)
    {
    utl::ProductUtlvM(*this, mat, vec);
    }

  void GetDiagonalMatrix(NDArray<T,2>& mat) const
    {
    mat.ReSize(this->Size(), this->Size());
    mat.Fill(0.0);
    mat.SetDiagonal(*this);
    }
  NDArray<T,2> GetDiagonalMatrix() const
    {
    NDArray<T,2> result(this->Size(), this->Size());
    this->GetDiagonalMatrix(result);
    return result;
    }

  void GetSubVector(const std::vector<int>& indexVec, NDArray<T,1>& vec) const
    {
    vec.ReSize(indexVec.size());
    for ( int i = 0; i < indexVec.size(); ++i ) 
      vec[i] = this->m_Data[ indexVec[i] ];
    }
  NDArray<T,1> GetSubVector(const std::vector<int>& indexVec) const
    {
    NDArray<T,1> vec;
    GetSubVector(indexVec, vec);
    return vec;
    }

  NDArray<T,1>& SetSubVector(const std::vector<int> indexVec, const NDArray<T,1>& vec)
    {
    for ( int i = 0; i < indexVec.size(); ++i ) 
      this->m_Data[ indexVec[i] ] = vec[i];
    return *this;
    }


protected:


private:


}; // -----  end of template class Vector  -----



template< typename T >
inline NDArray<T,1> 
operator*(const NDArray<T,2>& mat, const NDArray<T,1>& vec) 
{ 
  utlSAException(vec.Size()!=mat.Cols())(mat.Cols())(vec.Size()).msg("wrong size");
  NDArray<T,1> result;
  utl::ProductUtlMv(mat, vec, result);
  return result;
}

template< typename T >
inline NDArray<T,1> 
operator*(const NDArray<T,1>& vec, const NDArray<T,2>& mat) 
{ 
  utlSAException(vec.Size()!=mat.Rows())(mat.Rows())(vec.Size()).msg("wrong size");
  NDArray<T,1> result;
  utl::ProductUtlvM(vec, mat, result);
  return result;
}

template< typename T >
std::ostream & 
operator<<(std::ostream & os, const NDArray<T,1> & arr)
{
  utl::PrintContainer(arr.Begin(), arr.End(), "utl::NDArray<T,1>", " ", os);
  return os;
}



/** @}  */

}


#endif 
