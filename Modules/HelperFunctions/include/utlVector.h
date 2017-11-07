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

#include "utlNDArray.h"

namespace utl
{
/** \ingroup utlNDArray
* @{ */


template < class T, unsigned int Dim >
class NDArray;

template < class T, unsigned int Dim >
class NDArrayBase;

/**
 *   \class   NDArray<T,1>
 *   \brief   NDArray<T,1> is a vector class which uses blas mkl
 *
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
  typedef typename Superclass::ScalarValueType    ScalarValueType;
  
  typedef typename Superclass::SizeType           SizeType;
  typedef typename Superclass::ShapeType          ShapeType;
  
  typedef typename Superclass::Iterator           Iterator;
  typedef typename Superclass::ConstIterator      ConstIterator;
  typedef typename Superclass::Pointer            Pointer;
  typedef typename Superclass::ConstPointer       ConstPointer;
  typedef typename Superclass::Reference          Reference;
  typedef typename Superclass::ConstReference     ConstReference;

  using Superclass::Dimension;
  using Superclass::SubDimension;

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

  NDArray(const NDArray<T,1>& vec) : Superclass(vec) { }
  
  NDArray(NDArray<T,1>&& vec)   
  { 
    operator=(std::move(vec)); 
  }
  
  NDArray(const std::initializer_list<T>& r)
    {
    ReSize(r.size());
    utl::cblas_copy(this->Size(), r.begin(), 1, this->m_Data, 1);
    }
  
  template<typename EType>
  NDArray(const Expr<EType, typename EType::ValueType>& expr) : Superclass(expr)
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

  
  NDArray<T,1>& operator=(NDArray<T,1> & r)
    {
    Superclass::operator=(r);
    return *this;
    }

  NDArray<T,1>& operator=(NDArray<T,1> && r)
    {
    if ( this != &r ) 
      {
      this->Clear();
      this->Swap(r);
      }
    return *this;
    }


  /** rotate *this about axis by angle alpha anticlockwise in 3d space. 
   * https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula */
  NDArray<T,1> GetRotateAxis(const double alpha, const NDArray<T, 1>& axis)
    {
    utlException(this->Size()!=3, "wrong size");
    utlException(axis.Size()!=3, "wrong size");
    double cosA = std::cos(alpha);
    double sinA = std::sin(alpha);
    NDArray<T,1> p1 = *this % cosA;
    NDArray<T,1> p2(3);
    CrossProduct(axis, *this, p2);
    p2.Scale(sinA);
    NDArray<T,1> p3 = axis % (DotProduct(axis, *this, 3)*(1-cosA));
    return p1+p2+p3;
    }

  void RotateAxis(const double alpha, const NDArray<T, 1>& axis)
    {
    *this = GetRotateAxis(alpha, axis);
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




/** @}  */

}


#endif 
