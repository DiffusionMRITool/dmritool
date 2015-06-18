/**
 *       @file  utlNDArray.h
 *      @brief  
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlNDArray_h
#define __utlNDArray_h

#include <numeric>

#include "utlBlas.h"
#include "utlCore.h"
#include "utlExpression.h"
#include "utlCoreMacro.h"  
#include "utlSTDHeaders.h"

namespace utl
{
/** \ingroup utlNDArray
* @{ */

#define __utl_ndarray_alloc_blah(shape)                                                   \
do {                                                                                      \
  for ( int i = 0; i < this->Dimension; ++i )                                             \
    this->m_Shape[i] = shape[i];                                                          \
  this->ComputeOffSetTable();                                                             \
  SizeType size = this->GetSize();                                                        \
  this->m_Data = (size>0) ? (new T[size]) : NULL;                                         \
  this->m_IsShared = false;                                                               \
} while (false)                                                                           

template <class T, unsigned int Dim> class NDArray;
template <class T, unsigned int Dim> class NDArrayBase;

template <class T, unsigned int Dim, class EType>
utl_shared_ptr< NDArray<T,Dim> >
ToNDArray ( const Expr<EType>& expr )
{
  utl_shared_ptr< NDArray<T,Dim> > mat(new NDArray<T,Dim>(expr) );
  return mat;
}

template <class T>
inline void
PrintUtlMatrix(const NDArray<T,2>& mat, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  utl::PrintMatrix<NDArray<T,2> >(mat, mat.Rows(), mat.Columns(), str, separate, os);
}

template <class T>
inline void
PrintUtlVector(const NDArray<T,1>& vec, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  PrintContainer(vec.Begin(), vec.End(), str, separate, os);
}

template <class T, unsigned int Dim>
inline void
PrintUtlNDArray(const NDArrayBase<T,Dim>& arr, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  arr.Print(os<< str, separate );
}

template <class T, unsigned int Dim>
inline T
InnerProduct ( const NDArrayBase<T,Dim>& v1, const NDArrayBase<T,Dim>& v2 )
{
  utlSAException(v1.Size() != v2.Size())(v1.Size())(v2.Size()).msg("vector sizes mismatch");
  return utl::cblas_dot<T>(v1.Size(), v1.GetData(), 1, v2.GetData(), 1);
}

template <class T, unsigned int Dim>
inline T
DotProduct ( const NDArray<T,Dim>& v1, const NDArray<T,Dim>& v2 )
{
  return InnerProduct(v1, v2);
}

template< typename T, unsigned int Dim >
std::ostream & 
operator<<(std::ostream & os, const NDArray<T,Dim> & arr)
{
  arr.Print(os<< "utl::NDArray<T,Dim>" );
  return os;
}

/**
 *   \class   NDArray
 *   \brief   NDArray is a N-Dimensional array class (row-major, c version)
 *
 *   It uses blas, mkl, and expression templates.
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup utlNDArray
 */
template < class T, unsigned int Dim  >
class NDArray : public NDArrayBase<T,Dim>
{
public:
  typedef NDArray Self;
  typedef NDArrayBase<T,Dim>     Superclass;

  
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
  
  using Superclass::operator=;

public:
  NDArray() : Superclass()  {}

  explicit NDArray(const ShapeType& shape) : Superclass(shape) { }

  NDArray(const NDArray<T,Dim>& vec) : Superclass(vec) {  }
  
  template<typename EType>
  NDArray(const Expr<EType>& expr) : Superclass(expr)   {  }

/**
 * Constructor assumes input points to array of correct size.
 * Values are copied individually instead of with a binary copy.  This
 * allows the T's assignment operator to be executed.
 */
  NDArray(const T* vec, const ShapeType& shape) : Superclass(vec,shape)   {  }

  /**
   * Constructor to initialize entire array to one value.
   */
  NDArray(const ShapeType& shape, const T r) : Superclass(shape, r)   {  }

  template< typename TValue >
  NDArray(const NDArray<TValue, Dim> & r) : Superclass(r)  {  }
  
};


/**
 *   \class   NDArrayBase
 *   \brief   Base class for utl::NDArray
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup utlNDArray 
 */
template < class T, unsigned int Dim  >
class NDArrayBase : public Expr<NDArrayBase<T, Dim> > 
{
public:

  /** The element type stored at each location in the NDArrayBase. */
  typedef NDArrayBase Self;
  typedef Expr<NDArrayBase<T, Dim> >     Superclass;

  typedef typename Superclass::SizeType          SizeType;
  typedef typename Superclass::ShapeType         ShapeType;

  // typedef SizeType  ShapeType[Dim];

  enum { Dimension = Dim };
  
  typedef T ValueType;

  /** An iterator through the array. */
  typedef ValueType* Iterator;

  /** A const iterator through the array. */
  typedef const ValueType* ConstIterator;
  
  /** A pointer to the ValueType. */
  typedef ValueType* Pointer;
  typedef const ValueType* ConstPointer;

  /** A reference to the ValueType. */
  typedef ValueType & Reference;
  typedef const ValueType & ConstReference;

  class ConstReverseIterator;

  /** \class ReverseIterator
   * \brief A reverse iterator through an array.
   */
  class ReverseIterator
  {
  public:
    explicit ReverseIterator(Iterator i):m_Iterator(i) {}
    Iterator operator++()        { return --m_Iterator; }
    Iterator operator++(int)     { return m_Iterator--; }
    Iterator operator--()        { return ++m_Iterator; }
    Iterator operator--(int)     { return m_Iterator++; }
    Iterator operator->() const { return ( m_Iterator - 1 ); }
    ValueType & operator*() const { return *( m_Iterator - 1 ); }
    bool operator!=(const ReverseIterator & rit) const { return m_Iterator != rit.m_Iterator; }
    bool operator==(const ReverseIterator & rit) const { return m_Iterator == rit.m_Iterator; }

  private:
    Iterator m_Iterator;
    friend class ConstReverseIterator;
  };

  /** \class ConstReverseIterator
   * \brief A const reverse iterator through an array.
   */
  class ConstReverseIterator
  {
  public:
    explicit ConstReverseIterator(ConstIterator i):m_Iterator(i) {}
    ConstReverseIterator(const ReverseIterator & rit) { m_Iterator = rit.m_Iterator; }
    ConstIterator operator++()         { return --m_Iterator; }
    ConstIterator operator++(int)      { return m_Iterator--; }
    ConstIterator operator--()         { return ++m_Iterator; }
    ConstIterator operator--(int)      { return m_Iterator++; }
    ConstIterator operator->() const { return ( m_Iterator - 1 ); }
    const ValueType & operator*() const { return *( m_Iterator - 1 ); }
    bool operator!=(const ConstReverseIterator & rit) const { return m_Iterator != rit.m_Iterator; }
    bool operator==(const ConstReverseIterator & rit) const { return m_Iterator == rit.m_Iterator; }

  private:
    ConstIterator m_Iterator;
  };
  
  /**
   * Default constructor uses compiler's default initialization of memory.
   * For efficiency, no initialization to zero is done.
   */
  NDArrayBase() : m_IsShared(false), m_Data(NULL)
  {
  for ( int i = 0; i < Dimension; ++i )
    m_Shape[i]=0, m_OffSetTable[i]=0; 
  }

  explicit NDArrayBase(const ShapeType& shape) : m_IsShared(false), m_Data(NULL)
  {
  __utl_ndarray_alloc_blah(shape);
  }

  NDArrayBase(const NDArrayBase<T,Dim>& vec) : m_IsShared(false), m_Data(NULL)
  {
  for ( int i = 0; i < Dimension; ++i )
    m_Shape[i]=0, m_OffSetTable[i]=0; 
  operator=(vec);
  }
  
  template<typename EType>
  NDArrayBase(const Expr<EType>& expr) : m_IsShared(false), m_Data(NULL)
  {
  for ( int i = 0; i < Dimension; ++i )
    m_Shape[i]=0, m_OffSetTable[i]=0; 
  operator=(expr);
  }

/**
 * Constructor assumes input points to array of correct size.
 * Values are copied individually instead of with a binary copy.  This
 * allows the T's assignment operator to be executed.
 */
  NDArrayBase(const T* vec, const ShapeType& shape) : m_IsShared(false), m_Data(NULL)
  {
  __utl_ndarray_alloc_blah(shape);
  utl::cblas_copy<T>(this->Size(), vec, 1, m_Data, 1);
  }

  /**
   * Constructor to initialize entire array to one value.
   */
  NDArrayBase(const ShapeType& shape, const T r) : m_IsShared(false), m_Data(NULL)
  {
  __utl_ndarray_alloc_blah(shape);
  Fill(r);
  }

  template< typename TValue >
  NDArrayBase(const NDArrayBase<TValue, Dim> & r) : m_IsShared(false),  m_Data(NULL)
  {
  for ( int i = 0; i < Dimension; ++i )
    m_Shape[i]=0, m_OffSetTable[i]=0; 
  operator=(r);
  }
  

  /** This destructor is not virtual for performance reasons. However, this
   * means that subclasses cannot allocate memory.
   *
   * The destructor is PURPOSELY NOT DEFINED, in order to prevent inefficient
   * byte alignment of arrays of this object.
   *
   *
   * For a full discussion, see
   * http://www.itk.org/mailman/private/insight-developers/2008-June/010480.html
   *
   */
  ~NDArrayBase()
    { Clear(); }

  void Clear()
    {
    ClearData();
    ClearShape();
    m_Data = NULL;
    m_IsShared = false;
    }
  
  UTL_ALWAYS_INLINE double Eval( int i ) const
    {
    return m_Data[i];
    }
  
  /**
   * Get the size of the NDArrayBase.
   */
  UTL_ALWAYS_INLINE SizeType   Size() const
    {
    return m_OffSetTable[0];
    }
  UTL_ALWAYS_INLINE SizeType   GetSize() const
    {
    return m_OffSetTable[0];
    }

  UTL_ALWAYS_INLINE const ShapeType GetShape() const
  {
  return m_Shape;
  }
  
  UTL_ALWAYS_INLINE static SizeType  GetDimension() 
    { return Dimension; }

  UTL_ALWAYS_INLINE SizeType GetOffset(const ShapeType& shapeIndex) const
    {
    SizeType index = 0;
    for ( int i = 0; i < Dimension-1; ++i ) 
      index += shapeIndex[i]*m_OffSetTable[i+1];
    index += shapeIndex[Dimension-1];
    return index;
    }
    

#define __Array_Saver_Expr(saver, saverReal)                                                     \
template<typename EType>                                                                         \
inline NDArrayBase<T,Dim>& operator saver (const Expr<EType>& src){                              \
  SizeType srcDim = Expr<EType>::GetDimension();                                                 \
  utlSAException(srcDim>0 && srcDim!=Dimension)                                                  \
    (srcDim)(Dimension).msg("the expression has a difference size");                             \
  const EType &r = src.ConstRef();                                                               \
  const ShapeType rShape = r.GetShape();                                                         \
  if (srcDim>0)                                                                                  \
    this->ReSize(rShape);                                                                        \
  for( int i=0; i < this->Size(); ++i )                                                          \
    this->m_Data[i] saverReal r.Eval(i);                                                         \
  return *this;                                                                                  \
}                                                                                                


  __Array_Saver_Expr(=,  =)
  __Array_Saver_Expr(+=, +=)
  __Array_Saver_Expr(-=, -=)
  /** use %= for elementwse *=  */
  __Array_Saver_Expr(%=, *=)
  __Array_Saver_Expr(/=, /=)

#define __Array_Saver_Scalar(saver)                        \
inline NDArrayBase<T,Dim>& operator saver (const T val){   \
  for (Iterator p = Begin(); p!= End(); ++p)               \
    *p saver val;                                          \
  return *this;                                            \
}                                                      

  __Array_Saver_Scalar(+=)
  __Array_Saver_Scalar(-=)

  // these three operators use blas and mkl
  // __Array_Saver_Scalar(=)
  // __Array_Saver_Scalar(*=)
  // __Array_Saver_Scalar(/=)

  NDArrayBase<T,Dim>& operator=(const NDArrayBase<T,Dim> & r)
    {
    // std::cout << "copy matrix" << std::endl << std::flush;
    if ( this != &r ) 
      {
      this->ReSize(r.GetShape());
      utl::cblas_copy(Size(), r.Begin(), 1, m_Data, 1);
      }
    return *this;
    }
  
/** Operator= defined for a variety of types. */
  template< typename TValueType >
  NDArrayBase<T,Dim>& operator=(const NDArrayBase< TValueType, Dim > & r)
    {
    this->ReSize(r.GetShape());
    typename NDArrayBase< TValueType, Dim>::ConstIterator input = r.Begin();
    Iterator i = this->Begin();
    while ( i != this->End() )
      *i++ = static_cast< T >( *input++ );
    return *this;
    }

  /**
   * Assignment operator assumes input points to array of correct size.
   * Values are copied individually instead of with a binary copy.  This
   * allows the T's assignment operator to be executed.
   *
   * If Dimension==1, then ReSize is used. 
   * If Dimension>1, the sizes of the array and the input std::vector should be the same. 
   */
  NDArrayBase<T,Dim> & operator=(const std::vector<T>& r)
    {
    if (r.size()==Size())
      utl::cblas_copy(Size(), r.data(), 1, m_Data, 1);
    else if (Dimension==1)
      {
      SizeType shape[Dim];
      shape[0] = r.size();
      ReSize(shape);
      utl::cblas_copy(Size(), r.data(), 1, m_Data, 1);
      }
    else
      utlException(Dimension!=1, "should have only 1 dimension, or have the same size.");
    return *this;
    }

  NDArrayBase<T,Dim> & operator=(const T r)
    { Fill(r); return *this; }

  /** Operators == and != are used to compare whether two arrays are equal.
   * Note that arrays are equal when the number of components (size) is the
   * same, and each component value is equal. */
  bool operator==(const NDArrayBase<T,Dim> & r) const
    {
    if (!IsSameShape(r))
      return false;
    if (r.m_Data==m_Data)
      return true;
    ConstIterator i = this->Begin();
    ConstIterator j = r.Begin();

    while ( i != this->End() )
      {
      if ( *i != *j )
        return false;
      ++j;
      ++i;
      }
    return true;
    }

  bool operator!=(const NDArrayBase<T,Dim> & r) const
  { return !operator==(r); }

  bool IsEqual(const NDArrayBase<T,Dim>& r, const double eps) const
    {
    if (!IsSameShape(r))
      return false;
    if (r.m_Data==m_Data)
      return true;
    ConstIterator i = this->Begin();
    ConstIterator j = r.Begin();

    while ( i != this->End() )
      {
      if ( std::fabs(*i - *j) > eps )
        return false;
      ++j;
      ++i;
      }
    return true;
    }
  
  bool IsSameValues(const NDArrayBase<T,Dim>& r, const double eps) const
    {
    if (!IsSameSize(r))
      return false;
    if (r.m_Data==m_Data)
      return true;
    ConstIterator i = this->Begin();
    ConstIterator j = r.Begin();

    while ( i != this->End() )
      {
      if ( std::fabs(*i - *j) > eps )
        return false;
      ++j;
      ++i;
      }
    return true;
    }
  

  /** Allow the NDArrayBase to be indexed normally.  No bounds checking is done.
   * The separate versions are a work-around for an integer conversion bug in
   * Visual C++. */
  UTL_ALWAYS_INLINE Reference operator[](short index)                { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](short index) const { return m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator[](unsigned short index)       { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](unsigned short index) const { return m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator[](int index)                  { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](int index) const { return m_Data[index]; }
// false positive warnings with GCC 4.9
#if defined( __GNUC__ )
#if ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ == 9 )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
#endif
  UTL_ALWAYS_INLINE Reference operator[](unsigned int index)         { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](unsigned int index) const { return m_Data[index]; }
#if defined( __GNUC__ )
#if ( __GNUC__ == 4 ) && ( __GNUC_MINOR__ == 9 )
#pragma GCC diagnostic pop
#endif
#endif
  UTL_ALWAYS_INLINE Reference operator[](long index)                 { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](long index) const { return m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator[](unsigned long index)        { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](unsigned long index) const { return m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator[](long long index)                 { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](long long index) const { return m_Data[index]; }
  UTL_ALWAYS_INLINE Reference operator[](unsigned long long index)        { return m_Data[index]; }
  UTL_ALWAYS_INLINE ConstReference operator[](unsigned long long index) const { return m_Data[index]; }

  UTL_ALWAYS_INLINE Reference operator()(const ShapeType& shape)         { return m_Data[GetOffset(shape)]; }
  UTL_ALWAYS_INLINE ConstReference operator()(const ShapeType& shape) const { return m_Data[GetOffset(shape)]; }

  UTL_ALWAYS_INLINE Reference Back()            { return m_Data[Size()-1]; }
  UTL_ALWAYS_INLINE ConstReference Back() const { return m_Data[Size()-1]; }

  /**
   * Get an Iterator for the beginning of the NDArrayBase.
   */
  Iterator      Begin()
    {
    return Iterator(m_Data);
    }

  /**
   * Get a ConstIterator for the beginning of the NDArrayBase.
   */
  ConstIterator Begin() const
    {
    return ConstIterator(m_Data);
    }

  /**
   * Get an Iterator for the end of the NDArrayBase.
   */
  Iterator      End()
    {
    return Iterator(m_Data + Size());
    }

  /**
   * Get a ConstIterator for the end of the NDArrayBase.
   */
  ConstIterator End() const
    {
    return ConstIterator(m_Data + Size());
    }

  /**
   * Get a begin ReverseIterator.
   */
  ReverseIterator      rBegin()
    {
    return ReverseIterator(m_Data + Size());
    }

  /**
   * Get a begin ConstReverseIterator.
   */
  ConstReverseIterator rBegin() const
    {
    return ConstReverseIterator(m_Data + Size());
    }

  /**
   * Get an end ReverseIterator.
   */
  ReverseIterator      rEnd()
    {
    return ReverseIterator(m_Data);
    }

  /**
   * Get an end ConstReverseIterator.
   */
  ConstReverseIterator rEnd() const
    {
    return ConstReverseIterator(m_Data);
    }

  /** Set/Get element methods are more convenient in wrapping languages */
  void SetElement(unsigned short index, ConstReference value)
    { m_Data[index] = value; }
  ConstReference GetElement(unsigned short index) const 
    { return m_Data[index]; }

  /** Return a pointer to the data. */
  UTL_ALWAYS_INLINE T* GetData()
  {
    return m_Data; 
  }

  UTL_ALWAYS_INLINE const T* GetData() const
  {
    return m_Data; 
  }

  /** set m_Data from data, the ownership is in data. */
  void SetData( T* const data, const ShapeType& shape )
    {
    utlException(data==m_Data, "cannot set itself. Please use ReShape() for change the shape");
    Clear();
    m_Data = data;
    m_IsShared = true;
    for ( int i = 0; i < Dimension; ++i ) 
      m_Shape[i] = shape[i];
    ComputeOffSetTable();
    }

  /** copy m_Data from data  */
  void CopyData(T* const data, const ShapeType& shape)
    {
    ReSize(shape);
    utl::cblas_copy(Size(), data, 1, m_Data, 1);
    }

  /** If the current size is different from the size in shape, then data allocation happens. 
   * Otherwise, just change the shape using ReShape. */
  bool ReSize(const ShapeType& shape)
  {
  // if no change in shape or size, do not reallocate.
  if (IsSameShape(shape)) 
    return false;
  if (IsSameSize(shape)) 
    {
    ReShape(shape);
    return false;
    }
  if (this->m_Data) 
    {
    if (!m_IsShared)
      Clear();
    __utl_ndarray_alloc_blah(shape);
    }
  else 
    {
    // this happens if the array is default constructed.
    __utl_ndarray_alloc_blah(shape);
    }
  return true;
  }

  /** change the shape, no data allocation (same size)  */
  inline NDArrayBase<T,Dim>& ReShape(const ShapeType& shape)
  {
  SizeType newSize = std::accumulate(shape, shape+Dim, 0);
  utlException(newSize!=Size(), "need to keep the same size");
  for ( int i = 0; i < Dimension; ++i ) 
    m_Shape[i] = shape[i];
  ComputeOffSetTable();
  return *this;
  }
  
  template<typename EType>                                      
  UTL_ALWAYS_INLINE bool IsSameShape(const EType& src) const
  {
  SizeType srcDim = Expr<EType>::GetDimension();                                                 
  if (srcDim==0)
    return true;
  if (srcDim!=Dimension)
    return false;
  const EType &r = src.ConstRef();            
  const ShapeType rShape = r.GetShape();     
  for ( int i = 0; i < Dimension; ++i )                                                        
    {
    if (rShape[i]!=m_Shape[i])
      return false;
    }
  return true;
  }
  

  UTL_ALWAYS_INLINE bool IsSameShape(const ShapeType& shape) const
  {
  // utl::PrintContainer(m_Shape, m_Shape+Dimension, "m_Shape", " ");
  // utl::PrintContainer(shape, shape+Dimension, "shape", " ");
  for ( int i = 0; i < Dimension; ++i ) 
    {
    if (shape[i]!=m_Shape[i])
      return false;
    }
  return true;
  }
  
  UTL_ALWAYS_INLINE bool IsSameSize(const ShapeType& shape) const
  {
  return  Size() == std::accumulate(shape, shape+Dim, 0);
  }
  
  /** y = a+x, make sure vec and outVec have valid length. 
   * when outVec is null, the result is stored in *this. */
  NDArrayBase<T,Dim>& ElementAdd(T* const vec, T* outVec=NULL)
    { utl::vAdd(Size(), m_Data, vec, outVec?outVec:m_Data); return *this; }

  NDArrayBase<T,Dim>& ElementSubstract(T* const vec, T* outVec=NULL)
    { utl::vSub(Size(), m_Data, vec, outVec?outVec:m_Data); return *this; }

  NDArrayBase<T,Dim>& ElementMultiply(T* const vec, T* outVec=NULL)
    { utl::vMul(Size(), m_Data, vec, outVec?outVec:m_Data); return *this; }

  NDArrayBase<T,Dim>& ElementDivide(T* const vec, T* outVec=NULL)
    { utl::vDiv(Size(), m_Data, vec, outVec?outVec:m_Data); return *this; }

  NDArrayBase<T,Dim>& ElementExp(T* outVec=NULL)
    { utl::vExp(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementExp(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementExp(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementExp() const
    { return NDArrayBase<T,Dim>(*this).ElementExp(); }

  NDArrayBase<T,Dim>& ElementCos(T* outVec=NULL)
    { utl::vCos(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementCos(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementCos(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementCos() const
    { return NDArrayBase<T,Dim>(*this).ElementCos(); }

  NDArrayBase<T,Dim>& ElementSin(T* outVec=NULL)
    { utl::vSin(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementSin(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementSin(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementSin() const
    { return NDArrayBase<T,Dim>(*this).ElementSin(); }

  NDArrayBase<T,Dim>& ElementAbsolute(T* outVec=NULL)
    { utl::vAbs(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementAbsolute(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementAbsolute(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementAbsolute() const
    { return NDArrayBase<T,Dim>(*this).ElementAbsolute(); }

  NDArrayBase<T,Dim>& ElementInverse(T* outVec=NULL)
    { utl::vInv(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementInverse(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementInverse(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementInverse() const
    { return NDArrayBase<T,Dim>(*this).ElementInverse(); }

  NDArrayBase<T,Dim>& ElementSqrt(T* outVec=NULL)
    { utl::vSqrt(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementSqrt(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementSqrt(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementSqrt() const
    { return NDArrayBase<T,Dim>(*this).ElementSqrt(); }

  NDArrayBase<T,Dim>& ElementSquare(T* outVec=NULL)
    { utl::vSqr(Size(), m_Data, outVec?outVec:m_Data); return *this; }
  NDArrayBase<T,Dim>& ElementSquare(NDArrayBase<T,Dim>& vec)
    { vec.Resize(m_Shape); ElementSquare(vec.m_Data); return *this; }
  NDArrayBase<T,Dim> GetElementSquare() const
    { return NDArrayBase<T,Dim>(*this).ElementSquare(); }
  
  /** m_Data[i] = alpha*vec[i] + beta*m_Data[i] */
  NDArrayBase<T,Dim>& ElementAxpby(T* const vec, const T alpha, const T beta)
    { utl::cblas_axpby(Size(), alpha, vec, 1, beta, m_Data, 1); return *this; }

  NDArrayBase<T,Dim>& Scale(const T a) 
    { 
    if (std::fabs(a-1.0)>1e-10) 
      cblas_scal<T>(Size(),a,m_Data,1); 
    return *this; 
    }

  NDArrayBase<T,Dim>& operator+=(const NDArrayBase<T,Dim>& vec) 
    {
    utlException(!IsSameShape(vec), "wrong shape");
    utlSAException(vec.Size()!=Size())(Size())(vec.Size()).msg("wrong size");
    ElementAdd(vec.m_Data);
    return *this;
    }
  NDArrayBase<T,Dim>& operator-=(const NDArrayBase<T,Dim>& vec) 
    {
    utlException(!IsSameShape(vec), "wrong shape");
    utlSAException(vec.Size()!=Size())(Size())(vec.Size()).msg("wrong size");
    ElementSubstract(vec.m_Data);
    return *this;
    }
  // [>* use % for multiplication  <]
  NDArrayBase<T,Dim>& operator%=(const NDArrayBase<T,Dim>& vec) 
    {
    utlSAException(vec.Size()!=Size())(Size())(vec.Size()).msg("wrong size");
    ElementMultiply(vec.m_Data);
    return *this;
    }
  NDArrayBase<T,Dim>& operator/=(const NDArrayBase<T,Dim>& vec) 
    {
    utlSAException(vec.Size()!=Size())(Size())(vec.Size()).msg("wrong size");
    ElementDivide(vec.m_Data);
    return *this;
    }
  
  NDArrayBase<T,Dim>& operator%=(const T val)
    {
    Scale(val);
    return *this;
    }
  NDArrayBase<T,Dim>& operator/=(const T val)
    {
    Scale(1.0/val);
    return *this;
    }

  /**
   * Fill all elements of the array with the given value.
   */
  void Fill(const T& value)
    {
    std::fill(Begin(), End(), value);
    }

  /** Set elements from ptr (starting with m_Data + start )*/ 
  NDArrayBase<T,Dim>& CopyIn(T const * ptr, const int size, const int start=0)
    {
    utlException(size+start>Size(), "wrong size");
    utl::cblas_copy(size, ptr, 1, m_Data+start, 1);
    return *this;
    }

  /** Set elements to ptr (starting with m_Data + start )*/ 
  void CopyOut(T * ptr, const int size, const int start=0) const // from vector to array[].
    {
    utlException(size+start>Size(), "wrong size");
    utl::cblas_copy(size, m_Data+start, 1, ptr, 1);
    }

  
  /** Applies function to elements  */
  template<typename FuncT>
  void Apply(const FuncT& func, NDArrayBase<T,Dim>& vec) const
    {
    vec.ReSize(m_Shape);
    for (SizeType i = 0; i < Size(); ++i)
      vec[i] = func(m_Data[i]);
    }
  void Apply(T (*f)(T), NDArrayBase<T,Dim>& vec) const
    {
    vec.ReSize(m_Shape);
    for (SizeType i = 0; i < Size(); ++i)
      vec[i] = f(m_Data[i]);
    }
  void Apply(T (*f)(T const&), NDArrayBase<T,Dim>& vec) const
    {
    vec.ReSize(m_Shape);
    for (SizeType i = 0; i < Size(); ++i)
      vec[i] = f(m_Data[i]);
    }

  unsigned ArgMax() const
    {
    return utl::argmax(Begin(), End());
    }
  
  T MaxValue() const
    {
    return *std::max_element(Begin(), End());
    }
  
  unsigned ArgMin() const
    {
    return utl::argmin(Begin(), End());
    }
  
  T MinValue() const
    {
    return *std::min_element(Begin(), End());
    }

  unsigned ArgAbsoluteMax() const
    {
    return  utl::cblas_iamax<T>(Size(), m_Data, 1);
    }

  T AbsoluteMaxValue() const
    {
    unsigned index =  ArgAbsoluteMax();
    return m_Data[index];
    }
  
  unsigned ArgAbsoluteMin() const
    {
    return  utl::cblas_iamin(Size(), m_Data, 1);
    }
  T AbsoluteMinValue() const
    {
    unsigned index =  ArgAbsoluteMin();
    return m_Data[index];
    }
  
  /** Return largest absolute element value  */
  T GetInfNorm( ) const
    {
    double absMax = AbsoluteMaxValue();
    return absMax>=0 ? absMax : -absMax;
    }

  /**  sqrt of sum of squares of its elements */
  T GetTwoNorm( ) const
    {
    return utl::cblas_nrm2(Size(), m_Data, 1);
    }
  T GetSquaredTwoNorm( ) const
    {
    double norm = GetTwoNorm();
    return norm*norm;
    }
  T GetRootMeanSquares( ) const
    {
    return GetTwoNorm()/std::sqrt(T(Size()));
    }

  /** Return sum of absolute values of elements  */
  T GetOneNorm( ) const
    {
    return utl::cblas_asum(Size(), m_Data, 1);
    }
  
  /** number of non-zero values  */
  int GetZeroNorm( const double eps=1e-10 ) const
    {
    return NNZ(eps);
    }

  /** get sum  */
  T GetSum() const
    {
    T tot(0);
    for (ConstIterator p = Begin(); p!= End(); ++p)
      tot += *p;
    return tot;
    }
  
  /** get mean  */
  T GetMean() const
    {
    return GetSum()/(T)Size();
    }

  T GetMedian() const
    {
    NDArrayBase<T,Dim> vec(*this);
    const unsigned int N = vec.Size();
    std::nth_element(vec.Begin(),vec.Begin()+N/2,vec.End());
    return vec[N/2];
    }

  UTL_ALWAYS_INLINE void Swap(NDArrayBase<T,Dim>& vec)
    {
    std::swap(m_IsShared, vec.m_IsShared);
    std::swap(m_Shape, vec.m_Shape);
    std::swap(m_OffSetTable, vec.m_OffSetTable);
    std::swap(m_Data, vec.m_Data);
    }

  NDArrayBase<T,Dim>& Flip()
    {
    for (unsigned i=0; i<Size()/2; i++) 
      std::swap(m_Data[i], m_Data[Size()-1-i]);
    return *this;
    }

  bool IsZero() const
    {
    T const zero(0);
    for (unsigned i = 0; i < Size();++i)
      if ( !( m_Data[i] == zero) )
        return false;
    return true;
    }

  bool IsEmpty() const
    { return Size()==0; }

  /** number of non-zero values  */
  int NNZ(const double eps=1e-10) const 
    {
    int sum=0;
    for (int i = 0; i<Size(); ++i) 
      {
      if (std::fabs(m_Data[i])<eps) 
        ++sum;
      }
    return sum;
    }

  void SoftThreshold(const double threshold)
    {
    for (int i = 0; i<Size(); ++i) 
      {
      if (m_Data[i] > threshold) 
        m_Data[i] -= threshold;
      else if (m_Data[i] < -threshold)
        m_Data[i] += threshold;
      else 
        m_Data[i] = 0.0;
      }
    }

  void HardThreshold(const double threshold)
    {
    for (int i = 0; i<Size(); ++i) 
      {
      if (!(m_Data[i] > threshold || m_Data[i] < -threshold)) 
        m_Data[i] = 0;
      }
    }

  ValueType InnerProduct(const NDArrayBase<T,Dim>& vec) const
    { return utl::InnerProduct(*this, vec); }

  void 
  Print(std::ostream & os, const char* separate=" ") const
  {
  PrintVar1(true, Dimension, os);
  utl::PrintContainer(m_Shape, m_Shape+Dimension, "m_Shape", separate, os);
  utl::PrintContainer(m_OffSetTable, m_OffSetTable+Dimension,"m_OffSetTable", separate, os);
  utl::PrintContainer(m_Data, m_Data+Size(), "m_Data", separate, os);
  }


protected:

  /** Internal C array representation. */
  T* m_Data;

  SizeType m_OffSetTable[Dimension];
  SizeType m_Shape[Dimension];

  /** If m_IsShared is true, the memory is owned by other object  */
  bool m_IsShared;

  UTL_ALWAYS_INLINE void ClearData()
    {
    if (!m_IsShared && Size()>0)
      delete [] m_Data;
    }
  UTL_ALWAYS_INLINE void ClearShape()
    {
    for ( int i = 0; i < Dimension; ++i ) 
      m_Shape[i]=0, m_OffSetTable[i]=0;
    }

  UTL_ALWAYS_INLINE void ComputeOffSetTable()
    {
    for ( int i = Dimension-1; i >= 0; i-- ) 
      { 
      if (i==Dimension-1)
        m_OffSetTable[i]=m_Shape[i]; 
      else
        m_OffSetTable[i]=m_Shape[i]*m_OffSetTable[i+1]; 
      }
    }

};

/** @}  */

}

#include "utlVector.h"
#include "utlMatrix.h"


#endif 
