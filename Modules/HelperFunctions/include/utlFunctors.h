/**
 *       @file  utlFunctors.h
 *      @brief  functors in utl
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlFunctors_h
#define __utlFunctors_h

#include "utlCoreMacro.h"
#include "utlCore.h"
#include "utlExpression.h"
#include <cmath>


/** macro to define unary functors in utl::Functor */
#define utlUnaryFunctorBaseMacro(name, realFunc)                                                                         \
namespace Functor{                                                                                                       \
template<typename T>                                                                                                     \
struct name                                                                                                              \
{                                                                                                                        \
    UTL_ALWAYS_INLINE  T operator()(const T a) const                                                                     \
      { return realFunc(a); }                                                                                            \
    UTL_ALWAYS_INLINE bool operator!=(const name &) const                                                                \
      { return false; }                                                                                                  \
    UTL_ALWAYS_INLINE bool operator==(const name & other) const                                                          \
      { return !( *this != other ); }                                                                                    \
};                                                                                                                       \
}                                                                                                                        


/** macro to define unary functors in utl::Functor, then define functions in utl which can be used elementwisely for expressions. */
#define utlUnaryFunctorMacro(name, realFunc)                                                                             \
  utlUnaryFunctorBaseMacro(name, realFunc)                                                                               \
                                                                                                                         \
template<class ExprT>                                                                                                    \
auto name(const ExprT& expr) -> decltype(utl::F<utl::Functor::name<typename ExprT::ValueType> >(expr))                   \
  {                                                                                                                      \
  return utl::F<utl::Functor::name<typename ExprT::ValueType> >(expr);                                                   \
  }                                                                                                                      \
                                                                                                                         \
namespace Functor{                                                                                                       \
template<typename T, unsigned Dim>                                                                                       \
struct name<utl::NDArray<T,Dim> > : public VectorFunctorBase<utl::NDArray<T,Dim>, utl::NDArray<T,Dim> >                  \
{                                                                                                                        \
public:                                                                                                                  \
  bool operator!=(const name & other) const                                                                              \
  { return false; }                                                                                                      \
  bool operator==(const name & other) const                                                                              \
  { return !( *this != other ); }                                                                                        \
  void operator=(const name & other)  {}                                                                                 \
  inline utl::NDArray<T,Dim> operator()( const utl::NDArray<T,Dim> & A ) const                                           \
  {                                                                                                                      \
    return utl::name(A);                                                                                                 \
  }                                                                                                                      \
                                                                                                                         \
  int GetOutputDimension(const int inSize) const                                                                         \
    {  return inSize; }                                                                                                  \
};                                                                                                                       \
}



/** macro to define binary functors in utl::Functor */
#define utlBinaryFunctorBaseMacro(name, realFunc)                                                                                                     \
namespace Functor{                                                                                                                                    \
template<typename T>                                                                                                                                  \
struct name                                                                                                                                           \
{                                                                                                                                                     \
    UTL_ALWAYS_INLINE  T operator()(const T a, const T b) const                                                                                       \
      { return realFunc(a, b); }                                                                                                                      \
    UTL_ALWAYS_INLINE bool operator!=(const name &) const                                                                                             \
      { return false; }                                                                                                                               \
    UTL_ALWAYS_INLINE bool operator==(const name & other) const                                                                                       \
      { return !( *this != other ); }                                                                                                                 \
};                                                                                                                                                    \
}                                                                                                                                                     


/** macro to define binary functors in utl::Functor, then define functions in utl which can be used elementwisely for expressions. */
#define utlBinaryFunctorMacro(name, realFunc)                                                                                                         \
  utlBinaryFunctorBaseMacro(name, realFunc)                                                                                                           \
                                                                                                                                                      \
template<class TLeft, class TRight>                                                                                                                   \
auto name(const TLeft& lhs, const TRight& rhs) -> decltype(utl::F<utl::Functor::name<Expr2ValueType<TLeft,TRight>> >(lhs, rhs))                       \
  {                                                                                                                                                   \
  return utl::F<utl::Functor::name<Expr2ValueType<TLeft,TRight>> >(lhs, rhs);                                                                         \
  }                                                                                                                                                   \
template<class TLeft>                                                                                                                                 \
auto name(const TLeft& lhs, const ScalarExpr& rhs) -> decltype(utl::F<utl::Functor::name<Expr2ValueType<TLeft,ScalarExpr>> >(lhs, rhs))               \
  {                                                                                                                                                   \
  return utl::F<utl::Functor::name<Expr2ValueType<TLeft,ScalarExpr>> >(lhs, rhs);                                                                     \
  }                                                                                                                                                   \
template<class TRight>                                                                                                                                \
auto name(const ScalarExpr& lhs, const TRight& rhs) -> decltype(utl::F<utl::Functor::name<Expr2ValueType<ScalarExpr,TRight>> >(lhs, rhs))             \
  {                                                                                                                                                   \
  return utl::F<utl::Functor::name<Expr2ValueType<ScalarExpr,TRight>> >(lhs, rhs);                                                                    \
  }                                                                                                                                                   \


/** define functors from utl::Vector to double */
#define utlVectorToScalarFunctorMacro(name, funcName)           \
namespace Functor{                                              \
template< class TVector, class TOutput=double >                 \
class name : public VectorFunctorBase<TVector,TOutput>          \
{                                                               \
public:                                                         \
  bool operator!=(const name & other) const                     \
  { return false; }                                             \
  bool operator==(const name & other) const                     \
  { return !( *this != other ); }                               \
  void operator=(const name & other)  {}                        \
  inline TOutput operator()( const TVector & A ) const          \
  {                                                             \
    return A.funcName();                                        \
  }                                                             \
                                                                \
  int GetOutputDimension(const int) const                       \
    {  return 1; }                                              \
};                                                              \
}                                                             

namespace utl
{

template < class T, unsigned int Dim >
class NDArray;

/** 
*
* @addtogroup utlMath
*  @{
* */


namespace Functor 
{
template< class TVector, class TOutput=TVector >  
class VectorFunctorBase
{ 
public:
  typedef VectorFunctorBase Self;
  bool operator==(const Self & other) const  
  { return false; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 

  Self& operator=(const Self & other)
    {
    if (this!=&other)
      {
      m_LogLevel = other.m_LogLevel;
      }
    return *this;
    }

  int GetOutputDimension(const int inputSize) const
    {return -1;}
  int GetOutputDimension(const std::vector<int>& sizeVec) const
    {return -1;}

  void Initialize() {}
  void VerifyInputParameters(const int inputSize=-1) const {}

  utlSetGetMacro(LogLevel, int);
  
  void Print(std::ostream & os=std::cout) const
    {
    utlLogOSVar(os, m_LogLevel);
    }

protected:
  int m_LogLevel=1;
};


template< class TFunctor, class TVector=utl::NDArray<double,1>, class TOutput=TVector >  
class ScalarFunctorWrapper : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef ScalarFunctorWrapper Self;
  bool operator==(const Self & other) const  
  { return false; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  void operator=(const Self & other)   
    {m_Functor = other.m_Functor;}    

  ScalarFunctorWrapper()
    {
    m_Functor = TFunctor();
    }

  inline TOutput operator()( const TVector & A) const 
  {
  int inputSize = A.Size();
  TOutput out(inputSize);
  for ( int i = 0; i < inputSize; ++i ) 
    {
    out[i]= m_Functor(A[i]);
    }
  return out;
  }
  
  int GetOutputDimension(const int inputSize) const
    {
    return 1;
    } 
protected:
  TFunctor m_Functor;
};

template< class TFunctor, class TVector=utl::NDArray<double,1>, class TOutput=TVector >  
class VectorFunctorWrapper : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef VectorFunctorWrapper Self;
  bool operator==(const Self & other) const  
  { return false; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  void operator=(const Self & other)    
    {m_Functor = other.m_Functor;}    

  VectorFunctorWrapper()
    {
    m_Functor = TFunctor();
    }

  inline TOutput operator()( const TVector & A) const 
  {
  int inputSize = A.Size();
  TOutput out(inputSize);
  out= m_Functor(A);
  return out;
  }
  
  int GetOutputDimension(const int inputSize) const
    {
    return inputSize;
    } 
protected:
  TFunctor m_Functor;
};

template< class TFunction=std::function<double(double)> , class TVector=utl::NDArray<double,1>, class TOutput=TVector >  
class VectorUnaryFunctionWrapper : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef VectorUnaryFunctionWrapper Self;
  bool operator==(const Self & other) const  
  { return false; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  void operator=(const Self & other)  
    {m_Function = other.m_Function;}    

  VectorUnaryFunctionWrapper(TFunction func=nullptr)
    {
    m_Function = func; 
    }

  inline TOutput operator()( const TVector & A) const 
  {
  int inputSize = A.Size();
  TOutput out(inputSize);
  for ( int i = 0; i < inputSize; ++i ) 
    {
    out[i]= m_Function(A[i]);
    }
  return out;
  }
  
  int GetOutputDimension(const int inputSize) const
    {
    return inputSize;
    } 
protected:
  TFunction m_Function;
};

template< class TFunction=std::function<double(std::vector<double>)> , class TVector=utl::NDArray<double,1>, class TOutput=TVector >  
class VectorMultiVariableFunctionWrapper : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef VectorMultiVariableFunctionWrapper Self;
  bool operator==(const Self & other) const  
  { return false; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  void operator=(const Self & other)  
    {m_Function = other.m_Function;}    

  VectorMultiVariableFunctionWrapper(TFunction func=nullptr)
    {
    m_Function = func; 
    }

  inline TOutput operator()( const std::vector<TVector>& vec) const 
  {
  int N = vec.size();
  utlException(N==0, "cannot be empty vector");
  int inputSize = vec[0].Size();
  TOutput out(inputSize);
  std::vector<double> values;
  for ( int i = 0; i < inputSize; ++i ) 
    {
    values.clear();
    for ( int j = 0; j < N; ++j ) 
      {
      values.push_back(vec[j][i]);
      }
    out[i]= m_Function(values);
    }
  return out;
  }
  
  int GetOutputDimension(const int inputSize) const
    {
    return inputSize;
    } 
  int GetOutputDimension(const std::vector<int>& sizeVec) const
    {
    return sizeVec[0];
    } 
protected:
  TFunction m_Function;
};

}


/** Unary functions  */
utlUnaryFunctorMacro(Abs, std::abs);
utlUnaryFunctorMacro(Exp, std::exp);
utlUnaryFunctorMacro(Exp2, std::exp2);
utlUnaryFunctorMacro(Log, std::log);
utlUnaryFunctorMacro(Log10, std::log10);
utlUnaryFunctorMacro(Log2, std::log2);
utlUnaryFunctorMacro(Sqrt, std::sqrt);
utlUnaryFunctorMacro(Floor, std::floor);
utlUnaryFunctorMacro(Round, std::round);
utlUnaryFunctorMacro(LRound, std::lround);
utlUnaryFunctorMacro(Neg, -);
utlUnaryFunctorMacro(Sign, utl::sign);
utlUnaryFunctorMacro(Square, utl::square);
utlUnaryFunctorMacro(Cube, utl::cube);

utlUnaryFunctorMacro(Sin, std::sin);
utlUnaryFunctorMacro(Cos, std::cos);
utlUnaryFunctorMacro(Tan, std::tan);
utlUnaryFunctorMacro(Asin, std::asin);
utlUnaryFunctorMacro(Acos, std::acos);
utlUnaryFunctorMacro(Atan, std::atan);

// utlUnaryFunctorMacro(Real, std::real);  // return type: double
// utlUnaryFunctorMacro(Imag, std::imag);  // return type: double
utlUnaryFunctorMacro(Conj, std::conj);

/** Binary functions  */
utlBinaryFunctorMacro(Max, std::max);
utlBinaryFunctorMacro(Min, std::min);
utlBinaryFunctorMacro(Pow, std::pow);
utlBinaryFunctorMacro(Atan2, std::atan2);


/** define functors from utl::Vector to double */
utlVectorToScalarFunctorMacro(Mean, GetMean);
utlVectorToScalarFunctorMacro(Sum, GetSum);
utlVectorToScalarFunctorMacro(Median, GetMedian);
utlVectorToScalarFunctorMacro(TwoNorm, GetTwoNorm);
utlVectorToScalarFunctorMacro(OneNorm, GetOneNorm);
utlVectorToScalarFunctorMacro(InfNorm, GetInfNorm);
utlVectorToScalarFunctorMacro(ZeroNorm, GetZeroNorm);
utlVectorToScalarFunctorMacro(MaxValue, MaxValue);
utlVectorToScalarFunctorMacro(MinValue, MinValue);

namespace Functor
{


/** 
*   \class   Shred
*   \brief  obtain a part of a vector. 
*
* For a vector with size m_Offset + N*(m_ChunkSize + m_Space) + residual, 
* The output vector has size N*m_ChunkSize + std::min(residual, m_ChunkSize). 
* */
template< class TVector, class TOutput=TVector >  
class Shred : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef Shred Self;
  typedef VectorFunctorBase<TVector, TOutput> Superclass;

  bool operator==(const Self & other) const  
  { return m_Offset==other.m_Offset && m_ChunkSize==other.m_ChunkSize && m_Space==other.m_Space; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  
  Self& operator=(const Self & other)
    {
    Superclass::operator=(other);
    if (this!=&other)
      {
      m_Offset = other.m_Offset;
      m_ChunkSize = other.m_ChunkSize;
      m_Space = other.m_Space;
      }
    return *this;
    }

  inline TOutput operator()( const TVector & A ) const 
  {
  int inputSize = A.Size();
  TOutput out(GetOutputDimension(inputSize));
  int tIndex = 0;
  for (int t = m_Offset; t < inputSize; t++)
    {
    if (t >= m_Offset && (t - m_Offset) % (m_ChunkSize + m_Space) < m_ChunkSize)
      {
      out[tIndex] = A[t];
      tIndex++;
      }
    }
  return out;
  }
  
  int GetOutputDimension(const int inputSize) const
    {
    int outputSize = (inputSize - m_Offset) / (m_ChunkSize + m_Space) * m_ChunkSize;
    int residual = (inputSize - m_Offset) % (m_ChunkSize + m_Space);
    outputSize += std::min(m_ChunkSize, residual);
    return outputSize;
    } 

  void SetArguments(const std::vector<int>& vec)
    {
    utlException(vec.size()!=3, "wrong size of vec");
    m_Offset = vec[0];
    m_ChunkSize = vec[1];
    m_Space = vec[2];
    }

  std::vector<int> GetArguments()
    {
    std::vector<int> vec(3);
    vec[0] = m_Offset;
    vec[1] = m_ChunkSize;
    vec[2] = m_Space;
    return vec;
    }

  utlSetGetMacro(Offset, int);
  utlSetGetMacro(ChunkSize, int);
  utlSetGetMacro(Space, int);

  void Print(std::ostream & os=std::cout) const
    {
    Superclass::Print(os);
    utlLogOSVar(os, m_Offset, m_ChunkSize, m_Space);
    }

protected:
  int m_Offset;
  int m_ChunkSize;
  int m_Space;
};


template< class TVector, class TOutput=TVector >  
class DotProduct : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef DotProduct Self;
  bool operator==(const Self & other) const  
  { return true; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  void operator=(const Self & other)   {}    

  inline TOutput operator()( const TVector & A, const TVector & B ) const 
  {
  int inputSize = A.Size();
  TOutput out(GetOutputDimension(inputSize));
  out[0]=0;
  for ( int i = 0; i < inputSize; ++i ) 
    out[0] += A[i]*B[i];
  return out;
  }
  
  inline TOutput operator()(const std::vector<TVector>& vec) const 
  {
  utlException(vec.size()!=2, "need to set only 2 vectors");
  TVector A = vec[0];
  TVector B = vec[1];
  int inputSize = A.Size();
  utlSAException(B.Size()!=inputSize)(inputSize)(B.Size()).msg("two vectors have different sizes");
  TOutput out(GetOutputDimension(inputSize));
  out[0]=0;
  for ( int i = 0; i < inputSize; ++i ) 
    out[0] += A[i]*B[i];
  return out;
  }
  
  int GetOutputDimension(const int inputSize) const
    {
    return 1;
    } 
  int GetOutputDimension(const std::vector<int>& sizeVec) const
    {
    return 1;
    } 
};

template< class TVector, class TOutput=TVector >  
class Compose : public VectorFunctorBase<TVector, TOutput>
{ 
public:
  typedef Compose Self;
  bool operator==(const Self & other) const  
  { return true; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 
  void operator=(const Self & other)   
    {}    

  inline TOutput operator()(const std::vector<TVector>& vec) const 
  {
  std::vector<int> sizeVec;
  for ( int i = 0; i < vec.size(); ++i ) 
    sizeVec.push_back(vec[i].Size());
  int outDim = GetOutputDimension(sizeVec);

  TOutput out(outDim);
  int kk=0;
  for ( int i = 0; i < vec.size(); ++i ) 
    {
    for ( int j = 0; j < sizeVec[i]; ++j ) 
      {
      out[kk] = vec[i][j];
      kk++;
      }
    }
  return out;
  }
  
  int GetOutputDimension(const std::vector<int>& sizeVec) const
    {
    return std::accumulate(sizeVec.begin(), sizeVec.end(), 0);
    } 

};

}

/** @}  */


}

#endif 
