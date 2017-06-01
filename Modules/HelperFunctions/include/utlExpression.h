/**
 *       @file  utlExpression.h
 *      @brief  
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlExpression_h
#define __utlExpression_h

#include "utlTypeinfo.h"
#include "utlCoreMacro.h"
#include "utlSmartAssert.h"
#include "utlSTDHeaders.h"

namespace utl
{

/** 
* @addtogroup utlNDArray
*  @{
* */


template<class T1, class T2>
using SuperType = typename SuperFloatType<T1, T2>::type;

template<class TExpr>
using Expr1ValueType = typename TExpr::ValueType;

template<class TLeft, class TRight>
using Expr2ValueType = typename SuperFloatType<typename TLeft::ValueType, typename TRight::ValueType>::type;

template< typename SubType, typename ValueT> class Expr;
template<typename ValueT> class ScalarExprBase;

/** 
 * \typedef ScalarExpr 
 * \brief ScalarExprBase<double>
* */
typedef  ScalarExprBase<double>  ScalarExpr; 
/** 
 * \typedef ScalarComplexExpr
 * \brief ScalarExprBase<std::complex<double> >
 *  */
typedef  ScalarExprBase<std::complex<double> >  ScalarComplexExpr; 

/** 
 *
 *
 *  \class Expr
 *  \brief base class for expression
 *
 *  \tparam SubType inheritated class must put their type into this parameter
 *
 *  \ingroup utlNDArray
 */
template< typename SubType, typename ValueT>
class Expr
{
public:
  typedef ValueT                  ValueType;
  typedef unsigned int            SizeType;
  typedef SizeType const*         ShapeType;

  /** return reference of subtype instance of current class */
  UTL_ALWAYS_INLINE const SubType& ConstRef( void ) const
    {
    return *static_cast<const SubType*>(this);
    }
  /** return reference of subtype instance of current class */
  UTL_ALWAYS_INLINE SubType& Ref( void )
    {
    return *static_cast<SubType*>(this);
    }
  UTL_ALWAYS_INLINE static SizeType  GetDimension()
    { return SubType::GetDimension(); }

  UTL_ALWAYS_INLINE const ShapeType GetShape() const
    {
    return ConstRef().GetShape();
    }
};

/** \class ScalarExprBase
 *  \brief scalar expression base class for double and std::complex<double>
 *
 *  All scalar values will be converted to double or std::complex<double>.
 *
 *  \ingroup utlNDArray
 *  */
template<typename ValueT>
class ScalarExprBase: public Expr<ScalarExprBase<ValueT>, ValueT>
{
public:
  typedef Expr<ScalarExprBase<ValueT>, ValueT>      Superclass; 
  typedef ValueT                                ValueType;
  typedef typename Superclass::SizeType         SizeType;
  typedef typename Superclass::ShapeType        ShapeType;
  
  /** ScalarExpr has 0 dimension  */
  enum { Dimension = 0 };

  ValueT  m_Scalar;

  template<typename TValue>
  ScalarExprBase( TValue scalar ) : m_Scalar(static_cast<ValueType>(scalar)){
  }

  UTL_ALWAYS_INLINE static SizeType  GetDimension()
    { return 0; }
  UTL_ALWAYS_INLINE const ShapeType GetShape() const
    { return NULL; }
  UTL_ALWAYS_INLINE ValueType Eval( int i ) const
    {
    return m_Scalar;
    }
  /** return value  */
  UTL_ALWAYS_INLINE const ValueType ConstRef( void ) const
    {
    return m_Scalar;
    }
  /** return value  */
  UTL_ALWAYS_INLINE ValueType Ref( void )
    {
    return m_Scalar;
    }
};


/**
 *   \class   BinaryOpExpr 
 *   \brief   Binary operator expression class
 *   \ingroup utlNDArray
 */
template<typename OP,typename TLeft, typename TRight>
class BinaryOpExpr: public Expr<BinaryOpExpr<OP,TLeft,TRight>, typename SuperFloatType<typename TLeft::ValueType, typename TRight::ValueType>::type >
{
public:
    typedef Expr<BinaryOpExpr<OP,TLeft,TRight>, typename SuperFloatType<typename TLeft::ValueType, typename TRight::ValueType>::type >        Superclass;
    typedef typename Superclass::ValueType        ValueType;
    typedef typename Superclass::SizeType         SizeType;
    typedef typename Superclass::ShapeType        ShapeType;

    const TLeft& m_Left;
    const TRight& m_Right;
    OP m_OP;

    /** the two inputs should have the same shape, or at least one is a object of ScalarExpr  */
    BinaryOpExpr(const TLeft& lhs, const TRight& rhs): m_Left(lhs),m_Right(rhs), m_OP(OP())
    {
#if UTL_VERBOSITY>1
    SizeType lDim = TLeft::GetDimension();
    SizeType rDim = TRight::GetDimension();
    if (lDim>0 && rDim>0)
      {
      utlSAException(lDim>0 && rDim>0 && lDim!=rDim)(lDim)(rDim).msg("the dimensions of two sides are different");
      const ShapeType leftShape=m_Left.GetShape(), rightShape=m_Right.GetShape();
      for ( int i = 0; i < lDim; ++i ) 
        utlSAException(leftShape[i]!=rightShape[i])
          (i)(leftShape[i])(rightShape[i]).msg("the sizes of two sides are different");
      }
#endif
    }

    BinaryOpExpr(const TLeft&& lhs, const TRight&& rhs): m_Left(lhs),m_Right(rhs), m_OP(OP())
    {
#if UTL_VERBOSITY>1
    SizeType lDim = TLeft::GetDimension();
    SizeType rDim = TRight::GetDimension();
    if (lDim>0 && rDim>0)
      {
      utlSAException(lDim>0 && rDim>0 && lDim!=rDim)(lDim)(rDim).msg("the dimensions of two sides are different");
      const ShapeType leftShape=m_Left.GetShape(), rightShape=m_Right.GetShape();
      for ( int i = 0; i < lDim; ++i ) 
        utlSAException(leftShape[i]!=rightShape[i])
          (i)(leftShape[i])(rightShape[i]).msg("the sizes of two sides are different");
      }
#endif
    }
  
  BinaryOpExpr(const BinaryOpExpr&& other) : m_Left(other.m_Left),m_Right(other.m_Right), m_OP(OP())
    {
    }

    /** Determined by the non-scalar expression  */
  UTL_ALWAYS_INLINE static SizeType  GetDimension() 
    { 
    SizeType lDim = TLeft::GetDimension();
    SizeType rDim = TRight::GetDimension();
    utlSAException(lDim>0 && rDim>0 && lDim!=rDim)(lDim)(rDim).msg("the dimensions of two sides are different");
    if (lDim==0)
      return rDim;
    if (rDim==0)
      return lDim;
    return lDim;
    }
    /** Determined by the non-scalar expression  */
    UTL_ALWAYS_INLINE const ShapeType GetShape() const
      { 
      if (TLeft::GetDimension()==0)
        return m_Right.GetShape();
      if (TRight::GetDimension()==0)
        return m_Left.GetShape();
      return m_Left.GetShape();
      }
    UTL_ALWAYS_INLINE ValueType Eval( int i ) const
      {
      return m_OP( m_Left.Eval(i), m_Right.Eval(i) );
      }
};

template<typename OP,typename TLeft, typename ValueT>
class BinaryOpExpr<OP, TLeft, ScalarExprBase<ValueT>>: public Expr<BinaryOpExpr<OP,TLeft,ScalarExprBase<ValueT>>, Expr2ValueType<TLeft, ScalarExprBase<ValueT>> >
{
public:
    typedef Expr<BinaryOpExpr<OP,TLeft,ScalarExprBase<ValueT>>, Expr2ValueType<TLeft, ScalarExprBase<ValueT>> >        Superclass;
    typedef typename Superclass::ValueType        ValueType;
    typedef typename Superclass::SizeType         SizeType;
    typedef typename Superclass::ShapeType        ShapeType;

    const TLeft& m_Left;
    ValueType m_Scalar;
    OP m_OP;

    /** the two inputs should have the same shape, or at least one is a object of ScalarExpr  */
    BinaryOpExpr(const TLeft& lhs, const ScalarExprBase<ValueT>& rhs): m_Left(lhs), m_Scalar(rhs.m_Scalar), m_OP(OP())
    {
    }

    /** Determined by the non-scalar expression  */
  UTL_ALWAYS_INLINE static SizeType  GetDimension() 
    { 
    return TLeft::GetDimension();
    }
    /** Determined by the non-scalar expression  */
    UTL_ALWAYS_INLINE const ShapeType GetShape() const
      { 
      return m_Left.GetShape();
      }
    UTL_ALWAYS_INLINE ValueType Eval( int i ) const
      {
      return m_OP( m_Left.Eval(i), m_Scalar );
      }
};

template<typename OP,typename TRight, typename ValueT>
class BinaryOpExpr<OP, ScalarExprBase<ValueT>, TRight>: public Expr<BinaryOpExpr<OP,ScalarExprBase<ValueT>, TRight>, Expr2ValueType<TRight, ScalarExprBase<ValueT>> >
{
public:
    typedef Expr<BinaryOpExpr<OP,ScalarExprBase<ValueT>, TRight>, Expr2ValueType<ScalarExprBase<ValueT>, TRight> >        Superclass;
    typedef typename Superclass::ValueType        ValueType;
    typedef typename Superclass::SizeType         SizeType;
    typedef typename Superclass::ShapeType        ShapeType;

    const TRight& m_Right;
    ValueType m_Scalar;
    OP m_OP;

    /** the two inputs should have the same shape, or at least one is a object of ScalarExpr  */
    BinaryOpExpr(const ScalarExprBase<ValueT>& lhs, const TRight& rhs): m_Right(rhs), m_Scalar(lhs.m_Scalar), m_OP(OP())
    {
    }

    /** Determined by the non-scalar expression  */
  UTL_ALWAYS_INLINE static SizeType  GetDimension() 
    { 
    return TRight::GetDimension();
    }
    /** Determined by the non-scalar expression  */
    UTL_ALWAYS_INLINE const ShapeType GetShape() const
      { 
      return m_Right.GetShape();
      }
    UTL_ALWAYS_INLINE ValueType Eval( int i ) const
      {
      return m_OP( m_Scalar, m_Right.Eval(i) );
      }
};

template<typename OP,typename ValueT1, typename ValueT2>
class BinaryOpExpr<OP, ScalarExprBase<ValueT1>, ScalarExprBase<ValueT2>>: public Expr<BinaryOpExpr<OP,ScalarExprBase<ValueT1>, ScalarExprBase<ValueT2>>, SuperType<ValueT1, ValueT2> >
{
public:
    typedef Expr<BinaryOpExpr<OP,ScalarExprBase<ValueT1>, ScalarExprBase<ValueT2>>, SuperType<ValueT1, ValueT2> >        Superclass;
    typedef typename Superclass::ValueType        ValueType;
    typedef typename Superclass::SizeType         SizeType;
    typedef typename Superclass::ShapeType        ShapeType;

    ValueType m_Scalar1;
    ValueType m_Scalar2;
    OP m_OP;

    /** the two inputs should have the same shape, or at least one is a object of ScalarExpr  */
    BinaryOpExpr(const ScalarExprBase<ValueT1>& lhs, const ScalarExprBase<ValueT2>& rhs): m_Scalar1(lhs.m_Scalar), m_Scalar2(rhs.m_Scalar), m_OP(OP())
    {
    }

    /** Determined by the non-scalar expression  */
  UTL_ALWAYS_INLINE static SizeType  GetDimension() 
    { 
    return 0;
    }
    /** Determined by the non-scalar expression  */
    UTL_ALWAYS_INLINE const ShapeType GetShape() const
      { 
      return NULL;
      }
    UTL_ALWAYS_INLINE ValueType Eval( int i ) const
      {
      return m_OP( m_Scalar1, m_Scalar2 );
      }
};

/** Generate a BinaryOpExpr object  
 *  \ingroup utlNDArray
 * */
template<typename OP,typename TA, typename TB>
inline BinaryOpExpr<OP,TA,TB> 
MakeExpr( const Expr<TA, typename TA::ValueType> &lhs, const Expr<TB, typename TB::ValueType> &rhs )
{
  return BinaryOpExpr<OP,TA,TB>( lhs.ConstRef(), rhs.ConstRef() );
}

/** Generate a BinaryOpExpr object from two general expressions 
 *  \ingroup utlNDArray
 * */
template<typename OP,typename TLeft, typename TRight>
inline BinaryOpExpr<OP,TLeft,TRight> 
F(const Expr<TLeft, typename TLeft::ValueType>& lhs, const Expr<TRight, typename TRight::ValueType>& rhs)
{
  return BinaryOpExpr< OP,TLeft,TRight>(lhs.ConstRef(), rhs.ConstRef());
}

/** Generate a BinaryOpExpr object  */
template<typename OP,typename TA>
inline BinaryOpExpr<OP,TA,ScalarExpr > 
F( const Expr<TA, typename TA::ValueType> &lhs, const ScalarExpr &rhs ){
  return MakeExpr<OP>( lhs, rhs );
}
/** Generate a BinaryOpExpr object  */
template< typename OP,typename TB>
inline BinaryOpExpr<OP,ScalarExpr,TB> 
F( const ScalarExpr &lhs, const Expr<TB, typename TB::ValueType>& rhs ){
  return MakeExpr<OP>( lhs, rhs );
}
template<typename OP,typename TA>
inline BinaryOpExpr<OP,TA,ScalarComplexExpr > 
F( const Expr<TA, typename TA::ValueType> &lhs, const ScalarComplexExpr &rhs ){
  return MakeExpr<OP>( lhs, rhs );
}
/** Generate a BinaryOpExpr object  */
template< typename OP,typename TB>
inline BinaryOpExpr<OP,ScalarComplexExpr,TB> 
F( const ScalarComplexExpr &lhs, const Expr<TB, typename TB::ValueType>& rhs ){
  return MakeExpr<OP>( lhs, rhs );
}


/** \def __BinaryOpExpr_Op(name, op)
 *  \brief define expressions for binary math operators 
 *  \ingroup utlNDArray
 *  */
#define __BinaryOpExpr_Op(name, op)                                                                                            \
template<typename TLeft, typename TRight>                                                                                      \
inline BinaryOpExpr<name<Expr2ValueType<TLeft,TRight> >,TLeft,TRight>                                                          \
operator op (const Expr<TLeft, typename TLeft::ValueType>& lhs, const Expr<TRight, typename TRight::ValueType>& rhs){          \
    return MakeExpr<name<Expr2ValueType<TLeft,TRight> > >(lhs, rhs);                                                           \
}                                                                                                                              \
template<typename TRight>                                                                                                      \
inline BinaryOpExpr<name<Expr2ValueType<ScalarExpr,TRight> >,ScalarExpr,TRight>                                                \
operator op (const ScalarExpr &lhs, const Expr<TRight, typename TRight::ValueType>& rhs){                                      \
    return MakeExpr<name<Expr2ValueType<ScalarExpr,TRight> > >(lhs, rhs);                                                      \
}                                                                                                                              \
template<typename TLeft>                                                                                                       \
inline BinaryOpExpr<name<Expr2ValueType<TLeft,ScalarExpr> >, TLeft,ScalarExpr >                                                \
operator op (const Expr<TLeft, typename TLeft::ValueType>& lhs, const ScalarExpr& rhs){                                        \
    return MakeExpr<name<Expr2ValueType<TLeft,ScalarExpr> > >(lhs, rhs);                                                       \
}                                                                                                                              \
template<typename TRight>                                                                                                      \
inline BinaryOpExpr<name<Expr2ValueType<ScalarExpr,TRight> >,ScalarExpr,TRight>                                                \
operator op (const double lhs, const Expr<TRight, typename TRight::ValueType>& rhs){                                           \
    return MakeExpr<name<Expr2ValueType<ScalarExpr,TRight> > >(ScalarExpr(lhs), rhs);                                          \
}                                                                                                                              \
template<typename TLeft>                                                                                                       \
inline BinaryOpExpr<name<Expr2ValueType<TLeft,ScalarExpr> >, TLeft,ScalarExpr >                                                \
operator op (const Expr<TLeft, typename TLeft::ValueType>& lhs, const double rhs){                                             \
    return MakeExpr<name<Expr2ValueType<TLeft,ScalarExpr> > >(lhs, ScalarExpr(rhs));                                           \
}                                                                                                                              \
template<typename TRight>                                                                                                      \
inline BinaryOpExpr<name<Expr2ValueType<ScalarComplexExpr,TRight> >,ScalarComplexExpr,TRight>                                  \
operator op (const std::complex<double> lhs, const Expr<TRight, typename TRight::ValueType>& rhs){                             \
    return MakeExpr<name<Expr2ValueType<ScalarComplexExpr,TRight> > >(ScalarComplexExpr(lhs), rhs);                            \
}                                                                                                                              \
template<typename TLeft>                                                                                                       \
inline BinaryOpExpr<name<Expr2ValueType<TLeft,ScalarComplexExpr> >, TLeft,ScalarComplexExpr >                                  \
operator op (const Expr<TLeft, typename TLeft::ValueType>& lhs, const std::complex<double> rhs){                               \
    return MakeExpr<name<Expr2ValueType<TLeft,ScalarComplexExpr> > >(lhs, ScalarComplexExpr(rhs));                             \
}                                                                                                                              \


__BinaryOpExpr_Op(std::plus,        +)
__BinaryOpExpr_Op(std::minus,       -)
__BinaryOpExpr_Op(std::multiplies,  %)
__BinaryOpExpr_Op(std::divides,     /)


/**
 *   \class   UnaryOpExpr
 *   \brief   unary operator expression
 *   \ingroup utlNDArray
 */
template<typename OP, typename EType >
class UnaryOpExpr: public Expr< UnaryOpExpr<OP,EType>,  typename EType::ValueType >
{
public:
  typedef Expr<UnaryOpExpr<OP,EType>, typename EType::ValueType >          Superclass;
  typedef typename Superclass::ValueType        ValueType;
  typedef typename Superclass::SizeType         SizeType;
  typedef typename Superclass::ShapeType        ShapeType;

  const EType&  m_Expr;
  OP m_OP;
  
  UnaryOpExpr( const EType & expr):m_Expr(expr), m_OP(){}
  
  UTL_ALWAYS_INLINE static SizeType  GetDimension() 
    { 
    return EType::GetDimension();
    }
  UTL_ALWAYS_INLINE const ShapeType GetShape() const
    { 
    return m_Expr.GetShape();
    }
  UTL_ALWAYS_INLINE ValueType Eval( int i ) const
    {
    return m_OP( m_Expr.Eval(i) );
    }
};

/** make unary expression */
template<typename OP,typename EType>
inline UnaryOpExpr<OP,EType > MakeExpr( const Expr<EType, typename EType::ValueType> &expr )
{
  return UnaryOpExpr<OP,EType >( expr.ConstRef() );
}
/** make unary expression */
template<typename OP,typename EType>
inline UnaryOpExpr<OP,EType> F( const Expr<EType, typename EType::ValueType> &expr ){
  return MakeExpr<OP>(expr);
}

/** \def __UnaryOpExpr_Op(name, op)
 *  \brief define expressions for unary math operators */
#define __UnaryOpExpr_Op(name, op)                                        \
template<typename TExpr>                                                  \
inline UnaryOpExpr<name<typename TExpr::ValueType>,TExpr>                 \
operator op (const Expr<TExpr, typename TExpr::ValueType>& lhs){          \
    return MakeExpr<name<typename TExpr::ValueType>>(lhs);                \
}                                                                         \


__UnaryOpExpr_Op(std::negate, -)

/** @}  */

}


#endif 
