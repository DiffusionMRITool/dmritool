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

#include "utlCoreMacro.h"
#include "utlSmartAssert.h"
#include "utlFunctors.h"
#include "utlSTDHeaders.h"

namespace utl
{


/** 
 *  \defgroup utlNDArray 
 *  \brief utl::NDArray with supports for expression template, blas, mkl
 *
 *  @{
 *
 *  \class Expr
 *  \brief base class for expression
 *
 *  \tparam SubType inheritated class must put their type into this parameter
 *
 *  \ingroup utlNDArray
 */
template< typename SubType>
class Expr
{
public:
  typedef double                  ValueType;
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

/** \class ScalarExpr
 *  \brief scalar expression 
 *
 *  All scalar values will be converted to double type.
 *
 *  \ingroup utlNDArray
 *  */
class ScalarExpr: public Expr<ScalarExpr>
{
public:
  typedef Expr<ScalarExpr>             Superclass; 
  typedef Superclass::ValueType        ValueType;
  typedef Superclass::SizeType         SizeType;
  typedef Superclass::ShapeType        ShapeType;
  
  /** ScalarExpr has 0 dimension  */
  enum { Dimension = 0 };

  double  m_Scalar;

  template<typename TValue>
  ScalarExpr( TValue scalar ) : m_Scalar(static_cast<ValueType>(scalar)){
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
class BinaryOpExpr: public Expr<BinaryOpExpr<OP,TLeft,TRight> >
{
public:
    typedef Expr<BinaryOpExpr<OP,TLeft,TRight> >        Superclass;
    typedef typename Superclass::ValueType        ValueType;
    typedef typename Superclass::SizeType         SizeType;
    typedef typename Superclass::ShapeType        ShapeType;

    const TLeft& m_Left;
    const TRight& m_Right;
    OP m_OP;

    /** the two inputs should have the same shape, or at least one is a object of ScalarExpr  */
    BinaryOpExpr(const TLeft& lhs, const TRight& rhs): m_Left(lhs),m_Right(rhs), m_OP()
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

/** Generate a BinaryOpExpr object  
 *  \ingroup utlNDArray
 * */
template<typename OP,typename TA, typename TB>
inline BinaryOpExpr<OP,TA,TB> 
MakeExpr( const Expr<TA> &lhs, const Expr<TB> &rhs )
{
  return BinaryOpExpr<OP,TA,TB>( lhs.ConstRef(), rhs.ConstRef() );
}

/** Generate a BinaryOpExpr object from two general expressions 
 *  \ingroup utlNDArray
 * */
template<typename OP,typename TLeft, typename TRight>
inline BinaryOpExpr<OP,TLeft,TRight> 
F(const Expr<TLeft>& lhs, const Expr<TRight>& rhs)
{
  return BinaryOpExpr< OP,TLeft,TRight>(lhs.ConstRef(), rhs.ConstRef());
}

/** Generate a BinaryOpExpr object  */
template<typename OP,typename TA>
inline BinaryOpExpr<OP,TA,ScalarExpr > 
F( const Expr<TA> &lhs, const ScalarExpr &rhs ){
  return MakeExpr<OP>( lhs, rhs );
}
/** Generate a BinaryOpExpr object  */
template< typename OP,typename TB>
inline BinaryOpExpr<OP,ScalarExpr,TB> 
F( const ScalarExpr &lhs, const Expr<TB>& rhs ){
  return MakeExpr<OP>( lhs, rhs );
}


/** \def __BinaryOpExpr_Op(name, op)
 *  \brief define expressions for binary math operators 
 *  \ingroup utlNDArray
 *  */
#define __BinaryOpExpr_Op(name, op)                                      \
template<typename TLeft, typename TRight>                                \
inline BinaryOpExpr<name,TLeft,TRight>                                   \
operator op (const Expr<TLeft>& lhs, const Expr<TRight>& rhs){           \
    return MakeExpr<name>(lhs, rhs);                                     \
}                                                                        \
template<typename TRight>                                                \
inline BinaryOpExpr<name,ScalarExpr,TRight>                              \
operator op (const ScalarExpr &lhs, const Expr<TRight>& rhs){            \
    return MakeExpr<name>(lhs, rhs);                                     \
}                                                                        \
template<typename TLeft>                                                 \
inline BinaryOpExpr<name,TLeft,ScalarExpr >                              \
operator op (const Expr<TLeft>& lhs, const ScalarExpr& rhs){             \
    return MakeExpr<name>(lhs, rhs);                                     \
}                                                                        \

__BinaryOpExpr_Op(std::plus<double>,        +)
__BinaryOpExpr_Op(std::minus<double>,       -)
__BinaryOpExpr_Op(std::multiplies<double>,  %)
__BinaryOpExpr_Op(std::divides<double>,     /)


/**
 *   \class   UnaryOpExpr
 *   \brief   unary operator expression
 *   \ingroup utlNDArray
 */
template<typename OP, typename EType >
class UnaryOpExpr: public Expr< UnaryOpExpr<OP,EType> >
{
public:
  typedef Expr<UnaryOpExpr<OP,EType> >          Superclass;
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
inline UnaryOpExpr<OP,EType > MakeExpr( const Expr<EType> &expr )
{
  return UnaryOpExpr<OP,EType >( expr.ConstRef() );
}
/** make unary expression */
template<typename OP,typename EType>
inline UnaryOpExpr<OP,EType> F( const Expr<EType> &expr ){
  return MakeExpr<OP>(expr);
}

/** \def __UnaryOpExpr_Op(name, op)
 *  \brief define expressions for unary math operators */
#define __UnaryOpExpr_Op(name, op)                                       \
template<typename TExpr>                                                 \
inline UnaryOpExpr<name,TExpr>                                           \
operator op (const Expr<TExpr>& lhs){                                    \
    return MakeExpr<name>(lhs);                                          \
}                                                                        \

__UnaryOpExpr_Op(std::negate<double>, -)

/** @}  */

}


#endif 
