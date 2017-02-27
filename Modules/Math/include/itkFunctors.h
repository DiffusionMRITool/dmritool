/*=========================================================================
*
* Copyright Marius Staring, Stefan Klein, David Doria. 2011.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0.txt
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
*=========================================================================*/
#ifndef __itkFunctors_h
#define __itkFunctors_h

#include "vnl/vnl_math.h"
#include "vnl/vnl_erf.h"
#include "itkNumericTraits.h"


namespace itk {

namespace Functor {

/** @addtogroup Math
@{ */

/** Arithmetic functors which use m_Argument. */
#define __FunctorOneArgumentRight(Name, Op)                             \
template< class TInput, class TArgument=TInput, class TOutput=TInput >  \
class Name                                                              \
{                                                                       \
public:                                                                 \
  bool operator!=(const Name & other) const                             \
  { return m_Argument!=other.m_Argument; }                              \
  bool operator==(const Name & other) const                             \
  { return !( *this != other ); }                                       \
  inline TOutput operator()( const TInput & A ) const                   \
  {                                                                     \
    return static_cast<TOutput>( A Op this->m_Argument );               \
  }                                                                     \
  void SetArgument( TArgument arg ){ this->m_Argument = arg; }          \
  TArgument GetArgument( ) const { return this->m_Argument; }           \
private:                                                                \
  TArgument m_Argument;                                                 \
};                                                                      \


/** Arithmetic functors which use m_Argument. */
#define __FunctorOneArgumentLeft(Name, Op)                              \
template< class TInput, class TArgument=TInput, class TOutput=TInput >  \
class Name                                                              \
{                                                                       \
public:                                                                 \
  bool operator!=(const Name &) const                                   \
  { return m_Argument!=other.m_Argument; }                              \
  bool operator==(const Name & other) const                             \
  { return !( *this != other ); }                                       \
  inline TOutput operator()( const TInput & A ) const                   \
  {                                                                     \
    return static_cast<TOutput>( this->m_Argument Op A );               \
  }                                                                     \
  void SetArgument( TArgument arg ){ this->m_Argument = arg; };         \
  TArgument GetArgument( ) const { return this->m_Argument; }           \
private:                                                                \
  TArgument m_Argument;                                                 \
};                                                                      \

__FunctorOneArgumentRight(PLUS,   +)
__FunctorOneArgumentRight(RMINUS, -)
__FunctorOneArgumentRight(TIMES,  *)
__FunctorOneArgumentRight(RDIVIDE,/)
__FunctorOneArgumentRight(EQUAL, ==)

__FunctorOneArgumentRight(LMINUS, -)
__FunctorOneArgumentRight(LDIVIDE,/)


template< class TInput, class TArgument=TInput, class TOutput=TInput >
class NLOG
{
public:
  bool operator!=(const NLOG & other) const
  { return m_Argument!=other.m_Argument; } 
  bool operator==(const NLOG & other) const
  { return !( *this != other ); }
  inline TOutput operator()( const TInput & A ) const
  {
    return static_cast<TOutput>( vcl_log( static_cast<double>( A ) ) / vcl_log( static_cast<double>( this->m_Argument ) ) );
  }
  void SetArgument( TArgument arg ){ this->m_Argument = arg; };
  TArgument GetArgument( ) const { return this->m_Argument; } 
private:
  TArgument m_Argument;
};



/** In the following classes, the argument is always double */

template< class TInput, class TArgument=TInput, class TOutput=TInput >
class RPOWER
{
public:
  bool operator!=(const RPOWER & other) const
  { return m_Argument!=other.m_Argument; } 
  bool operator==(const RPOWER & other) const
  { return !( *this != other ); }
  inline TOutput operator()( const TInput & A ) const
  {
    return static_cast<TOutput>( vcl_pow( static_cast<double>( A ), static_cast<double>( this->m_Argument ) ) );
  }
  void SetArgument( TArgument arg ){ this->m_Argument = arg; };
  TArgument GetArgument( ) const { return this->m_Argument; } 
private:
  TArgument m_Argument;
};

template< class TInput, class TArgument=TInput, class TOutput=TInput >
class LPOWER
{
public:
  bool operator!=(const LPOWER & other) const
  { return m_Argument!=other.m_Argument; } 
  bool operator==(const LPOWER & other) const
  { return !( *this != other ); }
  inline TOutput operator()( const TInput & A ) const
  {
    return static_cast<TOutput>( vcl_pow( static_cast<double>( this->m_Argument ), static_cast<double>( A ) ) );
  }
  void SetArgument( TArgument arg ){ this->m_Argument = arg; };
  TArgument GetArgument( ) const { return this->m_Argument; } 
private:
  TArgument m_Argument;
};

// NOTE: use utl::Functor when there is no m_Argument

// [>* Funtions that don't use m_Argument. <]
// #define __FunctorNoArgument(Name, Func)                 \
// template< class TInput, class TOutput=TInput >          \
// class Name                                              \
// {                                                       \
// public:                                                 \
//   bool operator!=(const Name &) const                   \
//   { return false; }                                     \
//   bool operator==(const Name & other) const             \
//   { return !( *this != other ); }                       \
//   inline TOutput operator()( const TInput & A ) const   \
//   {                                                     \
//     return static_cast<TOutput>( Func(A) );             \
//   }                                                     \
// };                                                      \


// #define __FunctorNoArgument2(Name, Func, Func2)         \
// template< class TInput, class TOutput=TInput >          \
// class Name                                              \
// {                                                       \
// public:                                                 \
//   bool operator!=(const Name &) const                   \
//   { return false; }                                     \
//   bool operator==(const Name & other) const             \
//   { return !( *this != other ); }                       \
//   inline TOutput operator()( const TInput & A ) const   \
//   {                                                     \
//     return static_cast<TOutput>( Func(Func2(A)) );      \
//   }                                                     \
// };                                                      \

// __FunctorNoArgument(MATHSIGN,  vnl_math_sgn)
// __FunctorNoArgument(ABSOLUTE,   std::abs) // use utl::Functor::Abs

// __FunctorNoArgument2(FLOOR, std::floor,     static_cast<double>)
// __FunctorNoArgument2(CEIL,  std::ceil,      static_cast<double>)
// __FunctorNoArgument2(ROUND, vnl_math_rnd,   static_cast<double>)
// __FunctorNoArgument2(SQR,   vnl_math_sqr,   static_cast<double>)
// __FunctorNoArgument2(SQRT,  std::sqrt,      static_cast<double>)
// __FunctorNoArgument2(LN,    vcl_log,        static_cast<double>)
// __FunctorNoArgument2(LOG10, vcl_log10,      static_cast<double>)
// __FunctorNoArgument2(EXP,   std::exp,       static_cast<double>)
// __FunctorNoArgument2(SIN,   std::sin,       static_cast<double>)
// __FunctorNoArgument2(COS,   std::cos,       static_cast<double>)
// __FunctorNoArgument2(TAN,   std::tan,       static_cast<double>)
// __FunctorNoArgument2(ARCSIN,vcl_asin,       static_cast<double>)
// __FunctorNoArgument2(ARCCOS,vcl_acos,       static_cast<double>)
// __FunctorNoArgument2(ARCTAN,vcl_atan,       static_cast<double>)

// template< class TInput, class TOutput=TInput >
// class NEG
// {
// public:
//   bool operator!=(const NEG &) const
//   { return false; }
//   bool operator==(const NEG & other) const
//   { return !( *this != other ); }
//   inline TOutput operator()( const TInput & A ) const
//   {
//     return static_cast<TOutput>( -A ); 
//   }
// };

/** More complicated functors. */

template< class TInput, class TArgument=TInput, class TOutput=TInput >
class LINEAR
{
public:
  bool operator!=(const LINEAR & other) const
  { return m_Argument1!=other.m_Argument1 || m_Argument2!=other.m_Argument2; }
  bool operator==(const LINEAR & other) const
  { return !( *this != other ); }
  inline TOutput operator()( const TInput & A ) const
  {
    return static_cast<TOutput>( this->m_Argument1 * A + this->m_Argument2 );
  }
  void SetArgument1( TArgument arg ){ this->m_Argument1 = arg; };
  void SetArgument2( TArgument arg ){ this->m_Argument2 = arg; };
  TArgument GetArgument1( ) const { return this->m_Argument1; } 
  TArgument GetArgument2( ) const { return this->m_Argument2; } 
private:
  TArgument m_Argument1;
  TArgument m_Argument2;
};

    /** @} */

} // end namespace Functor


} // end namespace itk



#endif //#ifndef __itkFunctors_h
