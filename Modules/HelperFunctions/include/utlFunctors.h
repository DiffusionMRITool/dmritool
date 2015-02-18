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
#include "utlFunctors.h"

namespace utl
{

template<typename T>
struct Maximum
{
    UTL_ALWAYS_INLINE  T operator()(const T a, const T b) const
      { return a >= b ? a : b; }
};


template<typename T>
struct Minimum
{
    UTL_ALWAYS_INLINE  T operator()(const T a, const T b) const
      { return a <= b ? a : b; }
};

template<typename T>
struct Absolute
{
    UTL_ALWAYS_INLINE  T operator()(const T a) const
      { return std::abs(a); }
};


}

#endif 
