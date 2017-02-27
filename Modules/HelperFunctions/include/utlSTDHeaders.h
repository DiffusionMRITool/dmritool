/**
 *       @file  utlSTDHeaders.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-23-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlSTDHeaders_h
#define __utlSTDHeaders_h

/** @addtogroup utlHelperFunctions
@{ */


#if __cplusplus < 201103L

#include <tr1/memory>
#define utl_shared_ptr std::tr1::shared_ptr
#include <tr1/unordered_map>
#define utl_unordered_map std::tr1::unordered_map
#define utl_hash std::tr1::hash
  
#define __UTL_constexpr  const

/** http://stackoverflow.com/questions/23414270/c-complex-and-complex-h-in-the-same-file  */
extern "C" {
#include <complex.h>
#undef complex
}
#include <complex>


#else

#define __UTL_constexpr  constexpr

#include <memory>
#define utl_shared_ptr std::shared_ptr
#include <unordered_map>
#define utl_unordered_map std::unordered_map
#define utl_hash std::hash

#include <complex>
#include <type_traits>

#endif

// #include <functional>

#include <vector>

namespace utl
{

/** 
 * using boost::hash_combine template <class T> 
 * */
template <class T>
inline void 
hash_combine(std::size_t& seed, T const& v) 
{ 
  seed ^= utl_hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2); 
}

}

namespace std 
{ 

#if __cplusplus < 201103L
namespace tr1 
{
#endif

template<typename T> 
class hash<std::vector<T> > 
{ 
public:
typedef std::vector<T> argument_type; 
typedef std::size_t result_type; 
result_type operator()(argument_type const& in) const 
  { 
  size_t size = in.size(); 
  size_t seed = 0; 
  for (size_t i = 0; i < size; i++) //Combine the hash of the current vector with the hashes of the previous ones 
    utl::hash_combine(seed, in[i]); 
  return seed; 
  } 
}; 

#if __cplusplus < 201103L
}
#endif

}

    /** @} */

#endif 
