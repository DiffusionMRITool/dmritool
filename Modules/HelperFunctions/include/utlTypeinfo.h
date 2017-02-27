/**
 *       @file  utlTypeinfo.h
 *      @brief  
 *     Created  "07-04-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlTypeinfo_h
#define __utlTypeinfo_h


#include <type_traits>
#include <typeinfo>

#include <complex>
#include <iostream>
#include <memory>

#ifdef __GNUC__
#include <cstdlib>
#include <cxxabi.h>
#endif


namespace utl
{


template<bool B, class T, class F>
using conditional_t = typename ::std::conditional<B,T,F>::type;

template<class T>
using remove_const_t = typename ::std::remove_const<T>::type;

template<class T>
using remove_reference_t = typename ::std::remove_reference<T>::type;

template<class T>
using remove_reference_t = typename ::std::remove_reference<T>::type;

template<class... T>
using common_type_t = typename ::std::common_type<T...>::type;


/** 
 * \brief Remove complex:
 *
 * - std::complex<double> -> double
 * - std::complex<float> -> float
 * - other types: T -> T 
 *  */
template<typename T > 
struct remove_complex 
  { 
  typedef conditional_t
    <
    std::is_same<T, std::complex<double>>::value, double, 
    conditional_t<std::is_same<T, std::complex<float>>::value, float, T>
    > 
  type;
  };

template<class T>
using remove_complex_t = typename remove_complex<T>::type;

template<typename T, typename t> 
struct Superset 
  { 
  typedef conditional_t<sizeof(T) >= sizeof(t), T, t> type; 
  };

template<class T, typename t>
using superset_t = typename Superset<T,t>::type;

/** \brief float type which can be coverted to
 *
 * If it is scalar -> double. 
 * If it is complex -> std::complex<double> 
 * */
template<typename T, typename t> 
struct SuperFloatType
  { 
  static_assert(::std::is_scalar<T>::value 
    || ::std::is_same<T,std::complex<double>>::value
    || ::std::is_same<T,std::complex<float>>::value, "T must be a scalar or std::complex<double> or std::complex<float>");
  static_assert(::std::is_scalar<t>::value 
    || ::std::is_same<t,std::complex<double>>::value
    || ::std::is_same<t,std::complex<float>>::value, "t must be a scalar or std::complex<double> or std::complex<float>");
  typedef conditional_t<std::is_scalar<T>::value && std::is_scalar<t>::value, double, std::complex<double> > type; 
  };



#ifdef __GNUC__

/**
 * This implementation is adapted from the solution here:
 *  http://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname
 *
 * Note: this is also defined by clang and icc
 */
inline std::string Demangle(const char *name) 
{
    int status = -1;
    std::unique_ptr<char, void(*)(void*)> uptr { abi::__cxa_demangle(name, NULL, NULL, &status), std::free };
    return status == 0 ? uptr.get() : name;
}

#else
// TODO: support demangle for non-GCC compilers
inline std::string Demangle(const char *name) 
{
    return name;
}
#endif

template<class T>
inline std::string TypeName() 
{
    return Demangle(typeid(T).name());
}

template<class T>
inline std::string TypeName(const T&) 
{
    return TypeName<T>();
}


}


#endif 
