/**
 *       @file  utlCore11.h
 *      @brief  utl functions using c++11
 *     Created  "07-02-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlCore11_h
#define __utlCore11_h

#include "utlCoreMacro.h"
#include "utlTypeinfo.h"

/** @addtogroup utlHelperFunctions
@{ */

namespace utl 
{

template<std::size_t> 
struct Int_{};

template <class Tuple, size_t Pos>
inline std::ostream& 
PrintTuple(std::ostream& os, const Tuple& t, Int_<Pos> ) 
{
  os << std::get< std::tuple_size<Tuple>::value-Pos >(t) << ", ";
  return PrintTuple(os, t, Int_<Pos-1>());
}

/** print tuple  */
template <class Tuple>
inline std::ostream& 
PrintTuple(std::ostream& os, const Tuple& t, Int_<1> ) 
{
  return os << std::get<std::tuple_size<Tuple>::value-1>(t);
}


template <class... Args>
void
PrintTuple(const std::tuple<Args...>& t, const std::string& str="", std::ostream& os=std::cout)
{
  os << (str==""?"tuple":str) << " = " << t << std::endl;
}

/** Get number of arguments.  */
template<typename... Args>
inline auto 
GetNumberOfArgs(Args... args) -> decltype(sizeof...(args))
{
  return sizeof...(args);
}

// template <typename T>
// inline void 
// __PrintOS(std::ostream& os, const char* seperate, const T& t) 
// {
//     os  << t;
// }

// template<typename T, typename... Args>
// inline void 
// __PrintOS(std::ostream& os, const char* separate, const T& t, Args... args) // recursive variadic function
// {
//     os  << t << separate;
//     __PrintOS(os, separate, args...) ;
// }



/** Print variable number of arguments.  */
template<typename... Args>
inline void 
PrintOS(std::ostream& os, Args... args) 
{
  auto list = std::make_tuple(args...);
  os << list  << std::endl << std::flush;
}

/** Print variable number of arguments.  */
template<typename... Args>
inline void 
Print(Args... args) 
{
  PrintOS(std::cout, args...);
}


}

namespace std
{

/** print tuple  */
template <class... Args>
inline std::ostream& 
operator<<(std::ostream& os, const std::tuple<Args...>& t) 
{
  os << std::boolalpha; 
  int nn = sizeof...(Args);
  if (nn>1)
    os << '('; 
  utl::PrintTuple(os, t, utl::Int_<sizeof...(Args)>()); 
  if (nn>1)
    os << ')';
  os <<  std::noboolalpha << std::flush;
  return os;
}

}


/** @} */

#endif 
