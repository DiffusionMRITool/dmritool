/**
 *       @file  utlCoreMacro.h
 *      @brief  macros for utlCore
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "10-04-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlCoreMacro_h
#define __utlCoreMacro_h

/** @addtogroup utlHelperFunctions
@{ */
// Detect/configure OS variables.
//
// Define 'UTL_OS' to: '0' for an unknown OS (will try to minize library dependencies).
//                      '1' for a Unix-like OS (Linux, Solaris, BSD, MacOSX, Irix, ...).
//                      '2' for Microsoft Windows.
//                      (auto-detection is performed if 'UTL_OS' is not set by the user).
#ifndef UTL_OS
#if defined(unix)        || defined(__unix)      || defined(__unix__)       \
 || defined(linux)       || defined(__linux)     || defined(__linux__)      \
 || defined(sun)         || defined(__sun)                                  \
 || defined(BSD)         || defined(__OpenBSD__) || defined(__NetBSD__)     \
 || defined(__FreeBSD__) || defined __DragonFly__                           \
 || defined(sgi)         || defined(__sgi)                                  \
 || defined(__MACOSX__)  || defined(__APPLE__)                              \
 || defined(__CYGWIN__)                    
#define UTL_OS 1
#elif defined(_MSC_VER) || defined(WIN32)  || defined(_WIN32) || defined(__WIN32__) \
   || defined(WIN64)    || defined(_WIN64) || defined(__WIN64__)
#define UTL_OS 2
#else
#define UTL_OS 0
#endif
#elif !(UTL_OS==0 || UTL_OS==1 || UTL_OS==2)
#error UTL: Invalid configuration variable 'UTL_OS'.
#error (correct values are '0 = unknown OS', '1 = Unix-like OS', '2 = Microsoft Windows').
#endif

#define UTL_COMMA ,


#ifndef M_EPS
#define M_EPS 1e-9
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


/** 1/(4*pi^2)  */
#ifndef ONE_OVER_4_PI_2
#define ONE_OVER_4_PI_2  0.025330295910584442860969865802431910 
#endif

#if defined(__BORLANDC__)
  #define __utl_LOCATION__ __FUNC__
#elif defined(_WIN32) && !defined(__MINGW32__) && !defined(__CYGWIN__) && !defined(CSWIG)
  #define __utl_LOCATION__ __FUNCSIG__
#elif defined(__GNUC__)
  #define __utl_LOCATION__ __PRETTY_FUNCTION__
#else
  #define __utl_LOCATION__ __FUNCTION__
#endif

#define UTL_ALWAYS_INLINE inline __attribute__((always_inline))

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#if UTL_OS==1
#include <unistd.h>
#endif

#include "utlSTDHeaders.h"

// #include "utlSmartAssert.h"

enum {
  COLOR_NORMAL=0,
  COLOR_BOLD,
  COLOR_BLACK,
  COLOR_WHITE,
  COLOR_PURPLE,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW,
  COLOR_BLUE,
  COLOR_CYAN
};

enum {
 ROW_MAJOR=0,
 COLUMN_MAJOR 
};

namespace utl
{

inline std::string 
GetColoredString(const std::string& str, const int color)
{
#if UTL_OS==1
  static const std::string t_normal("\033[0;0;0m");
  static const std::string t_black("\033[0;30;59m");
  static const std::string t_white("\033[0;37;59m");
  static const std::string t_bold("\033[1m");
  static const std::string t_red("\033[1;31;59m");
  static const std::string t_green("\033[0;32;59m");
  static const std::string t_yellow("\033[0;33;59m");
  static const std::string t_blue("\033[0;34;59m");
  static const std::string t_purple("\033[0;35;59m");
  static const std::string t_cyan("\033[0;36;59m");
  bool isTerminal = isatty(fileno(stdout)) != 0;
#else
  static const bool isTerminal = false;
#endif

  if (!isTerminal)
    return str;

  switch ( color )
    {
    case COLOR_RED : return t_red+str+t_normal;
    case COLOR_BOLD :  return t_bold+str+t_normal;
    case COLOR_PURPLE :   return t_purple+str+t_normal;
    case COLOR_GREEN :   return t_green+str+t_normal;
    case COLOR_YELLOW :   return t_yellow+str+t_normal;
    case COLOR_CYAN :   return t_cyan+str+t_normal;
    case COLOR_BLUE :   return t_blue+str+t_normal;
    case COLOR_BLACK :   return t_black+str+t_normal;
    case COLOR_WHITE :   return t_white+str+t_normal;
    default :  return t_normal+str+t_normal;
    }
}

}


#if UTL_OS==1

#define __UTL_FATAL_STRING  utl::GetColoredString("Fatal", COLOR_RED)
#define __UTL_ERROR_STRING  utl::GetColoredString("Error", COLOR_RED)
#define __UTL_WARNING_STRING  utl::GetColoredString("Warning", COLOR_RED)
#define __UTL_DEBUG_STRING  utl::GetColoredString("Debug", COLOR_CYAN)
#define __UTL_LOG_STRING  utl::GetColoredString("Log", COLOR_BOLD)
#define __UTL_BOLD(str) utl::GetColoredString(str, COLOR_BOLD)
#define __UTL_EXPSTR(str) utl::GetColoredString(str, COLOR_GREEN)

// #define __UTL_FATAL_STRING  "\033[1;31;59mFatal\033[0;0;0m"
// #define __UTL_ERROR_STRING  "\033[1;31;59mError\033[0;0;0m"
// #define __UTL_WARNING_STRING  "\033[1;31;59mWarning\033[0;0;0m"
// #define __UTL_DEBUG_STRING  "\033[0;36;59mDebug\033[0;0;0m"
// #define __UTL_LOG_STRING  "\033[0;36;59mLog\033[0;0;0m"
// #define __UTL_BOLD(str)   std::string("\033[1m").append(str).append("\033[0;0;0m")
// #define __UTL_EXPSTR(str) std::string("\033[0;32;59m").append(str).append("\033[0;0;0m")

#else

#define __UTL_FATAL_STRING  "Fatal"
#define __UTL_ERROR_STRING  "Error"
#define __UTL_WARNING_STRING  "Warning"
#define __UTL_DEBUG_STRING  "Debug"
#define __UTL_LOG_STRING  "Log"
#define __UTL_BOLD(str) str
#define __UTL_EXPSTR(str) str

#endif


#define __UTL_Print_LOCATION  "In "<< __UTL_BOLD("File") << ": " <<__FILE__<< ", "<<__UTL_BOLD("Line")<<": " << __LINE__   << ", "<<__UTL_BOLD("Function")<<": " << __utl_LOCATION__ <<  "\n" 


enum{
  /// no log
  LOG_MUTE=0,   
  /// normal log
  LOG_NORMAL=1, 
  /// used for debug information. this->GetDebug()
  LOG_DEBUG=2,  
  /// log for large matrix or vectors. 
  LOG_LARGE=3,   
  /// for all possible logs. 
  LOG_ALL=100000000   
};

namespace utl
{
/** global Log verbosity level. Some class may have its own local LogLevel which override the global one.  */
static int  LogLevel=LOG_NORMAL; 

inline bool IsLogMute(const int level=utl::LogLevel)
  {
  return level<=LOG_MUTE;
  }
inline bool IsLogNormal(const int level=utl::LogLevel)
  {
  return LOG_NORMAL<=level;
  }
inline bool IsLogDebug(const int level=utl::LogLevel)
  {
  return LOG_DEBUG<=level;
  }
inline bool IsLogLarge(const int level=utl::LogLevel)
  {
  return LOG_LARGE<=level;
  }
inline bool IsLogAll(const int level=utl::LogLevel)
  {
  return LOG_ALL<=level;
  }
}


#if __cplusplus >= 201103L

#include "utlCore11.h"

// namespace utl
// {
// template<typename... Args>
// auto GetNumberOfArgs(Args... args) -> decltype(sizeof... (args));
// }

#define utlNumberOfArgs(...)   utl::GetNumberOfArgs(__VA_ARGS__) 


/** Conditional print variable number of arguments to os */
#define PrintVar(cond, os, ...)                                                                          \
do                                                                                                       \
{ if ((cond))                                                                                            \
    {                                                                                                    \
    if (utlNumberOfArgs(__VA_ARGS__)>1)                                                                  \
      { os << "(" << #__VA_ARGS__ << ") = ";  utl::PrintOS(os, __VA_ARGS__);   }                         \
    else                                                                                                 \
      { os << #__VA_ARGS__ << " = ";  utl::PrintOS(os, __VA_ARGS__);   }                                 \
    }                                                                                                    \
} while(0)


#define utlVLogOS_IF(cond, level, os)                                     \
do                                                                        \
  {                                                                       \
  if ((level<=utl::LogLevel) && (cond))                                   \
    { os << std::endl << std::flush; }                                    \
  } while (0)
  

#define utlVLogOS(level, os)                            utlVLogOS_IF(true, level, os)
#define utlVLog_IF(cond, level, expr)                   utlVLogOS_IF(cond, level, std::cout << expr)
#define utlVLog(level, expr)                            utlVLogOS_IF(true, level, std::cout << expr)

#define utlLogOS_IF(cond, os)                           utlVLogOS_IF(cond, 0, os)
#define utlLogOS(os)                                    utlVLogOS_IF(true, 0, os)
#define utlLog_IF(cond, expr)                           utlVLogOS_IF(cond, 0, std::cout << expr)
#define utlLog(expr)                                    utlVLogOS_IF(true, 0, std::cout << expr)


#define utlVLogOSVar_IF(cond, level, os, ...)                                                            \
do                                                                                                       \
{ if ( level<=utl::LogLevel )                                                                            \
     {  PrintVar(cond, os, __VA_ARGS__);  }                                                              \
} while(0)


#define utlVLogOSVar(level, os, ...)                    utlVLogOSVar_IF(true, level, os,        __VA_ARGS__)
#define utlVLogVar_IF(cond, level, ...)                 utlVLogOSVar_IF(cond, level, std::cout, __VA_ARGS__)
#define utlVLogVar(level, ...)                          utlVLogOSVar_IF(true, level, std::cout, __VA_ARGS__)


#define utlLogOSVar_IF(cond, os, ...)                   utlVLogOSVar_IF(cond, 0, os,        __VA_ARGS__)
#define utlLogOSVar(os, ...)                            utlVLogOSVar_IF(true, 0, os,        __VA_ARGS__)
#define utlLogVar_IF(cond, ...)                         utlVLogOSVar_IF(cond, 0, std::cout, __VA_ARGS__)
#define utlLogVar(...)                                  utlVLogOSVar_IF(true, 0, std::cout, __VA_ARGS__)


  
#define utlVLogOSPosition_IF(cond,level, os)                                                                                   \
do                                                                                                                             \
{  if ((level<=utl::LogLevel) && (cond))                                                                                       \
    { os << "\n"<<__UTL_BOLD("Work Flow")<<": "<< __UTL_Print_LOCATION << std::flush;                                          \
    }                                                                                                                          \
} while(0)

#define utlVLogPosition(level)   utlVLogOSPosition_IF(true, level, std::cout)    


#if UTL_VERBOSITY>0

#define utlOSPrintVar(cond, os, ...)                          PrintVar(cond,      os,   __VA_ARGS__)
#define utlPrintVar(cond, ...)                                PrintVar(cond, std::cout, __VA_ARGS__)

#else

#define utlOSPrintVar(cond, os, ...)                          
#define utlPrintVar(cond, ...)            


#endif


#endif 


template <typename T, size_t N>
char(&__utlArraySizeHelper(T(&array)[N]))[N];
#ifndef COMPILER_MSVC
template <typename T, size_t N>
char(&__utlArraySizeHelper(const T(&array)[N]))[N];
#endif
#define utlArraySize(array) (sizeof(__utlArraySizeHelper(array)))


// #define utlArraySize(arr) (sizeof(arr)/sizeof(*(arr)))


/**
 * \brief  macros which are used for debug with or without matlab mex files 
 * \author Jian Cheng
 */
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#define utlAbort(expout) do { mexErrMsgTxt(expout); } while(0)
#else
#define utlAbort(expout) do { std::cerr << expout <<"\n" << std::flush; abort(); } while(0)
#endif

/**
 * \brief  macros which are used for debug 
 * \author Jian Cheng
 */
#define __utlConditionFailPrint(cond)      __UTL_Print_LOCATION <<__UTL_BOLD("Expression")<<": '" << __UTL_EXPSTR(#cond) << "' failed. " << "\n"
#define __utlConditionSucceedPrint(cond)   __UTL_Print_LOCATION <<__UTL_BOLD("Expression")<<": '" << __UTL_EXPSTR(#cond) << "' satisfied. " << "\n"

#define utlSASetLog(log) Assert::set_log(log) 

#define utlSAGlobalPrintIf(expr)     \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert(#expr).log().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A 


#define utlSAGlobalPrint utlSAGlobalPrintIf("") 

#define utlSAGlobalWarning(expr)        \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert(#expr).warn().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A 


#define utlSAGlobalException(expr)        \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert(#expr).error().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A 


#define utlSAGlobalAssert(expr)        \
    if ( (expr) ) ; \
    else ::smart_assert::make_assert(#expr).error().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_assert).SMART_ASSERT_A 


#define utlGlobalException(cond,expout)                                                                                                                 \
do                                                                                                                                                      \
{                                                                                                                                                       \
  if ((cond))                                                                                                                                           \
    { std::cerr << "\n"<<__UTL_ERROR_STRING<<": " << __utlConditionSucceedPrint(cond) << __UTL_BOLD("msg")<<": '"<<expout << "'\n" << std::flush;       \
    utlAbort(""); }                                                                                                                                     \
} while (0)                                                                                               


#define utlGlobalAssert(cond,expout)                                                                                                                    \
do                                                                                                                                                      \
{                                                                                                                                                       \
  if (!(cond))                                                                                                                                          \
    { std::cerr << "\n"<<__UTL_ERROR_STRING<<": " << __utlConditionFailPrint(cond) << __UTL_BOLD("msg")<<": '"<<expout << "\n" << std::flush;           \
    utlAbort(""); }                                                                                                                                     \
} while (0)


#define utlOSGlobalWarning(cond,expout,os)                                                                                                              \
do                                                                                                                                                      \
{  if ((cond))                                                                                                                                          \
    { os << "\n"<<__UTL_WARNING_STRING<<": "<< __utlConditionSucceedPrint(cond) <<  __UTL_BOLD("msg")<<": '"<<expout << "\n" << std::flush;   }         \
} while(0)


#define utlGlobalWarning(cond,expout)                         utlOSGlobalWarning(cond,expout,std::cout)               


#define PrintEnum1(cond,var,val1,os)                                                                     \
do                                                                                                       \
{ if ((cond) && ((var)==(val1)))                                                                         \
    { os << #var <<" = " << #val1 << std::endl << std::flush; }                                          \
} while(0)


#define PrintEnum2(cond,var,val1,val2,os)                                                                \
do                                                                                                       \
{ PrintEnum1(cond,var,val1,os);                                                                          \
  if ((cond) && ((var)==(val2)))                                                                         \
    { os << #var <<" = " << #val2 << std::endl << std::flush; }                                          \
} while(0)


#define PrintEnum3(cond,var,val1,val2,val3,os)                                                           \
do                                                                                                       \
{ PrintEnum2(cond,var,val1,val2,os);                                                                     \
  if ((cond) && ((var)==(val3)))                                                                         \
    { os << #var <<" = " << #val3 << std::endl << std::flush; }                                          \
} while(0)


#define PrintEnum4(cond,var,val1,val2,val3,val4,os)                                                      \
do                                                                                                       \
{ PrintEnum3(cond,var,val1,val2,val3,os);                                                                \
  if ((cond) && ((var)==(val4)))                                                                         \
    { os << #var <<" = " << #val4 << std::endl << std::flush; }                                          \
} while(0)


#define PrintEnum5(cond,var,val1,val2,val3,val4,val5,os)                                                 \
do                                                                                                       \
{ PrintEnum4(cond,var,val1,val2,val3,val4,os);                                                           \
  if ((cond) && ((var)==(val5)))                                                                         \
    { os << #var <<" = " << #val5 << std::endl << std::flush; }                                          \
} while(0)


#define PrintEnum6(cond,var,val1,val2,val3,val4,val5,val6,os)                                            \
do                                                                                                       \
{ PrintEnum5(cond,var,val1,val2,val3,val4,val5,os);                                                      \
  if ((cond) && ((var)==(val6)))                                                                         \
    { os << #var <<" = " << #val6 << std::endl << std::flush; }                                          \
} while(0)


#define PrintVar1(cond,var,os)                                                                           \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << std::boolalpha << #var <<" = " <<  (var) << std::endl << std::flush << std::noboolalpha; }   \
} while(0)


#define PrintVar2(cond,var1,var2,os)                                                                     \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << std::boolalpha << "("<<#var1<<", "<< #var2<<") = (" <<  (var1) << ", " << (var2) << ")"      \
    << std::endl << std::flush << std::noboolalpha; }                                                    \
} while(0)


#define PrintVar3(cond,var1,var2,var3,os)                                                                                     \
do                                                                                                                            \
{ if ((cond))                                                                                                                 \
    { os << std::boolalpha << "("<<#var1<<", "<<#var2<<", "<<#var3<<") = ("                                                   \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ")" << std::endl << std::flush << std::noboolalpha; }                   \
} while(0)


#define PrintVar4(cond,var1,var2,var3,var4,os)                                                                                \
do                                                                                                                            \
{ if ((cond))                                                                                                                 \
    { os << std::boolalpha << "("<<#var1<<", "<<#var2<<", "<<#var3<<", "<<#var4<<") = ("                                      \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ", " << (var4) << ")" << std::endl << std::flush << std::noboolalpha; } \
} while(0)


#define PrintVar5(cond,var1,var2,var3,var4,var5,os)                                                      \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << std::boolalpha << "("<<#var1<<", "<<#var2<<", "<<#var3<<", "<<#var4<<", "<<#var5<<") = ("    \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ", " << (var4) << ", " << (var5) <<")"             \
    << std::endl << std::flush << std::noboolalpha; }                                                    \
} while(0)


#define PrintVar6(cond,var1,var2,var3,var4,var5,var6,os)                                                                  \
do                                                                                                                        \
{ if ((cond))                                                                                                             \
    { os << std::boolalpha << "("<<#var1<<", "<<#var2<<", "<<#var3<<", "<<#var4<<", "<<#var5<<", "<<#var6<<") = ("        \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ", " << (var4) << ", " << (var5) <<", " << (var6) <<")"             \
    << std::endl << std::flush << std::noboolalpha; }                                                                     \
} while(0)



#if UTL_VERBOSITY>0

#define utlSAPrintIf(expr)           utlSAGlobalPrintIf(expr)
#define utlSAPrint                   utlSAGlobalPrint
#define utlSAWarning(expr)           utlSAGlobalWarning(expr)
#define utlSAException(expr)         utlSAGlobalException(expr)
#define utlSAAssert(expr)            utlSAGlobalAssert(expr) 

// #define utlSAShowPosition(expr)      if (!(expr)) ; else utlSAPrint

#define utlException(cond,expout)    utlGlobalException(cond,expout)
#define utlAssert(cond,expout)       utlGlobalAssert(cond,expout)
#define utlOSWarning(cond,expout,os) utlOSGlobalWarning(cond,expout,os)
#define utlWarning(cond,expout)      utlGlobalWarning(cond,expout)


#define utlDebug(cond,expout)                                                                                                                          \
do                                                                                                                                                     \
{  if ((cond))                                                                                                                                         \
    { std::cerr << "\n"<<__UTL_DEBUG_STRING<<": " << __utlConditionSucceedPrint(cond) << __UTL_BOLD("msg")<<": '"<<expout << "\n" << std::flush; }     \
} while(0)


#define utlOSShowPosition(cond,os)    utlVLogOSPosition_IF(cond, 0,  os)                                                          

#define utlShowPosition(cond)   utlOSShowPosition(cond,std::cout)                                                                     

#define utlPrintVar1(cond,var)                                              PrintVar1(cond,var,std::cout)
#define utlPrintVar2(cond,var1,var2)                                        PrintVar2(cond,var1,var2,std::cout)
#define utlPrintVar3(cond,var1,var2,var3)                                   PrintVar3(cond,var1,var2,var3,std::cout) 
#define utlPrintVar4(cond,var1,var2,var3,var4)                              PrintVar4(cond,var1,var2,var3,var4,std::cout) 
#define utlPrintVar5(cond,var1,var2,var3,var4,var5)                         PrintVar5(cond,var1,var2,var3,var4,var5,std::cout)
#define utlPrintVar6(cond,var1,var2,var3,var4,var5,var6)                    PrintVar6(cond,var1,var2,var3,var4,var5,var6,std::cout)

#define utlOSPrintVar1(cond,var,os)                                         PrintVar1(cond,var,os)
#define utlOSPrintVar2(cond,var1,var2,os)                                   PrintVar2(cond,var1,var2,os)
#define utlOSPrintVar3(cond,var1,var2,var3,os)                              PrintVar3(cond,var1,var2,var3,os) 
#define utlOSPrintVar4(cond,var1,var2,var3,var4,os)                         PrintVar4(cond,var1,var2,var3,var4,os) 
#define utlOSPrintVar5(cond,var1,var2,var3,var4,var5,os)                    PrintVar5(cond,var1,var2,var3,var4,var5,os)
#define utlOSPrintVar6(cond,var1,var2,var3,var4,var5,var6,os)               PrintVar6(cond,var1,var2,var3,var4,var5,var6,os)

#else

#define utlSAPrintIf(expr)    utlSAGlobalPrintIf(false)
#define utlSAPrint            utlSAPrintIf("")
#define utlSAWarning(expr)    utlSAGlobalWarning(false)
#define utlSAException(expr)  utlSAGlobalException(false)
#define utlSAAssert(expr)     utlSAGlobalAssert(true)

// #define utlSAShowPosition(expr)     if (true) ; else utlSAPrint

#define utlException(cond,expout) 
#define utlDebug(cond,expout) 
#define utlAssert(cond,expout) 
#define utlOSWarning(cond,expout,os) 
#define utlWarning(cond,expout) 
#define utlOSShowPosition(cond,os)
#define utlShowPosition(cond)

#define utlPrintVar1(cond,var)                                                                       
#define utlPrintVar2(cond,var1,var2)                                                                 
#define utlPrintVar3(cond,var1,var2,var3)                                                            
#define utlPrintVar4(cond,var1,var2,var3,var4)                                                       
#define utlPrintVar5(cond,var1,var2,var3,var4, var5)                                                 
#define utlPrintVar6(cond,var1,var2,var3,var4, var5, var6)                                           

#define utlOSPrintVar1(cond,var)                                                                       
#define utlOSPrintVar2(cond,var1,var2)                                                                 
#define utlOSPrintVar3(cond,var1,var2,var3)                                                            
#define utlOSPrintVar4(cond,var1,var2,var3,var4)                                                       
#define utlOSPrintVar5(cond,var1,var2,var3,var4, var5)                                                 
#define utlOSPrintVar6(cond,var1,var2,var3,var4, var5, var6)                                           

#endif

#define utlSetMacro(name, type)                     \
  virtual void Set##name (const type _arg)          \
    {                                               \
    if ( this->m_##name != _arg )                   \
      {                                             \
      this->m_##name = _arg;                        \
      }                                             \
    }


#define utlGetMacro(name, type)                     \
  virtual type Get##name () const                   \
    {                                               \
    return this->m_##name;                          \
    }


#define utlSetGetMacro(name, type) \
  utlSetMacro(name, type);         \
  utlGetMacro(name, type);

    /** @} */

#endif 
