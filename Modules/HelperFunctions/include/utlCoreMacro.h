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
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_E
#define M_E             2.7182818284590452354
#endif
#ifndef M_SQRTPI
#define M_SQRTPI   1.7724538509055160272981674833411452 
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


#include "utlSmartAssert.h"


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


// #define __UTL_ERROR_STRING  utl::GetColoredString("Error", COLOR_RED)
// #define __UTL_WARNING_STRING  utl::GetColoredString("Warning", COLOR_RED)
// #define __UTL_DEBUG_STRING  utl::GetColoredString("Debug", COLOR_RED)
// #define __UTL_BOLD(str) utl::GetColoredString(str, COLOR_BOLD)
// #define __UTL_EXPSTR(str) utl::GetColoredString(str, COLOR_GREEN)

#define __UTL_ERROR_STRING  "Error"
#define __UTL_WARNING_STRING  "Warning"
#define __UTL_DEBUG_STRING  "Debug"
#define __UTL_BOLD(str) str
#define __UTL_EXPSTR(str) str

/**
 * \brief  macros which are used for debug 
 * \author Jian Cheng
 */
#define __utlConditionFailPrint(cond)      "In File: " <<__FILE__<< ", Line: " << __LINE__   << ", Function: " << __utl_LOCATION__ <<  "\nExpression: '" << #cond << "' failed. " << "\n"
#define __utlConditionSucceedPrint(cond)   "In File: " <<__FILE__<< ", Line: " << __LINE__   << ", Function: " << __utl_LOCATION__ <<  "\nExpression: '" << #cond << "' satisfied. " << "\n"

#define utlSASetLog(log) Assert::set_log(log) 

#define utlSAGlobalPrintIf(expr)     \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert(#expr).log().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A \

#define utlSAGlobalPrint utlSAGlobalPrintIf("") 

#define utlSAGlobalWarning(expr)        \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert(#expr).warn().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A \

#define utlSAGlobalException(expr)        \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert(#expr).error().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A \

#define utlSAGlobalAssert(expr)        \
    if ( (expr) ) ; \
    else ::smart_assert::make_assert(#expr).error().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_assert).SMART_ASSERT_A \

#define utlGlobalException(cond,expout)                                                                                      \
do                                                                                                                           \
{                                                                                                                            \
  if ((cond))                                                                                                                \
    { std::cerr << "\n"<<__UTL_ERROR_STRING<<": " << __utlConditionSucceedPrint(cond) << expout << "\n" << std::flush;       \
    utlAbort(""); }                                                                                                          \
} while (0)                                                                                               

#define utlGlobalAssert(cond,expout)                                                                                         \
do                                                                                                                           \
{                                                                                                                            \
  if (!(cond))                                                                                                               \
    { std::cerr << "\n"<<__UTL_ERROR_STRING<<": " << __utlConditionFailPrint(cond) << expout << "\n" << std::flush;          \
    utlAbort(""); }                                                                                                          \
} while (0)

#define utlOSGlobalWarning(cond,expout,os)                                                                                   \
do                                                                                                                           \
{  if ((cond))                                                                                                               \
    { os << "\n"<<__UTL_WARNING_STRING<<": "<< __utlConditionSucceedPrint(cond) <<  expout << "\n" << std::flush;   }        \
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
    { os << #var <<" = " <<  (var) << std::endl << std::flush; }                                         \
} while(0)

#define PrintVar2(cond,var1,var2,os)                                                                     \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << "("<<#var1<<", "<< #var2<<") = (" <<  (var1) << ", " << (var2) << ")"                        \
    << std::endl << std::flush; }                                                                        \
} while(0)

#define PrintVar3(cond,var1,var2,var3,os)                                                                \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << "("<<#var1<<", "<<#var2<<", "<<#var3<<") = ("                                                \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ")" << std::endl << std::flush; }                  \
} while(0)

#define PrintVar4(cond,var1,var2,var3,var4,os)                                                           \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << "("<<#var1<<", "<<#var2<<", "<<#var3<<", "<<#var4<<") = ("                                   \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ", " << (var4) << ")" << std::endl << std::flush; }\
} while(0)

#define PrintVar5(cond,var1,var2,var3,var4,var5,os)                                                      \
do                                                                                                       \
{ if ((cond))                                                                                            \
    { os << "("<<#var1<<", "<<#var2<<", "<<#var3<<", "<<#var4<<", "<<#var5<<") = ("                      \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ", " << (var4) << ", " << (var5) <<")"             \
    << std::endl << std::flush; }                                                                        \
} while(0)

#define PrintVar6(cond,var1,var2,var3,var4,var5,var6,os)                                                         \
do                                                                                                               \
{ if ((cond))                                                                                                    \
    { os << "("<<#var1<<", "<<#var2<<", "<<#var3<<", "<<#var4<<", "<<#var5<<", "<<#var6<<") = ("                 \
    <<  (var1) << ", " << (var2) << ", " << (var3) << ", " << (var4) << ", " << (var5) <<", " << (var6) <<")"    \
    << std::endl << std::flush; }                                                                                \
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

#define utlDebug(cond,expout)                                                                                                \
do                                                                                                                           \
{  if ((cond))                                                                                                               \
    { std::cerr << "\n"<<__UTL_DEBUG_STRING<<": " << __utlConditionSucceedPrint(cond) << expout << "\n" << std::flush; }     \
} while(0)

#define utlOSShowPosition(cond,os)                                                                                             \
do                                                                                                                             \
{  if ((cond))                                                                                                                 \
    { os << "\n"<<"Work Flow"<<": In File: " <<__FILE__<< ", Line: " << __LINE__ << "\n"                                       \
    << "Function: " << __utl_LOCATION__ << "\n" << std::flush;                                                                 \
    }                                                                                                                          \
} while(0)

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


    /** @} */

#endif 
