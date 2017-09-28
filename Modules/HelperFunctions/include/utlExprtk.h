/**
 *       @file  utlExprtk.h
 *      @brief  
 *     Created  "04-05-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlExprtk_h
#define __utlExprtk_h

#include "exprtk_lib.h"
#include <functional>
#include <vector>

namespace utl
{

/** 
* convert a string to std::function which maps T to T.  
*  Example: "cos(x)+ 2*x+ x^3",  "if (x>3) 2; else 1;"  
*  
*  NOTE: it is slow if it is called for many times. It is better to put symbol_table, expression, and parser definition outside of the loop.
*  */
template <class T=double>
inline std::function<T(T)> 
GetScalarFunctionFromString ( const std::string funcStr )
{
  auto func = 
    [funcStr](T xval)
      {
      T x(xval);
      exprtk::symbol_table<T> symbol_table;
      symbol_table.add_variable("x",x);
      symbol_table.add_constants();

      exprtk::expression<T> expression;
      expression.register_symbol_table(symbol_table);

      exprtk::parser<T> parser;
      parser.compile(funcStr ,expression);

      T result = expression.value();
      return result;
      }; 
  return func;
}


/** 
* convert a string to std::function which maps std::vector<T> to T.  
*  Example: "cos(x[0])+ 2*x[1]",  "max(x[0], x[1])",  "if (x[0]>1.0) 1.0; else x[1]"  "x[0]>x[1]"
*
*  NOTE: it is slow if it is called for many times. It is better to put symbol_table, expression, and parser definition outside of the loop.
*  */
template <class T=double>
inline std::function<T(std::vector<T>)> 
GetVectorFunctionFromString ( const std::string funcStr )
{
  auto func = 
    [funcStr](const std::vector<T>& xval)
      {
      std::vector<T> x(xval);
      exprtk::symbol_table<T> symbol_table;
      symbol_table.add_vector("x",x);
      symbol_table.add_constants();

      exprtk::expression<T> expression;
      expression.register_symbol_table(symbol_table);

      exprtk::parser<T> parser;
      parser.compile(funcStr,expression);

      T result = expression.value();
      return result;
      }; 
  return func;
}

}

#endif 
