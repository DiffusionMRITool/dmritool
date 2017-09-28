/**
 *       @file  test_exprtk.cxx
 *      @brief  
 *     Created  "04-04-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "exprtk_lib.h"
// #include "exprtk.hpp"

template <typename T>
inline T myotherfunc(T v0, T v1, T v2)
{
   return std::abs(v0 - v1) * v2;
}


template <typename T>
void custom_function()
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;

   std::string expression_string =
                  "otherfunc(3 * y, x / 2, abs(-x) * y)";

   T x = T(1);
   T y = T(2);
   // myfunc<T> mf;

   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);
   symbol_table.add_variable("y",y);
   // symbol_table.add_function("myfunc",mf);
   symbol_table.add_function("otherfunc",myotherfunc);
   symbol_table.add_function("std::abs",std::abs);
   symbol_table.add_constants();

   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);

   T result = expression.value();
   printf("Result: %10.5f\n",result);
}

template <typename T>
void vector_function()
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;

   std::string expression_string =
                  "x[0]+x[1]";
                  // " for (var i := 0; i < min(x[],y[],z[]); i += 1) "
                  // " {                                              "
                  // "    z[i] := 3sin(x[i]) + 2log(y[i]);            "
                  // " }                                              ";

   std::vector<T> x({ T(1.1), T(2.2), T(3.3), T(4.4), T(5.5) });
   // T* x = &vec.front();
   // T x[] = { T(1.1), T(2.2), T(3.3), T(4.4), T(5.5) };
   T y[] = { T(2.1), T(2.2), T(3.3), T(4.4), T(5.5) };
   // T z[] = { T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) };

   symbol_table_t symbol_table;
   symbol_table.add_vector("x",x);
   // symbol_table.add_vector("x",x);
   // symbol_table.add_vector("y",y);
   // symbol_table.add_vector("z",z);
   // symbol_table.add_variable("f",f);
   symbol_table.add_constants();

   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);

   T result = expression.value();
   printf("Result: %10.5f\n", result);
}

int 
main (int argc, char const* argv[])
{
   custom_function<double>();

   vector_function<double>();

   typedef exprtk::symbol_table<double> symbol_table_t;
   typedef exprtk::expression<double>     expression_t;
   typedef exprtk::parser<double>             parser_t;

     {
     std::vector<double> x({ 1.1, 2.2, 3.3, 4.4, 5.5 });
     std::vector<double> x1({ 2.1, 2.2, 3.3, 4.4, 5.5 });
     double y[] = {2.1, 2.2, 3.3, 4.4, 5.5};

     symbol_table_t symbol_table;
     symbol_table.add_vector("x",x);

     std::string expression_string =
       "x[0]+x[1]";

     expression_t expression;

     parser_t parser;

     expression.register_symbol_table(symbol_table);

     double result;

     parser.compile(expression_string,expression);
     result = expression.value();
     printf("Result: %10.5f\n", result);

     symbol_table.get_vector("x")->operator=(x1);
     parser.compile(expression_string,expression);
     result = expression.value();
     printf("Result: %10.5f\n", result);

     std::vector<double> x2({1.1, 2.3});
     symbol_table.get_vector("x")->operator=(x2);
     parser.compile(expression_string,expression);
     result = expression.value();
     printf("Result: %10.5f\n", result);
     }
  
  return 0;
}
