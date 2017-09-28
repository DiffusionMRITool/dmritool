/**
 *       @file  exprtk_lib.h
 *      @brief  
 *     Created  "04-04-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __exprtk_lib_h
#define __exprtk_lib_h

/** http://partow.net/programming/exprtk/index.html  */
#include "exprtk.hpp"


// use extern to build libs, so that it is fast in building
// Use extern template instantiations (since C++11) 
// http://stackoverflow.com/questions/373142/what-techniques-can-be-used-to-speed-up-c-compilation-times 

extern template class exprtk::expression<double>;
//extern template class exprtk::symbol_table<double>;
extern template class exprtk::parser<double>;

#endif 


