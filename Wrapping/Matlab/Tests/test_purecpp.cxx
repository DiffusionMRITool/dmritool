/**
 *       @file  test_purecpp.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-07-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "utl.h"

/**
 * \brief  test pure cpp code
 */
int 
main (int argc, char const* argv[])
{
  std::vector<std::vector<std::string> > strMatrix; 
  utl::ReadLines("aa.txt", strMatrix);
  utlGlobalException(true, "adad");
  return 0;
}
