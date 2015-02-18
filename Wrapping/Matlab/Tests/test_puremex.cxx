/**
 *       @file  test_puremex.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-09-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include <iostream>
#include "mex.h" 
// #include "mexutils.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  // utlGlobalException(nrhs!=1, "Bad number of inputs arguments");
  // double num = mxGetScalar(prhs[0]);
  std::cout << "num = " << 15 << std::endl << std::flush;
}
