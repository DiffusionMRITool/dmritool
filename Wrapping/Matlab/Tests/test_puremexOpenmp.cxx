/**
 *       @file  test_puremexOpenmp.cxx
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

#include "mex.h"
#include <omp.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  // omp_set_dynamic(1);
  // omp_set_num_threads(3);
#pragma omp parallel
    {
    mexPrintf("Max num threads %d.\n", omp_get_max_threads());
#pragma omp for
    for (int i = 0; i < 10; i++)
      {
      mexPrintf("Num threads %d, thread ID %d.\n", omp_get_num_threads(), omp_get_thread_num());
      }
    }
}

