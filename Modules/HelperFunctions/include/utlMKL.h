/**
 *       @file  utlMKL.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-18-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlMKL_h
#define __utlMKL_h

#include <iostream>

#include "mkl_service.h"
#include "mkl.h"

namespace utl
{

/** @addtogroup utlHelperFunctions
@{ */

/** in-place copy or transpose  */
template <class T> inline void 
mkl_imatcopy(const char ordering, const char trans, const int rows, const int cols, const T alpha, T* A, const int lda, const int ldb);

/** out-place copy or transpose  */
template <class T> inline void 
mkl_omatcopy(const char ordering, const char trans, const int rows, const int cols, const T alpha, const T* A, const int lda, T* B, const int ldb);

template <> inline void 
mkl_imatcopy<double>(const char ordering, const char trans, const int rows, const int cols, const double alpha, double* A, const int lda, const int ldb)
{
  mkl_dimatcopy(ordering,trans,rows,cols,alpha,A,lda,ldb);
}
template <> inline void 
mkl_imatcopy<float>(const char ordering, const char trans, const int rows, const int cols, const float alpha, float* A, const int lda, const int ldb)
{
  mkl_simatcopy(ordering,trans,rows,cols,alpha,A,lda,ldb);
}

template <> inline void 
mkl_omatcopy<double>(const char ordering, const char trans, const int rows, const int cols, const double alpha, const double* A, const int lda, double* B, const int ldb)
{
  mkl_domatcopy(ordering,trans,rows,cols,alpha,A,lda,B,ldb);
}
template <> inline void 
mkl_omatcopy<float>(const char ordering, const char trans, const int rows, const int cols, const float alpha, const float* A, const int lda, float* B, const int ldb)
{
  mkl_somatcopy(ordering,trans,rows,cols,alpha,A,lda,B,ldb);
}

    /** @} */

}

#endif 
