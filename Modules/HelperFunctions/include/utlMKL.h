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
#include <complex>

#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl_service.h"
#include "mkl.h"

namespace utl
{

/** @addtogroup utlMath
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
mkl_imatcopy<std::complex<double> >(const char ordering, const char trans, const int rows, const int cols, const std::complex<double> alpha, std::complex<double>* A, const int lda, const int ldb)
{
  mkl_zimatcopy(ordering,trans,rows,cols,alpha,A,lda,ldb);
}
template <> inline void 
mkl_imatcopy<std::complex<float> >(const char ordering, const char trans, const int rows, const int cols, const std::complex<float> alpha, std::complex<float>* A, const int lda, const int ldb)
{
  mkl_cimatcopy(ordering,trans,rows,cols,alpha,A,lda,ldb);
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
template <> inline void 
mkl_omatcopy<std::complex<double> >(const char ordering, const char trans, const int rows, const int cols, const std::complex<double> alpha, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb)
{
  mkl_zomatcopy(ordering,trans,rows,cols,alpha,A,lda,B,ldb);
}
template <> inline void 
mkl_omatcopy<std::complex<float> >(const char ordering, const char trans, const int rows, const int cols, const std::complex<float> alpha, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb)
{
  mkl_comatcopy(ordering,trans,rows,cols,alpha,A,lda,B,ldb);
}

    /** @} */

}

#endif 
