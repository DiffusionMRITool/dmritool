/**
 *       @file  utlBlas.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __utlBlas_h
#define __utlBlas_h

/** @addtogroup utlMath
@{ */

#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus

#ifdef UTL_USE_MKL
  #include <mkl_cblas.h>
#else
  #include <utl_cblas.h>
#endif

#ifdef __cplusplus
}
#endif //__cplusplus

#include <vector>
#include <cmath>
#include <algorithm>

namespace utl
{

#if defined (MKL_INT)
  typedef MKL_INT INTT;
#elif defined (OPENBLAS_CONFIG_H)
  typedef blasint INTT;
#elif defined (INT_64BITS)
  typedef long long int INTT;
#else
  typedef int INTT;
#endif


/**************************************************************************************/
/** template function definitions. Level 1 */
/**************************************************************************************/

template <class T> inline CBLAS_INDEX 
cblas_iamax(const INTT N, const T *X, const INTT incX);
template <class T> inline T
cblas_nrm2(const INTT N, const T *X, const INTT incX);
template <class T> inline T
cblas_asum(const INTT N, const T *X, const INTT incX);

// #ifdef UTL_USE_MKL
template <class T> inline CBLAS_INDEX 
cblas_iamin(const INTT N, const T *X, const INTT incX);
// #endif

/** dot product between two vectors  */
template <class T> inline T
cblas_dot( const INTT N, const T *X, const INTT incX, const T *Y, const INTT incY);
/** copy from one vector to another vector  */
template <class T> inline void
cblas_copy( const INTT N, const T *X, const INTT incX, T *Y, const INTT incY);
/** swap two vectors  */
template <class T> inline void
cblas_swap( const INTT N, T *X, const INTT incX, T *Y, const INTT incY);
/** Computes the product of a vector by a scalar. x = a*x */
template <class T> inline void
cblas_scal (const INTT N, const T alpha, T *X, const INTT incX);


/**************************************************************************************/
/** template function definitions. Level 2 */
/**************************************************************************************/

/** out-product of two vectors. \f$ A = alpha*x*y'+ A, \f$   */
template <class T> inline void 
cblas_ger( const CBLAS_ORDER order, const INTT M, const INTT N, const T alpha, const T *X, const INTT incX, const T *Y, const INTT incY, T *A, const INTT lda);

/** out-product of one vector. \f$ A = alpha*x*x' + A \f$  */
template <class T> inline void 
cblas_syr (const CBLAS_ORDER order, const CBLAS_UPLO Uplo, const INTT N, const T alpha, const T *X, const INTT incX, T *A, const INTT lda);

/** Matrix-vector product of general matrices. \f$ y := alpha*A*x + beta*y \f$   */
template <class T> inline void 
cblas_gemv( const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const INTT M, const INTT N, const T alpha, const T *A, const INTT lda, const T *X, const INTT incX, const T beta, T *Y, const INTT incY);


/**************************************************************************************/
/** template function definitions. Level 3 */
/**************************************************************************************/

/** Matrix-matrix product of general matrices   */
template <class T> inline void 
cblas_gemm( const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const INTT M, const INTT N, const INTT K, const T alpha, const T *A, const INTT lda, const T *B, const INTT ldb, const T beta, T *C, const INTT ldc);

/** Rank-k update—multiplies a symmetric matrix by its transpose and adds a second matrix.  */
template <class T> inline void
cblas_syrk( CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, INTT N,  INTT K, T alpha, T *A,  INTT lda, T beta, T*C, INTT ldc);



template <class T> inline void 
cblas_axpby (const INTT N, const T alpha, const T *X, const INTT incX, const T beta, T *Y, const INTT incY);

/**************************************************************************************/
/** function implementations using double and float  */

template <> inline CBLAS_INDEX
cblas_iamax<double>(const INTT N, const double *X, const INTT incX)
{
  return cblas_idamax(N,X,incX);
}
template <> inline CBLAS_INDEX
cblas_iamax<float>(const INTT N, const float *X, const INTT incX)
{
  return cblas_isamax(N,X,incX);
}

#ifdef UTL_USE_MKL
template <> inline CBLAS_INDEX
cblas_iamin<double>(const INTT N, const double *X, const INTT incX)
{
  return cblas_idamin(N,X,incX);
}
template <> inline CBLAS_INDEX
cblas_iamin<float>(const INTT N, const float *X, const INTT incX)
{
  return cblas_isamin(N,X,incX);
}
#else
template <typename T> inline CBLAS_INDEX
cblas_iamin(const INTT N, const T *X, const INTT incX)
{
  int iter;
  std::vector<T> vec(N);
  for ( int i = 0; i < N; i = i+incX ) 
    vec[i]=std::fabs(X[i]);
  iter = std::min_element(vec.begin(), vec.end());
  int n=0;
  for ( int ii = vec.begin();  ii!=iter; ++ii ) 
    n++;
  return n;
}
#endif

template <> inline double
cblas_nrm2<double>( const INTT N, const double *X, const INTT incX)
{
  return cblas_dnrm2(N,X,incX);
}
template <> inline float
cblas_nrm2<float>( const INTT N, const float *X, const INTT incX)
{
  return cblas_snrm2(N,X,incX);
}

template <> inline double
cblas_asum<double>( const INTT N, const double *X, const INTT incX)
{
  return cblas_dasum(N,X,incX);
}
template <> inline float
cblas_asum<float>( const INTT N, const float *X, const INTT incX)
{
  return cblas_sasum(N,X,incX);
}


template <> inline double
cblas_dot<double>( const INTT N, const double *X, const INTT incX, const double *Y, const INTT incY)
{
  return cblas_ddot(N,X,incX,Y,incY);
}
template <> inline float
cblas_dot<float>( const INTT N, const float *X, const INTT incX, const float *Y, const INTT incY)
{
  return cblas_sdot(N,X,incX,Y,incY);
}

template <> inline void
cblas_copy<double>( const INTT N, const double *X, const INTT incX, double *Y, const INTT incY)
{
  cblas_dcopy(N,X,incX,Y,incY);
}
template <> inline void
cblas_copy<float>( const INTT N, const float *X, const INTT incX, float *Y, const INTT incY)
{
  cblas_scopy(N,X,incX,Y,incY);
}

template <> inline void
cblas_swap<double>( const INTT N, double *X, const INTT incX, double *Y, const INTT incY)
{
  cblas_dswap(N,X,incX,Y,incY);
}
template <> inline void
cblas_swap<float>( const INTT N, float *X, const INTT incX, float *Y, const INTT incY)
{
  cblas_sswap(N,X,incX,Y,incY);
}

template <> inline void
cblas_scal<double>(const INTT N, const double alpha, double *X, const INTT incX)
{
  cblas_dscal(N, alpha, X, incX);
}
template <> inline void
cblas_scal<float>(const INTT N, const float alpha, float *X, const INTT incX)
{
  cblas_sscal(N, alpha, X, incX);
}




template <> inline void 
cblas_ger<float>( const CBLAS_ORDER order, const INTT M, const INTT N, const float alpha, const float *X, const INTT incX, const float *Y, const INTT incY, float *A, const INTT lda)
{
  cblas_sger(order,M,N,alpha,X,incX,Y,incY,A,lda);
}
template <> inline void 
cblas_ger<double>( const CBLAS_ORDER order, const INTT M, const INTT N, const double alpha, const double *X, const INTT incX, const double *Y, const INTT incY, double *A, const INTT lda)
{
  cblas_dger(order,M,N,alpha,X,incX,Y,incY,A,lda);
}

template <> inline void 
cblas_syr<float>(const CBLAS_ORDER order, const CBLAS_UPLO Uplo, const INTT N, const float alpha, const float *X, const INTT incX, float *A, const INTT lda)
{
  cblas_ssyr(order,Uplo,N,alpha,X,incX,A,lda);
}
template <> inline void 
cblas_syr<double>(const CBLAS_ORDER order, const CBLAS_UPLO Uplo, const INTT N, const double alpha, const double *X, const INTT incX, double *A, const INTT lda)
{
  cblas_dsyr(order,Uplo,N,alpha,X,incX,A,lda);
}

template <> inline void 
cblas_gemv<double>( const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const INTT M,  const INTT N,  const double alpha, const double *A,  const INTT lda, const  double *X,  const INTT incX, const double beta, double *Y,  const INTT incY)
{
   cblas_dgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}
template <> inline void 
cblas_gemv<float>( const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const INTT M, const INTT N, const float alpha, const float *A,  const INTT lda, const  float *X,  const INTT incX, const float beta, float *Y,  const INTT incY)
{
   cblas_sgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}





template <> inline void 
cblas_gemm<double>( const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const INTT M,  const INTT N,  const INTT K, const double alpha, const double *A,  const INTT lda,  const double *B, const  INTT ldb, const double beta, double *C, const INTT ldc) 
{
   cblas_dgemm(order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
}
template <> inline void 
cblas_gemm<float>( const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,  const CBLAS_TRANSPOSE TransB, const INTT M,  const INTT N,  const INTT K, const float alpha, const float *A,  const INTT lda,  const float *B, const  INTT ldb, const float beta, float *C,  const INTT ldc) 
{
   cblas_sgemm(order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
}

template <> inline void
cblas_syrk<double>( CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, INTT N,  INTT K, double alpha, double *A,  INTT lda, double beta, double *C, INTT ldc)
{
  cblas_dsyrk(order,Uplo,Trans,N,K,alpha,A,lda,beta,C,ldc);
}
template <> inline void
cblas_syrk<float>( CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, INTT N,  INTT K, float alpha, float *A,  INTT lda, float beta, float *C, INTT ldc)
{
  cblas_ssyrk(order,Uplo,Trans,N,K,alpha,A,lda,beta,C,ldc);
}


#ifdef UTL_USE_MKL
template <> inline void 
cblas_axpby<double>(const INTT N, const double alpha, const double *X, const INTT incX, const double beta, double *Y, const INTT incY)
{
  cblas_daxpby(N,alpha,X,incX,beta,Y,incY);
}
template <> inline void 
cblas_axpby<float>(const INTT N, const float alpha, const float *X, const INTT incX, const float beta, float *Y, const INTT incY)
{
  cblas_saxpby(N,alpha,X,incX,beta,Y,incY);
}
#else
template <typename T> inline void 
cblas_axpby(const INTT N, const T alpha, const T *X, const INTT incX, const T beta, T *Y, const INTT incY)
{
  for ( int i=0, j=0; i<N, j<N; i+=incX, j+=incY ) 
    Y[j] = alpha*X[i] + beta*Y[j];
}

#endif

}


/** 
*  macro to define gemm for row-major matrix product
*
* http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ C = alpha * A * B + beta * C \f$
*  \note: C should be pre-allocated
* */
#define __utl_gemm_MatrixTimesMatrix(T, FuncName, FuncHelperName, RowMajorMatrixName, GetRowsFuncName, GetColsFuncName, GetDataFuncName, ReSizeFuncName) \
template <class T> \
inline bool \
FuncName(const bool bATrans, const bool bBTrans, const T alpha, const RowMajorMatrixName& A, const RowMajorMatrixName& B, const T beta, RowMajorMatrixName& C)  \
{ \
  CBLAS_TRANSPOSE transA, transB; \
  unsigned int M = 0, N = 0, K = 0; \
  if(!bATrans && !bBTrans) \
    { \
    transA = CblasNoTrans; \
    transB = CblasNoTrans; \
    M = A.GetRowsFuncName(); \
    N = B.GetColsFuncName(); \
    K = A.GetColsFuncName(); \
    utlSAException(A.GetColsFuncName() != B.GetRowsFuncName())(A.GetColsFuncName())(B.GetRowsFuncName()).msg("matrix dimension mismatch"); \
    } \
  else if(bATrans && !bBTrans)  \
    { \
    transA = CblasTrans; \
    transB = CblasNoTrans; \
    M = A.GetColsFuncName(); \
    N = B.GetColsFuncName(); \
    K = A.GetRowsFuncName(); \
    utlSAException(A.GetRowsFuncName() != B.GetRowsFuncName())(A.GetRowsFuncName())(B.GetRowsFuncName()).msg("matrix dimension mismatch"); \
    } \
  else if(!bATrans && bBTrans)  \
    { \
    transA = CblasNoTrans; \
    transB = CblasTrans; \
    M = A.GetRowsFuncName(); \
    N = B.GetRowsFuncName(); \
    K = A.GetColsFuncName(); \
    utlSAException(A.GetColsFuncName() != B.GetColsFuncName())(A.GetColsFuncName())(B.GetColsFuncName()).msg("matrix dimension mismatch"); \
    } \
  else \
    { \
    transA = CblasTrans; \
    transB = CblasTrans; \
    M = A.GetColsFuncName(); \
    N = B.GetRowsFuncName(); \
    K = A.GetRowsFuncName(); \
    utlSAException(A.GetRowsFuncName() != B.GetColsFuncName())(A.GetRowsFuncName())(B.GetColsFuncName()).msg("matrix dimension mismatch"); \
    } \
  utl::cblas_gemm(CblasRowMajor, transA, transB, M, N, K, alpha, A.GetDataFuncName(), A.GetColsFuncName(), B.GetDataFuncName(), B.GetColsFuncName(), beta, C.GetDataFuncName(), C.GetColsFuncName()); \
  return true; \
} \
\
template <class T>  \
inline void \
Product##FuncHelperName##MM(const RowMajorMatrixName& A, const RowMajorMatrixName& B, RowMajorMatrixName& C, const double alpha=1.0, const double beta=0.0)  \
{ \
  C.ReSizeFuncName(A.GetRowsFuncName(), B.GetColsFuncName()); \
  FuncName<T>(false, false, alpha, A, B, beta, C); \
} \
 \
template <class T>  \
inline void \
Product##FuncHelperName##MM(const RowMajorMatrixName& A1, const RowMajorMatrixName& A2, const RowMajorMatrixName& A3, RowMajorMatrixName& C)  \
{ \
  RowMajorMatrixName tmp; \
  Product##FuncHelperName##MM<T>(A1, A2, tmp); \
  Product##FuncHelperName##MM<T>(tmp, A3, C); \
} \
 \
template <class T>  \
inline void \
Product##FuncHelperName##MM(const RowMajorMatrixName& A1, const RowMajorMatrixName& A2, const RowMajorMatrixName& A3, const RowMajorMatrixName& A4, RowMajorMatrixName& C)  \
{ \
  RowMajorMatrixName tmp; \
  Product##FuncHelperName##MM<T>(A1, A2, A3, tmp); \
  Product##FuncHelperName##MM<T>(tmp, A4, C); \
} \
 \
template <class T>  \
inline void \
Product##FuncHelperName##MM(const RowMajorMatrixName& A1, const RowMajorMatrixName& A2, const RowMajorMatrixName& A3, const RowMajorMatrixName& A4, const RowMajorMatrixName& A5, RowMajorMatrixName& C)  \
{ \
  RowMajorMatrixName tmp; \
  Product##FuncHelperName##MM<T>(A1, A2, A3, A4, tmp); \
  Product##FuncHelperName##MM<T>(tmp, A5, C); \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##MtM(const RowMajorMatrixName& A, const RowMajorMatrixName& B, RowMajorMatrixName& C, const double alpha=1.0, const double beta=0.0)  \
{ \
  C.ReSizeFuncName(A.GetColsFuncName(), B.GetColsFuncName()); \
  FuncName<T>(true, false, alpha, A, B, beta, C); \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##MMt(const RowMajorMatrixName& A, const RowMajorMatrixName& B, RowMajorMatrixName& C, const double alpha=1.0, const double beta=0.0)  \
{ \
  C.ReSizeFuncName(A.GetRowsFuncName(), B.GetRowsFuncName()); \
  FuncName<T>(false, true, alpha, A, B, beta, C); \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##MtMt(const RowMajorMatrixName& A, const RowMajorMatrixName& B, RowMajorMatrixName& C, const double alpha=1.0, const double beta=0.0)  \
{ \
  C.ReSizeFuncName(A.GetColsFuncName(), B.GetRowsFuncName()); \
  FuncName<T>(true, true, alpha, A, B, beta, C); \
} \


/** 
*  macro to define gemm for row-major matrix and vector product
*
*  http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ Y = alpha * A * X + beta * Y \f$
*  \note: Y should be pre-allocated
* */

#define __utl_gemv_MatrixTimesVector(T, FuncName, FuncHelperName, RowMajorMatrixName, GetRowsFuncName, GetColsFuncName, MatrixGetDataFuncName, VectorName, GetSizeFuncName, VectorGetDataFuncName, ReSizeFuncName) \
template <class T>  \
inline bool  \
FuncName(const bool bATrans, const T alpha, const RowMajorMatrixName& A, const VectorName& X, const T beta, VectorName& Y)  \
{ \
  CBLAS_TRANSPOSE TransA; \
  if(bATrans)  \
    { \
    TransA = CblasTrans; \
    utlSAException(A.GetRowsFuncName() != X.GetSizeFuncName())(A.GetRowsFuncName())(X.GetSizeFuncName()).msg("matrix and vector dimension mismatch"); \
    } \
  else  \
    { \
    TransA = CblasNoTrans; \
    utlSAException(A.GetColsFuncName() != X.GetSizeFuncName())(A.GetColsFuncName())(X.GetSizeFuncName()).msg("matrix and vector dimension mismatch"); \
    } \
  utl::cblas_gemv<T>(CblasRowMajor, TransA, A.GetRowsFuncName(), A.GetColsFuncName(), alpha, A.MatrixGetDataFuncName(), A.GetColsFuncName(), X.VectorGetDataFuncName(), 1, beta, Y.VectorGetDataFuncName(), 1); \
  return true; \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##Mv(const RowMajorMatrixName& A, const VectorName& b, VectorName& c, const double alpha=1.0, const double beta=0.0)  \
{ \
  c.ReSizeFuncName(A.GetRowsFuncName()); \
  FuncName<T>(false, alpha, A, b, beta, c); \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##Mtv(const RowMajorMatrixName& A, const VectorName& b, VectorName& c, const double alpha=1.0, const double beta=0.0)  \
{ \
  c.ReSizeFuncName(A.GetColsFuncName()); \
  FuncName<T>(true, alpha, A, b, beta, c); \
} \
 \


/** 
*  \f$ Y = alpha X^T * A + beta * Y \f$ 
*  \note: Y should be pre-allocated 
* */ 
#define __utl_gevm_MatrixTimesVector(T, FuncName, FuncHelperName, RowMajorMatrixName, GetRowsFuncName, GetColsFuncName, MatrixGetDataFuncName, VectorName, GetSizeFuncName, VectorGetDataFuncName, ReSizeFuncName) \
template <class T>  \
inline bool  \
FuncName(const bool bATrans, const T alpha, const VectorName& X, const RowMajorMatrixName& A, const T beta, VectorName& Y)  \
{ \
  CBLAS_TRANSPOSE transA; \
  unsigned int M = 1, N = 0, K = 0; \
  if(bATrans)  \
    { \
    transA = CblasTrans; \
    N = A.GetRowsFuncName(); \
    K = A.GetColsFuncName(); \
    utlSAException(A.GetColsFuncName() != X.GetSizeFuncName())(A.GetColsFuncName())(X.GetSizeFuncName()).msg("matrix and vector dimension mismatch"); \
    } \
  else  \
    { \
    N = A.GetColsFuncName(); \
    K = A.GetRowsFuncName(); \
    transA = CblasNoTrans; \
    utlSAException(A.GetRowsFuncName() != X.GetSizeFuncName())(A.GetRowsFuncName())(X.GetSizeFuncName()).msg("matrix and vector dimension mismatch"); \
    } \
  utl::cblas_gemm<T>(CblasRowMajor, CblasNoTrans, transA, M, N, K, alpha, X.VectorGetDataFuncName(), X.GetSizeFuncName(), A.MatrixGetDataFuncName(), A.GetColsFuncName(), beta, Y.VectorGetDataFuncName(), Y.GetSizeFuncName()); \
  return true; \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##vM(const VectorName& b, const RowMajorMatrixName& A, VectorName& c, const double alpha=1.0, const double beta=0.0)  \
{ \
  c.ReSizeFuncName(A.GetColsFuncName()); \
  FuncName<T>(false, alpha, b, A, beta, c); \
} \
 \
template <class T>  \
inline void  \
Product##FuncHelperName##vMt(const VectorName& b, const RowMajorMatrixName& A, VectorName& c, const double alpha=1.0, const double beta=0.0)  \
{ \
  c.ReSizeFuncName(A.GetRowsFuncName()); \
  FuncName<T>(true, alpha, b, A, beta, c); \
} 


/**
 * \brief macro for syrk and row-major matrix
 *
 * \param trans If false, then C := alpha* A*A' + beta* C; If true, then C := alpha* A'*A + beta* C
 * \param alpha 
 * \param A  MxN matrix
 * \param beta
 * \param C MxM or NxN symmetric matrix
 */

#define __utl_syrk_Matrix(T, FuncName, FuncHelperName, RowMajorMatrixName, GetRowsFuncName, GetColsFuncName, GetDataFuncName, ReSizeFuncName) \
template <class T> \
inline void  \
FuncName( const bool trans, const T alpha, const RowMajorMatrixName& A, const T beta, RowMajorMatrixName& C ) \
{ \
  int size = C.GetRowsFuncName()*C.GetColsFuncName(); \
  utlSAException(size>0 && C.GetRowsFuncName()!=C.GetColsFuncName())(C.GetRowsFuncName())(C.GetColsFuncName()).msg("C should be symmetric matrix or empty matrix"); \
  CBLAS_TRANSPOSE transBlas; \
  int K, N; \
  if (!trans) \
    { \
    utlSAException(size>0 && A.GetRowsFuncName()!=C.GetRowsFuncName())(A.GetRowsFuncName())(C.GetRowsFuncName()).msg("wrong size"); \
    transBlas=CblasNoTrans; \
    N = A.GetRowsFuncName(); \
    K = A.GetColsFuncName(); \
    } \
  else  \
    { \
    utlSAException(size>0 && A.GetColsFuncName()!=C.GetRowsFuncName())(A.GetColsFuncName())(C.GetRowsFuncName()).msg("wrong size"); \
    transBlas=CblasTrans; \
    N = A.GetColsFuncName(); \
    K = A.GetRowsFuncName(); \
    } \
  utlException(size>0 && (C.GetRowsFuncName()!=N || C.GetColsFuncName()!=N), "wrong size of C"); \
  if (size==0) \
    { \
    utlSAException(std::fabs(beta)>1e-10)(beta)(C.GetRowsFuncName())(C.GetColsFuncName()).msg("when C is empty matrix, beta should be zero"); \
    C.ReSizeFuncName(N,N); \
    } \
  RowMajorMatrixName B=A; \
  utl::cblas_syrk<T>(CblasRowMajor, CblasUpper, transBlas, N, K, alpha, B.GetDataFuncName(), B.GetColsFuncName(), beta, C.GetDataFuncName(), N); \
 \
  T* data = C.GetDataFuncName(); \
  for ( int i = 0; i < C.GetRowsFuncName(); ++i )  \
    for ( int j = 0; j < i; ++j )  \
      data[i*N+j] = data[j*N+i]; \
} \
 \
template <class T> \
void  \
Product##FuncHelperName##XXt ( const RowMajorMatrixName& A, RowMajorMatrixName& C, const double alpha=1.0, const double beta=0.0 ) \
{ \
  FuncName<T>(false, alpha, A, beta, C); \
} \
template <class T> \
void  \
Product##FuncHelperName##XtX ( const RowMajorMatrixName& A, RowMajorMatrixName& C, const double alpha=1.0, const double beta=0.0 ) \
{ \
  FuncName<T>(true, alpha, A, beta, C); \
}

    /** @} */

#endif 

