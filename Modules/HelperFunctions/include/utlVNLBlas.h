/**
 *       @file  utlVNLBlas.h
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

#ifndef __utlVNLBlas_h
#define __utlVNLBlas_h

#include "utlCore.h"
#include "utlVNL.h"
#include "utlBlas.h"

#ifdef UTL_USE_MKL
#include "utlMKL.h"
#endif

namespace utl 
{

/** @addtogroup utlMath
@{ */

/**
 * \brief MatrixCopy. A := alpha * op(A)  
 *
 * \tparam T
 * \param mat input matrix 
 * \param matOut output matrix
 * \param alpha scale factor
 * \param trans 'N' or 'n': op(A) is A;  'T' or 't': op(A) is transpose of A; 'C' or 'c': conjugate transpose;  'R' or 'r': conjugate
 */
template <class T> inline void 
MatrixCopy(const vnl_matrix<T>& mat, vnl_matrix<T>& matOut, const T alpha, const char trans='N')
{
#ifdef UTL_USE_MKL
  if (trans=='N' || trans=='n' || trans=='R' || trans=='r')
    matOut.set_size(mat.rows(), mat.cols());
  else if (trans=='T' || trans=='t' || trans=='C' || trans=='c')
    matOut.set_size(mat.cols(), mat.rows());
  else
    utlException(true, "wrong trans");
  if ((trans=='N' || trans=='n') && std::fabs(alpha-1.0)<1e-10 )
    utl::cblas_copy(mat.rows()*mat.cols(), mat.data_block(), 1, matOut.data_block(), 1);
  else
    utl::mkl_omatcopy<T>('R', trans, mat.rows(), mat.cols(), alpha, mat.data_block(), mat.cols(), matOut.data_block(), matOut.cols()); 
#else
  if (trans=='N' || trans=='n')
    matOut = mat;
  else if (trans=='T' || trans=='t')
    matOut = mat.transpose();
  else if (trans=='C' || trans=='c')
    matOut = mat.conjugate_transpose();
  else if (trans=='R' || trans=='r')
    {
    matOut = mat;
    vnl_c_vector<T>::conjugate(matOut.begin(),matOut.begin(),matOut.size());  
    }

  if (std::fabs(alpha-1)>1e-10)
    matOut *= alpha;
#endif
}


/** http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ C = alpha * A * B + beta * C \f$
*  \note: C should be pre-allocated
*
*  define several functions. 
* 
*  template \<class T\> inline bool 
*  gemm_VnlMatrixTimesMatrix(const bool bATrans, const bool bBTrans, const T alpha, const vnl_matrix<T>& a, const vnl_matrix<T>& b, const T beta, vnl_matrix<T>& c);
*
*  \f$ C = alpha * A * B + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMM(const vnl_matrix<T>& A, const vnl_matrix<T>& B, vnl_matrix<T>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A * B^T + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMMt(const vnl_matrix<T>& A, const vnl_matrix<T>& B, vnl_matrix<T>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A^T * B + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMtM(const vnl_matrix<T>& A, const vnl_matrix<T>& B, vnl_matrix<T>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A^T * B^T + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMtMt(const vnl_matrix<T>& A, const vnl_matrix<T>& B, vnl_matrix<T>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A1*A2*A3 \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMM(const vnl_matrix<T>& A1, const vnl_matrix<T>& A2, const vnl_matrix<T>& A3, vnl_matrix<T>& C);
*
*  \f$ C = alpha * A1*A2*A3*A4 \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMM(const vnl_matrix<T>& A1, const vnl_matrix<T>& A2, const vnl_matrix<T>& A3, const vnl_matrix<T>& A4, vnl_matrix<T>& C);
*
*  \f$ C = alpha * A1*A2*A3*A5 \f$ \n
*  template \<class T\> inline void 
*  ProductVnlMM(const vnl_matrix<T>& A1, const vnl_matrix<T>& A2, const vnl_matrix<T>& A3, const vnl_matrix<T>& A4, const vnl_matrix<T>& A5, vnl_matrix<T>& C);
*
* */
__utl_gemm_MatrixTimesMatrix(T, gemm_VnlMatrixTimesMatrix, Vnl, vnl_matrix<T>, rows, cols, data_block, set_size);


/** http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ Y = alpha * A * X + beta * Y \f$
*  \note: Y should be pre-allocated
*
*  define several functions. 
*
*  template \<class T\>  inline bool  
*  gemv_VnlMatrixTimesVector(const bool bATrans, const T alpha, const vnl_matrix<T>& A, const vnl_vector<T>& X, const T beta, vnl_vector<T>& Y);
*
*  \f$ Y = alpha * A * X + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductVnlMv(const vnl_matrix<T>& A, const vnl_vector<T>& b, vnl_vector<T>& c, const double alpha=1.0, const double beta=0.0);
*
*  \f$ Y = alpha * A^T * X + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductVnlMtv(const vnl_matrix<T>& A, const vnl_vector<T>& b, vnl_vector<T>& c, const double alpha=1.0, const double beta=0.0);
*
* */
__utl_gemv_MatrixTimesVector(T, gemv_VnlMatrixTimesVector, Vnl, vnl_matrix<T>, rows, cols, data_block, vnl_vector<T>, size, data_block, set_size);


/**
*  \f$ Y = alpha X^T * A + beta * Y \f$
*  \note: Y should be pre-allocated
*
*  define several functions. 
*
*  template \<class T\>  inline bool  
*  gemm_VnlVectorTimesMatrix(const bool bATrans, const T alpha, const vnl_vector<T>& X, const vnl_matrix<T>& A, const T beta, vnl_vector<T>& Y)  
*
*  \f$ Y = alpha * b * A + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductVnlvM(const vnl_vector<T>& b, const vnl_matrix<T>& A, vnl_vector<T>& c, const double alpha=1.0, const double beta=0.0)  
*
*  \f$ Y = alpha * b * A^T + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductVnlvMt(const vnl_vector<T>& b, const vnl_matrix<T>& A, vnl_vector<T>& c, const double alpha=1.0, const double beta=0.0)  
*
* */
__utl_gevm_MatrixTimesVector(T, gemm_VnlVectorTimesMatrix, Vnl, vnl_matrix<T>, rows, cols, data_block, vnl_vector<T>, size, data_block, set_size);

/**
 * \brief syrk_VnlMatrix 
 *
 *  define several functions. 
 *
 * \param trans If false, then C := alpha* A*A' + beta* C; If true, then C := alpha* A'*A + beta* C
 * \param alpha 
 * \param A  MxN matrix
 * \param beta
 * \param C MxM or NxN symmetric matrix
 *  
 *  template \<class T\> inline void  
 *  syrk_VnlMatrix( const bool trans, const T alpha, const vnl_matrix<T>& A, const T beta, vnl_matrix<T>& C ) 
 *
 *  \f$ C = alpha *A*A^T  + beta*C \f$ \n
 *  template \<class T\> void  
 *  ProductVnlXXt ( const vnl_matrix<T>& A, vnl_matrix<T>& C, const double alpha=1.0, const double beta=0.0 ) 
 *
 *  \f$ C = alpha *A^T*A  + beta*C \f$ \n
 *  template \<class T\> void  
 *  ProductVnlXtX ( const vnl_matrix<T>& A, vnl_matrix<T>& C, const double alpha=1.0, const double beta=0.0 ) 
 */
__utl_syrk_Matrix(T, syrk_VnlMatrix, Vnl, vnl_matrix<T>, rows, cols, data_block, set_size); 

template <class T> 
inline T 
InnerProduct(const vnl_vector<T>& v1, const vnl_vector<T>& v2)
{
  utlSAException(v1.size() != v2.size())(v1.size())(v2.size()).msg("vector sizes mismatch");
  return utl::cblas_dot<T>(v1.size(), v1.data_block(), 1, v2.data_block(), 1);
}

/** \f$ A = alpha*x*y'+ A \f$  */
template <class T> 
inline void 
OuterProduct(const vnl_vector<T>& v1, const vnl_vector<T>& v2, vnl_matrix<T>& mat, const double alpha=1.0)
{
  int M = v1.size(), N = v2.size();
  if (M!=mat.rows() || N!=mat.cols())
    {
    mat.set_size(M, N);
    mat.fill(0.0);
    }
  utl::cblas_ger<T>(CblasRowMajor, M, N, alpha, v1.data_block(), 1, v2.data_block(), 1, mat.data_block(), mat.cols());
}

/** \f$ A = alpha*x*x'+ A \f$  */
template <class T> 
inline void 
OuterProduct(const vnl_vector<T>& v1, vnl_matrix<T>& mat, const double alpha=1.0)
{
  int M = v1.size();
  if (M!=mat.rows() || M!=mat.cols())
    {
    mat.set_size(M, M);
    mat.fill(0.0);
    }
  utl::cblas_syr<T>(CblasRowMajor, CblasUpper, M, alpha, v1.data_block(), 1, mat.data_block(), mat.cols());
  T* data = mat.data_block();
  for ( int i = 0; i < M; ++i ) 
    for ( int j = 0; j < i; ++j ) 
      data[i*M+j] = data[j*M+i];
}

template <class T> 
inline void 
GetRow(const vnl_matrix<T>& mat, const int index, vnl_vector<T>& v1)
{
  if (v1.size()!=mat.cols())
    v1.set_size(mat.cols());
  utl::cblas_copy<T>(v1.size(), mat.data_block()+mat.cols()*index, 1, v1.data_block(), 1);
}

template <class T> 
inline void 
GetColumn(const vnl_matrix<T>& mat, const int index, vnl_vector<T>& v1)
{
  if (v1.size()!=mat.rows())
    v1.set_size(mat.rows());
  utl::cblas_copy<T>(v1.size(), mat.data_block()+index, mat.cols(), v1.data_block(), 1);
}


    /** @} */

}


#endif 

