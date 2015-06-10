/**
 *       @file  utlVNLLapack.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-13-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlVNLLapack_h
#define __utlVNLLapack_h


#include "utlCore.h"
#include "utlVNL.h"
#include "utlLapack.h"
#include "utlVNLBlas.h"

namespace utl
{

/** @addtogroup utlHelperFunctions
@{ */

/**
 * \brief syev_VnlMatrix 
 *  eigen-decomposition for symmetric matrix.
 *  http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f.html
 *
 * \param mat  a symmetric matrix (only upper triangular matrix is used)
 * \param eigenValues Eigen-values are in \b increasing order. 
 * \param eigenVectors Eigen-vectors. each row is an eigen-vector
 */
template <class T> inline void 
syev_VnlMatrix ( const vnl_matrix<T>& mat, vnl_vector<T>& eigenValues, vnl_matrix<T>& eigenVectors)
{
  utlException(mat.rows()!=mat.columns(), "The matrix should be symmetric");
  char JOBZ='V', UPLO='U';
  int N = mat.rows(),INFO=0, LWORK=-1;

  utl::MatrixCopy(mat, eigenVectors, 1.0, 'N');
  eigenValues.set_size(N);
  // T query;
  // // utl::syev<T>(&JOBZ,&UPLO,&N,eigenVectors.data_block(), &N, eigenValues.data_block(),&query,&LWORK,&INFO);
  // utl::syev<T>(LAPACK_ROW_MAJOR, JOBZ,UPLO,N,eigenVectors.data_block(), N, eigenValues.data_block(),query,LWORK,&INFO);
  // LWORK=static_cast<int>(query); 

  // T *const WORK = new T[LWORK];
  // utl::syev<T>(&JOBZ,&UPLO,&N,eigenVectors.data_block(), &N, eigenValues.data_block(),WORK,&LWORK,&INFO);
  // delete[] WORK;
  // utlGlobalException(INFO, "LAPACK library function dsyev_() returned error code INFO=" << INFO);

  INFO = utl::syev<T>(LAPACK_COL_MAJOR, JOBZ,UPLO,N,eigenVectors.data_block(), N, eigenValues.data_block());
  utlGlobalException(INFO, "LAPACK library function dsyev_() returned error code INFO=" << INFO);
}

/**
 * \brief syevd_VnlMatrix 
 *  eigen-decomposition for symmetric matrix. dsyevd is faster than dsyev
 *  http://www.netlib.org/lapack/explore-html/d1/da2/dsyevd_8f.html
 *
 * \param mat  a symmetric matrix (only upper triangular matrix is used)
 * \param eigenValues Eigen-values are in \b increasing order. 
 * \param eigenVectors Eigen-vectors. each row is an eigen-vector
 */
template <class T> inline void 
syevd_VnlMatrix ( const vnl_matrix<T>& mat, vnl_vector<T>& eigenValues, vnl_matrix<T>& eigenVectors)
{
  utlException(mat.rows()!=mat.columns(), "The matrix should be symmetric");
  char JOBZ='V', UPLO='U';
  int N = mat.rows(),INFO=0, LWORK=-1, LIWORK=-1;

  utl::MatrixCopy(mat, eigenVectors, 1.0, 'N');
  eigenValues.set_size(N);

  // T query; 
  // int queryI;
  // utl::syevd<T>(&JOBZ,&UPLO,&N,eigenVectors.data_block(), &N, eigenValues.data_block(),&query,&LWORK,&queryI,&LIWORK,&INFO);
  // LWORK=static_cast<int>(query); 
  // LIWORK=static_cast<int>(queryI); 

  // T *const WORK = new T[LWORK];
  // int *const IWORK = new int[LIWORK];
  // utl::syevd<T>(&JOBZ,&UPLO,&N,eigenVectors.data_block(), &N, eigenValues.data_block(),WORK,&LWORK,IWORK,&LIWORK,&INFO);
  // delete[] WORK;
  // delete[] IWORK;
  // utlGlobalException(INFO, "LAPACK library function dsyevd_() returned error code INFO=" << INFO);

  INFO = utl::syevd<T>(LAPACK_COL_MAJOR, JOBZ,UPLO,N,eigenVectors.data_block(), N, eigenValues.data_block());
  utlGlobalException(INFO, "LAPACK library function dsyevd_() returned error code INFO=" << INFO);
}


/**
 * \brief dgesvd_VnlMatrix 
 *
 * \param mat matrix with size MxN.
 * \param U left singular vectors. If format is 'A', U is MxM matrix. If format is 'S', U size is M x min(M,N)
 * \param s singular values with size min(M,N). Sored in \b decreasing order.
 * \param V right singular vectors. If format is 'A', V size is NxN. If format is 'S', V size is N x min(M,N)
 * \param format 'S' or 'A'. 'A' means full size, 'S' means reduced size. 
 */
template <class T> inline void 
gesvd_VnlMatrix(const vnl_matrix<T>& mat, vnl_matrix<T>& U, vnl_vector<T>& s, vnl_matrix<T>& V, char format='S' )
{
  utlException(format!='S' && format!='A', "wrong format. format should be 'A' or 'S' ");
  int M=mat.rows(), N=mat.columns(), INFO=0;
  int min_MN = utl::min(M,N), LDVT;
  s.set_size(min_MN);
  if (format=='S')
    {
    U.set_size(min_MN, M);
    V.set_size(N, min_MN);
    LDVT = min_MN;
    }
  else
    {
    U.set_size(M, M);
    V.set_size(N, N);
    LDVT = N;
    }
  char formatU=format, formatV=format;
  int LDA=M, LDU=M, LWORK=-1;
  vnl_matrix<T> matTmp;
  utl::MatrixCopy(mat, matTmp, 1.0, 'T');

  // T query;
  // // dgesdd_(&format, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, &query, &LWORK, IWORK, &INFO);
  // utl::gesvd<T>(&formatU, &formatV, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, &query, &LWORK, &INFO);
  // LWORK=static_cast<int>(query); 

  // T* WORK = new T[LWORK];
  // // dgesdd_(&format, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, WORK, &LWORK, IWORK, &INFO);
  // utl::gesvd<T>(&formatU, &formatV, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, WORK, &LWORK, &INFO);
  // delete[] WORK;

  T* superb = new T[min_MN];
  INFO=utl::gesvd<T>(LAPACK_COL_MAJOR, formatU, formatV, M, N, matTmp.data_block(), LDA, s.data_block(), U.data_block(), LDU, V.data_block(), LDVT, superb);
  delete[] superb;
  U.inplace_transpose();

  utlGlobalException(INFO, "LAPACK library function dgesvd_() returned error code INFO=" << INFO);
}

/**
 * \brief dgesdd_VnlMatrix 
 *  dgesdd is faster than dgesvd. 
 *  http://www.netlib.org/lapack/explore-html/db/db4/dgesdd_8f.html
 *  http://www.netlib.org/lapack/lug/node71.html
 *
 * \param mat matrix with size MxN.
 * \param U left singular vectors. If format is 'A', U is MxM matrix. If format is 'S', U size is M x min(M,N)
 * \param s singular values with size min(M,N). Sored in \b decreasing order.
 * \param V right singular vectors. If format is 'A', V size is NxN. If format is 'S', V size is N x min(M,N)
 * \param format 'S' or 'A'. 'A' means full size, 'S' means reduced size. 
 */
template <class T> inline void 
gesdd_VnlMatrix(const vnl_matrix<T>& mat, vnl_matrix<T>& U, vnl_vector<T>& s, vnl_matrix<T>& V, char format='S' )
{
  utlException(format!='S' && format!='A', "wrong format. format should be 'A' or 'S' ");
  int M=mat.rows(), N=mat.columns(), INFO=0;
  int min_MN = utl::min(M,N), LDVT;
  s.set_size(min_MN);
  if (format=='S')
    {
    U.set_size(min_MN, M);
    V.set_size(N, min_MN);
    LDVT = min_MN;
    }
  else
    {
    U.set_size(M, M);
    V.set_size(N, N);
    LDVT = N;
    }

  vnl_matrix<T> matTmp;
  utl::MatrixCopy(mat, matTmp, 1.0, 'T');
  char formatU=format, formatV=format;
  int LDA=M, LDU=M, LWORK;

  // int* IWORK = new int[8*min_MN];
  // T query;
  // LWORK=-1;
  // utl::gesdd<T>(&format, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, &query, &LWORK, IWORK, &INFO);
  // // dgesvd_(&formatU, &formatV, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, &query, &LWORK, &INFO);
  // LWORK=static_cast<int>(query); 

  // // LWORK = min_MN*(6+4*min_MN)+utl::max(M,N);
  // T* WORK = new T[LWORK];
  // utl::gesdd<T>(&format, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, WORK, &LWORK, IWORK, &INFO);
  // // dgesvd_(&formatU, &formatV, &M, &N, matTmp.data_block(), &LDA, s.data_block(), U.data_block(), &LDU, V.data_block(), &LDVT, WORK, &LWORK, &INFO);
  // delete[] WORK;
  // delete[] IWORK;

  INFO=utl::gesdd<T>(LAPACK_COL_MAJOR, format, M, N, matTmp.data_block(), LDA, s.data_block(), U.data_block(), LDU, V.data_block(), LDVT);
  U.inplace_transpose();

  utlGlobalException(INFO, "LAPACK library function dgesdd_() returned error code INFO=" << INFO);
}

/**
 * \brief  EigenDecompositionSymmetricVnlMatrix
 *  eigen-decomposition for symmetric matrix.
 *
 * \param mat  a symmetric matrix (only upper triangular matrix is used)
 * \param eigenValues Eigen-values are in \b increasing order. 
 * \param eigenVectors Eigen-vectors. each row is an eigen-vector
 */
template <class T> inline void 
EigenDecompositionSymmetricVnlMatrix ( const vnl_matrix<T>& mat, vnl_vector<T>& eigenValues, vnl_matrix<T>& eigenVectors )
{
#ifdef UTL_USE_FASTLAPACK
  // NOTE: fast but may be wrong for big matrix and multi-thread if openblas is not correctly built
  syevd_VnlMatrix<T>(mat, eigenValues, eigenVectors);
#else
  // slow but robust
  syev_VnlMatrix<T>(mat, eigenValues, eigenVectors);
#endif
}

/**
 * \brief SVDVnlMatrix 
 *
 * \param mat matrix with size MxN.
 * \param U left singular vectors. If format is 'A', U is MxM matrix. If format is 'S', U size is M x min(M,N)
 * \param s singular values with size min(M,N). Sored in \b decreasing order.
 * \param V right singular vectors. If format is 'A', V size is NxN. If format is 'S', V size is N x min(M,N)
 * \param format 'S' or 'A'. 'A' means full size, 'S' means reduced size. 
 */
template <class T> inline void 
SVDVnlMatrix ( const vnl_matrix<T>& mat, vnl_matrix<T>& U, vnl_vector<T>& s, vnl_matrix<T>& V, char format='S' )
{
#ifdef UTL_USE_FASTLAPACK
  // NOTE: fast but may be wrong for big matrix and multi-thread if openblas is not correctly built
  utl::gesdd_VnlMatrix<T>(mat, U, s, V, format);
#else
  // slow but robust
  utl::gesvd_VnlMatrix<T>(mat, U, s, V, format);
#endif
}

/** inverse of a symmetric matrix (non-singular). If the matrix is singular, it stop with an error. 
 * It is fast than PInverseVnlMatrix, but only works for non-singular matrix. */
template <class T> inline void 
InverseSymmericVnlMatrix( const vnl_matrix<T>& mat, vnl_matrix<T>& result, const T eps=1e-8 )
{
  int n = mat.rows();
  utl::MatrixCopy(mat, result, 1.0, 'N');

  char uplo='U';
  int lwork=-1, INFO=0;
  int* ipiv= new int[n];

  // T *query, *work;
  // query = new T[1];
  // utl::sytrf<T>(&uplo,&n,result.data_block(),&n,ipiv,query,&lwork,&INFO);
  // lwork=static_cast<INTT>(*query); 
  // delete[] query;
  // work = new T[lwork];
  // utl::sytrf<T>(&uplo,&n,result.data_block(),&n,ipiv,work,&lwork,&INFO);
  // delete[] work;
  // work = new T[2*n];
  // utl::sytri<T>(&uplo,&n,result.data_block(),&n,ipiv,work,&INFO);
  // delete[] work;
  
  utl::sytrf<T>(LAPACK_COL_MAJOR, uplo,n,result.data_block(),n,ipiv);
  INFO=utl::sytri<T>(LAPACK_COL_MAJOR, uplo,n,result.data_block(),n,ipiv);
  delete[] ipiv;

  utlGlobalException(INFO>0, "The matrix is singular and its inverse could not be computed. \
    LAPACK library function sytrid_() returned error code INFO=" << INFO);
  utlGlobalException(INFO<0, "LAPACK library function sytrid_() returned error code INFO=" << INFO);

  T* data = result.data_block();
  for ( int i = 0; i < n; ++i ) 
    for ( int j = 0; j < i; ++j ) 
       data[j*n+i] = data[i*n+j];
}

/** pseudo-inverse of a symmetric matrix which can be singular. If the matrix is not singular, it returns the inverse.  */
template <class T> inline void 
PInverseSymmericVnlMatrix( const vnl_matrix<T>& mat, vnl_matrix<T>& result, const T eps=1e-8 )
{
  int N = mat.rows();
  vnl_matrix<T> eigenVectors, tmp, S(N,N,0.0);
  vnl_vector<T> eigenValues, v;
  EigenDecompositionSymmetricVnlMatrix(mat, eigenValues, eigenVectors);
  for (unsigned i=0; i<N; ++i)
    {
    if ( eigenValues[i]>eps || eigenValues[i]<-eps ) 
      S(i,i) = 1.0/eigenValues[i];
    else
      S(i,i)=0.0;
    }
  utl::ProductVnlMtM(eigenVectors, S, tmp);
  utl::ProductVnlMM(tmp, eigenVectors, result);
}

/** pseudo-inverse of a general matrix.  */
template <class T> inline void 
PInverseVnlMatrix( const vnl_matrix<T>& mat, vnl_matrix<T>& result, const T eps=1e-8 )
{
  vnl_matrix<T> U,V, tmp;
  vnl_vector<T> S, uVec, vVec;
  SVDVnlMatrix(mat, U, S, V, 'S');

  vnl_matrix<T> diag(S.size(), S.size(),0.0);
  for ( int i = 0; i < S.size(); i += 1 )
    {
    if ( S[i]>eps || S[i]<-eps ) 
      diag(i,i) = 1.0/S[i];
    else
      diag(i,i) = 0.0;
    }
  utl::ProductVnlMM(V, diag, tmp);
  utl::ProductVnlMMt(tmp, U, result);

  // result.set_size(mat.cols(), mat.rows());
  // result.fill(0.0);
  // vnl_vector<T> vV(V.rows()), vU(U.cols());
  // for ( int i = 0; i < S.size(); ++i ) 
  //   {
  //   if ( S[i]>eps || S[i]<-eps ) 
  //     {
  //     utl::GetColumn(V,i,vV);
  //     utl::GetColumn(U,i,vU);
  //     utl::OuterProductvv(vV, vU, result, 1.0/S[i]);
  //     }
  //   }
}

/** The projection onto the plane Aeq^T*x=beq  */
template <class T>
void 
GetEqualityConstraintProjection ( const vnl_matrix<T>& Aeq, const vnl_vector<T>& beq, const vnl_matrix<T>& QInverse, 
  vnl_matrix<T>& projMatrix, vnl_vector<T>& projVector )
{
  int n = QInverse.rows();
  utlException(Aeq.rows()!=n, "wrong size! Aeq.rows()="<<Aeq.rows()<<", n="<<n);
  vnl_matrix<T> AeqT=Aeq.transpose();
  projMatrix.set_size(n,n);
  projMatrix.set_identity();
  vnl_matrix<T> temp, tmp2, tmp3;
  ProductVnlMM(AeqT, QInverse, Aeq, tmp2);
  PInverseSymmericVnlMatrix(tmp2, temp);
  ProductVnlMM( QInverse, Aeq, temp, tmp2);
  ProductVnlMM( tmp2, AeqT, tmp3);
  projMatrix -= tmp3;

  ProductVnlMv( tmp2, beq, projVector);
}
    /** @} */

}


#endif 
