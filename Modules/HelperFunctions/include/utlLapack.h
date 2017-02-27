/**
 *       @file  utlLapack.h
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

#ifndef __utlLapack_h
#define __utlLapack_h

/** @addtogroup utlMath
@{ */

#include "utlCoreMacro.h"

#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus

#ifdef UTL_USE_MKL
  #include "utlMKL.h"
  #include <mkl_lapacke.h>
  #include <mkl_lapack.h>
  #include <mkl.h>
  #include <complex>
#else

  #define lapack_complex_double std::complex<double>
  #define lapack_complex_float std::complex<float>
  #include <lapacke.h>
#endif


// #if defined(UTL_USE_MKL)
//   #include <mkl_lapack.h>
// #else
// #if !defined(_MKL_LAPACK_H_) && !defined(_LAPACKE_H_)

//   [>* http:www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f.html  <]
//   extern void dsyev_(char* JOBZ, char* UPLO, int* N, double *A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);
//   extern void ssyev_(char* JOBZ, char* UPLO, int* N, float *A, int* LDA, float* W, float* WORK, int* LWORK, int* INFO);
//   [>* http:www.netlib.org/lapack/explore-html/d1/da2/dsyevd_8f.html  <]
//   extern void dsyevd_(char* JOBZ, char* UPLO, int* N, double *A, int* LDA, double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
//   extern void ssyevd_(char* JOBZ, char* UPLO, int* N, float *A, int* LDA, float* W, float* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

//   [>* http:www.netlib.org/lapack/explore-html/d8/d2d/dgesvd_8f.html  <] 
//   extern void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double *A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);
//   extern void sgesvd_(char* JOBU, char* JOBVT, int* M, int* N, float *A, int* LDA, float* S, float* U, int* LDU, float* VT, int* LDVT, float* WORK, int* LWORK, int* INFO);
//   [>* http:www.netlib.org/lapack/explore-html/db/db4/dgesdd_8f.html  <]
//   extern void dgesdd_(char* JOBZ, int* M, int* N, double *A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO);
//   extern void sgesdd_(char* JOBZ, int* M, int* N, float *A, int* LDA, float* S, float* U, int* LDU, float* VT, int* LDVT, float* WORK, int* LWORK, int* IWORK, int* INFO);

//   [>* http:www.netlib.org/lapack/explore-html/d5/d8f/dsytri_8f.html  <]
//   extern void dsytri_(char* UPLO, int* N, double *A, int* LDA, int* IPIV, double* WORK, int* INFO);
//   extern void ssytri_(char* UPLO, int* N, float *A, int* LDA, int* IPIV, float* WORK, int* INFO);

//   [>* http:www.netlib.org/lapack/explore-html/dd/df4/dsytrf_8f.html  <]
//   extern void dsytrf_(char* UPLO, int* N, double *A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
//   extern void ssytrf_(char* UPLO, int* N, float *A, int* LDA, int* IPIV, float* WORK, int* LWORK, int* INFO);
// #endif

// #endif

float LAPACKE_slange( int matrix_order, char norm, lapack_int m, lapack_int n, const float* a, lapack_int lda );
double LAPACKE_dlange( int matrix_order, char norm, lapack_int m, lapack_int n, const double* a, lapack_int lda );
double LAPACKE_zlange( int matrix_order, char norm, lapack_int m, lapack_int n, const std::complex<double>* a, lapack_int lda );
float LAPACKE_clange( int matrix_order, char norm, lapack_int m, lapack_int n, const std::complex<float>* a, lapack_int lda );


#ifdef __cplusplus
  }
#endif /* __cplusplus */

namespace utl
{

// [>************************************************************************************<]
// [>* template function definitions  <]

// template <class T> inline void 
// syev(char* JOBZ, char* UPLO, int* N, T *A, int* LDA, T* W, T* WORK, int* LWORK, int* INFO);
// template <class T> inline void 
// syevd(char* JOBZ, char* UPLO, int* N, T *A, int* LDA, T* W, T* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

// template <class T> inline void 
// gesvd(char* JOBU, char* JOBVT, int* M, int* N, T *A, int* LDA, T* S, T* U, int* LDU, T* VT, int* LDVT, T* WORK, int* LWORK, int* INFO);
// template <class T> inline void 
// gesdd(char* JOBZ, int* M, int* N, T *A, int* LDA, T* S, T* U, int* LDU, T* VT, int* LDVT, T* WORK, int* LWORK, int* IWORK, int* INFO);

// template <class T> inline void 
// sytri(char* UPLO, int* N, T *A, int* LDA, int* IPIV, T* WORK, int* INFO);

// template <class T> inline void 
// sytrf(char* UPLO, int* N, T* A, int* LDA, int* IPIV, T* WORK, int* LWORK, int* INFO); 	

// template <class T> inline T
// lange(char* NORM, int* M, int* N, const T* A, int* LDA, T* WORK); 	



// // ************************************************************************************
// [>* function implementations using double and float  <]

// template <> inline void 
// syev<double>(char* JOBZ, char* UPLO, int* N, double *A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO)
// { dsyev_(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO);}
// template <> inline void 
// syev<float>(char* JOBZ, char* UPLO, int* N, float *A, int* LDA, float* W, float* WORK, int* LWORK, int* INFO)
// { ssyev_(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO);}

// template <> inline void 
// syevd<double>(char* JOBZ, char* UPLO, int* N, double *A, int* LDA, double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO)
// { dsyevd_(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO);}
// template <> inline void 
// syevd<float>(char* JOBZ, char* UPLO, int* N, float *A, int* LDA, float* W, float* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO)
// { ssyevd_(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO);}

// template <> inline void 
// gesvd<double>(char* JOBU, char* JOBVT, int* M, int* N, double *A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO)
// { dgesvd_(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO);}
// template <> inline void 
// gesvd<float>(char* JOBU, char* JOBVT, int* M, int* N, float *A, int* LDA, float* S, float* U, int* LDU, float* VT, int* LDVT, float* WORK, int* LWORK, int* INFO)
// { sgesvd_(JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO);}

// template <> inline void 
// gesdd<double>(char* JOBZ, int* M, int* N, double *A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO)
// { dgesdd_(JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,IWORK,INFO);}
// template <> inline void 
// gesdd<float>(char* JOBZ, int* M, int* N, float *A, int* LDA, float* S, float* U, int* LDU, float* VT, int* LDVT, float* WORK, int* LWORK, int* IWORK, int* INFO)
// { sgesdd_(JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,IWORK,INFO);}

// template <> inline void 
// sytri<double>(char* UPLO, int* N, double *A, int* LDA, int* IPIV, double* WORK, int* INFO)
// { dsytri_(UPLO,N,A,LDA,IPIV,WORK,INFO);}
// template <> inline void 
// sytri<float>(char* UPLO, int* N, float *A, int* LDA, int* IPIV, float* WORK, int* INFO)
// { ssytri_(UPLO,N,A,LDA,IPIV,WORK,INFO);}

// template <class T> inline void 
// sytrf(char* UPLO, int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO)
// { dsytrf_(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO);}

// template <class T> inline void 
// sytrf(char* UPLO, int* N, float* A, int* LDA, int* IPIV, float* WORK, int* LWORK, int* INFO)	
// { ssytrf_(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO);}

// template <> inline double 
// lange<double>( char* norm, int* m, int* n, const double* a, int* lda, double* work )
// { return dlange_(norm,m,n,a,lda,work); }
// template <> inline float
// lange<float>( char* norm, int* m, int* n, const float* a, int* lda, float* work )
// { return slange_(norm,m,n,a,lda,work); }


/**************************************************************************************/
/** template function definitions  */

template <class T> inline int 
syev(int matrix_order, char JOBZ, char UPLO, int N, T *A, int LDA, T* W);
template <class T> inline int 
syevd(int matrix_order, char JOBZ, char UPLO, int N, T *A, int LDA, T* W);
template <class T> inline int 
geev( int matrix_layout, char jobvl, char jobvr, int n, T* a, int lda, T* wr, T* wi, T* vl, int ldvl, T* vr, int ldvr );

template <class T> inline int 
gesvd(int matrix_order, char JOBU, char JOBVT, int M, int N, T *A, int LDA, T* S, T* U, int LDU, T* VT, int LDVT, T* superb);
template <class T> inline int 
gesdd(int matrix_order, char JOBZ, int M, int N, T *A, int LDA, T* S, T* U, int LDU, T* VT, int LDVT);

template <class T> inline int 
sytri(int matrix_order, char UPLO, int N, T *A, int LDA, const int* IPIV);
template <class T> inline int
sytrf(int matrix_order, char UPLO, int N, T* A, int LDA, int* IPIV); 	

template <class T> inline T
lange(int matrix_order, char norm, int m, int n, const T* A, int LDA);


template <class T> inline int
getrf (int matrix_layout , int m , int n , T * a , int lda, int * ipiv );

template <class T> inline int
getri(int matrix_layout , int n , T * a , int lda , const int * ipiv);

/**************************************************************************************/
/** function implementations using double and float  */

template <> inline int 
syev<double>(int matrix_order, char JOBZ, char UPLO, int N, double *A, int LDA, double* W)
{ return LAPACKE_dsyev(matrix_order,JOBZ,UPLO,N,A,LDA,W);}
template <> inline int 
syev<float>(int matrix_order, char JOBZ, char UPLO, int N, float *A, int LDA, float* W)
{ return LAPACKE_ssyev(matrix_order,JOBZ,UPLO,N,A,LDA,W);}

template <> inline int 
syevd<double>(int matrix_order, char JOBZ, char UPLO, int N, double *A, int LDA, double* W)
{ return LAPACKE_dsyevd(matrix_order,JOBZ,UPLO,N,A,LDA,W);}
template <> inline int
syevd<float>(int matrix_order, char JOBZ, char UPLO, int N, float *A, int LDA, float* W)
{ return LAPACKE_ssyevd(matrix_order,JOBZ,UPLO,N,A,LDA,W);}
template <> inline int 
syevd<std::complex<double> >(int matrix_order, char JOBZ, char UPLO, int N, std::complex<double> *A, int LDA, std::complex<double>* W)
{ utlGlobalException(true, "connot use complex values");  return 0;}
template <> inline int 
syevd<std::complex<float> >(int matrix_order, char JOBZ, char UPLO, int N, std::complex<float> *A, int LDA, std::complex<float>* W)
{ utlGlobalException(true, "connot use complex values");  return 0;}

template <> inline int 
geev<double>( int matrix_layout, char jobvl, char jobvr, int n, double* a, int lda, double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr )
{ return LAPACKE_dgeev( matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr ); }
template <> inline int 
geev<float>( int matrix_layout, char jobvl, char jobvr, int n, float* a, int lda, float* wr, float* wi, float* vl, int ldvl, float* vr, int ldvr )
{ return LAPACKE_sgeev( matrix_layout, jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr ); }
template <> inline int 
geev<std::complex<double> >( int matrix_layout, char jobvl, char jobvr, int n, std::complex<double>* a, int lda, std::complex<double>* wr, std::complex<double>* wi, std::complex<double>* vl, int ldvl, std::complex<double>* vr, int ldvr )
{ utlGlobalException(true, "connot use complex values"); return 0; }
template <> inline int 
geev<std::complex<float> >( int matrix_layout, char jobvl, char jobvr, int n, std::complex<float>* a, int lda, std::complex<float>* wr, std::complex<float>* wi, std::complex<float>* vl, int ldvl, std::complex<float>* vr, int ldvr )
{ utlGlobalException(true, "connot use complex values"); return 0; }
inline int 
geev( int matrix_layout, char jobvl, char jobvr, int n, std::complex<double>* a, int lda, std::complex<double>* w, std::complex<double>* vl, int ldvl, std::complex<double>* vr, int ldvr )
{ return LAPACKE_zgeev( matrix_layout, jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr ); }

template <> inline int 
gesvd<double>(int matrix_order, char JOBU, char JOBVT, int M, int N, double *A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double* superb)
{ return LAPACKE_dgesvd(matrix_order,JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,superb);}
template <> inline int 
gesvd<float>(int matrix_order, char JOBU, char JOBVT, int M, int N, float *A, int LDA, float* S, float* U, int LDU, float* VT, int LDVT, float* superb)
{ return LAPACKE_sgesvd(matrix_order,JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,superb);}
  
template <> inline int 
gesdd<double>(int matrix_order, char JOBZ, int M, int N, double *A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT)
{ return LAPACKE_dgesdd(matrix_order,JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT);}
template <> inline int 
gesdd<float>(int matrix_order, char JOBZ, int M, int N, float *A, int LDA, float* S, float* U, int LDU, float* VT, int LDVT)
{ return LAPACKE_sgesdd(matrix_order,JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT);}
inline int 
gesdd(int matrix_order, char JOBZ, int M, int N, std::complex<double> *A, int LDA, double* S, std::complex<double>* U, int LDU, std::complex<double>* VT, int LDVT)
{ return LAPACKE_zgesdd(matrix_order,JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT);}
inline int 
gesdd(int matrix_order, char JOBZ, int M, int N, std::complex<float> *A, int LDA, float* S, std::complex<float>* U, int LDU, std::complex<float>* VT, int LDVT)
{ return LAPACKE_cgesdd(matrix_order,JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT);}

template <> inline int 
sytri<double>(int matrix_order, char UPLO, int N, double *A, int LDA, const int* IPIV)
{ return LAPACKE_dsytri(matrix_order,UPLO,N,A,LDA,IPIV);}
template <> inline int 
sytri<float>(int matrix_order, char UPLO, int N, float *A, int LDA, const int* IPIV)
{ return LAPACKE_ssytri(matrix_order,UPLO,N,A,LDA,IPIV); }
template <> inline int 
sytri<std::complex<double> >(int matrix_order, char UPLO, int N, std::complex<double> *A, int LDA, const int* IPIV)
{ return LAPACKE_zsytri(matrix_order,UPLO,N,A,LDA,IPIV);}
template <> inline int 
sytri<std::complex<float> >(int matrix_order, char UPLO, int N, std::complex<float> *A, int LDA, const int* IPIV)
{ return LAPACKE_csytri(matrix_order,UPLO,N,A,LDA,IPIV);}

template <> inline int
sytrf<double>(int matrix_order, char UPLO, int N, double* A, int LDA, int* IPIV)
{ return LAPACKE_dsytrf(matrix_order,UPLO,N,A,LDA,IPIV); }
template <> inline int
sytrf<float>(int matrix_order, char UPLO, int N, float* A, int LDA, int* IPIV)
{ return LAPACKE_ssytrf(matrix_order,UPLO,N,A,LDA,IPIV); }
template <> inline int
sytrf<std::complex<double> >(int matrix_order, char UPLO, int N, std::complex<double>* A, int LDA, int* IPIV)
{ return LAPACKE_zsytrf(matrix_order,UPLO,N,A,LDA,IPIV); }
template <> inline int
sytrf<std::complex<float> >(int matrix_order, char UPLO, int N, std::complex<float>* A, int LDA, int* IPIV)
{ return LAPACKE_csytrf(matrix_order,UPLO,N,A,LDA,IPIV); }

template <> inline double
lange<double>(int matrix_order, char norm, int m, int n, const double* A, int LDA)
{ return LAPACKE_dlange(matrix_order,norm,m,n,A,LDA); }
template <> inline float
lange<float>(int matrix_order, char norm, int m, int n, const float* A, int LDA)
{ return LAPACKE_slange(matrix_order,norm,m,n,A,LDA); }
inline double
lange(int matrix_order, char norm, int m, int n, const std::complex<double>* A, int LDA)
{ return LAPACKE_zlange(matrix_order,norm,m,n,A,LDA); }
inline float
lange(int matrix_order, char norm, int m, int n, const std::complex<float>* A, int LDA)
{ return LAPACKE_clange(matrix_order,norm,m,n,A,LDA); }

template <> inline int
getrf<double>(int matrix_layout , int m , int n , double * a , int lda, int * ipiv )
{ return LAPACKE_dgetrf(matrix_layout, m, n, a, lda, ipiv); }
template <> inline int
getrf<float>(int matrix_layout , int m , int n , float * a , int lda, int * ipiv )
{ return LAPACKE_sgetrf(matrix_layout, m, n, a, lda, ipiv); }
template <> inline int
getrf<std::complex<double> >(int matrix_layout , int m , int n , std::complex<double> * a , int lda, int * ipiv )
{ return LAPACKE_zgetrf(matrix_layout, m, n, a, lda, ipiv); }
template <> inline int
getrf<std::complex<float> >(int matrix_layout , int m , int n , std::complex<float> * a , int lda, int * ipiv )
{ return LAPACKE_cgetrf(matrix_layout, m, n, a, lda, ipiv); }

template <> inline int
getri<double>(int matrix_layout , int n , double * a , int lda , const int * ipiv)
{ return LAPACKE_dgetri(matrix_layout, n, a, lda, ipiv); }
template <> inline int
getri<float>(int matrix_layout , int n , float * a , int lda , const int * ipiv)
{ return LAPACKE_sgetri(matrix_layout, n, a, lda, ipiv); }
template <> inline int
getri<std::complex<double> >(int matrix_layout , int n , std::complex<double> * a , int lda , const int * ipiv)
{ return LAPACKE_zgetri(matrix_layout, n, a, lda, ipiv); }
template <> inline int
getri<std::complex<float> >(int matrix_layout , int n , std::complex<float> * a , int lda , const int * ipiv)
{ return LAPACKE_cgetri(matrix_layout, n, a, lda, ipiv); }


#define __utl_getri_Matrix(T, FuncName, RowMajorMatrixName, GetRowsFuncName, GetColsFuncName, MatrixGetDataFuncName, ReSizeFuncName)                                                \
template <class T> inline void                                                                                                                                                      \
FuncName(const RowMajorMatrixName& mat, RowMajorMatrixName& result)                                                                                                                 \
{                                                                                                                                                                                   \
  int N = mat.GetRowsFuncName(),INFO=0;                                                                                                                                             \
  int ipiv[N];                                                                                                                                                                      \
  result = mat;                                                                                                                                                                     \
  INFO = getrf<T>(LAPACK_ROW_MAJOR, N,N,result.MatrixGetDataFuncName(),N,ipiv);                                                                                                     \
  INFO = getri<T>(LAPACK_ROW_MAJOR, N,result.MatrixGetDataFuncName(),N,ipiv);                                                                                                       \
}


/** 
 *
 * \brief geev_VnlMatrix 
 *  Calculate non-symmetric eigen-decomposition. 
 *  Define 3 functions.
 *  1) calculate only eigenvalues. 
 *  2) eigenvalues, right eigenvectors.
 *  3) eigenvalues, right eigenvectors, left eigenvectors.
 *
 * \param mat matrix with size NxN.
 * \param valReal real part of right eigen-values.
 * \param valImg imginary part of right eigen-values.
 * \param vecRealR real part of right eigen-vectors.
 * \param vecImgR part of right eigen-vectors.
 * \param vecRealL real part of left eigen-vectors.
 * \param vecImgL part of left eigen-vectors.
 *
 * http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga8ec1625302675b981eb34ed024b27a47.html 
 * http://www.netlib.org/lapack/lug/node31.html
 *
 *
 * */
#define __utl_geev_Matrix(T, FuncName, RowMajorMatrixName, GetRowsFuncName, GetColsFuncName, MatrixGetDataFuncName, VectorName, VectorGetDataFuncName, ReSizeFuncName)                                                                        \
template <class T> inline void                                                                                                                                                                                                                \
FuncName ( const RowMajorMatrixName& mat, VectorName& valReal, VectorName& valImg, RowMajorMatrixName& vecRealR, RowMajorMatrixName& vecImgR, RowMajorMatrixName& vecRealL, RowMajorMatrixName& vecImgL)                                      \
{                                                                                                                                                                                                                                             \
  utlException(mat.GetRowsFuncName()!=mat.GetColsFuncName(), "The matrix should be square");                                                                                                                                                  \
  char JOBVL='V', JOBVR='V';                                                                                                                                                                                                                  \
  int N = mat.GetRowsFuncName(),INFO=0;                                                                                                                                                                                                       \
  int LDA=N, LDVL=N, LDVR=N;                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                              \
  RowMajorMatrixName matCopy(mat), vecR(N,N), vecL(N,N);                                                                                                                                                                                      \
  valReal.ReSizeFuncName(N);                                                                                                                                                                                                                  \
  valImg.ReSizeFuncName(N);                                                                                                                                                                                                                   \
  vecRealR.ReSizeFuncName(N,N);                                                                                                                                                                                                               \
  vecImgR.ReSizeFuncName(N,N);                                                                                                                                                                                                                \
  vecRealL.ReSizeFuncName(N,N);                                                                                                                                                                                                               \
  vecImgL.ReSizeFuncName(N,N);                                                                                                                                                                                                                \
                                                                                                                                                                                                                                              \
  INFO = geev<T>( LAPACK_ROW_MAJOR, JOBVL, JOBVR, N, matCopy.MatrixGetDataFuncName(), LDA, valReal.VectorGetDataFuncName(), valImg.VectorGetDataFuncName(), vecL.MatrixGetDataFuncName(), LDVL, vecR.MatrixGetDataFuncName(), LDVR );         \
                                                                                                                                                                                                                                              \
  T const * vecR_data=vecR.MatrixGetDataFuncName();                                                                                                                                                                                           \
  T* vecRealR_data=vecRealR.MatrixGetDataFuncName();                                                                                                                                                                                          \
  T* vecImgR_data=vecImgR.MatrixGetDataFuncName();                                                                                                                                                                                            \
  T const * vecL_data=vecL.MatrixGetDataFuncName();                                                                                                                                                                                           \
  T* vecRealL_data=vecRealL.MatrixGetDataFuncName();                                                                                                                                                                                          \
  T* vecImgL_data=vecImgL.MatrixGetDataFuncName();                                                                                                                                                                                            \
  for ( int i = 0; i < N; ++i )                                                                                                                                                                                                               \
    {                                                                                                                                                                                                                                         \
    int j=0;                                                                                                                                                                                                                                  \
    int k0 = i*LDVR;                                                                                                                                                                                                                          \
    while( j < N )                                                                                                                                                                                                                            \
      {                                                                                                                                                                                                                                       \
      int kk=k0+j;                                                                                                                                                                                                                            \
      if( valImg[j] == (double)0.0 )                                                                                                                                                                                                          \
        {                                                                                                                                                                                                                                     \
        vecRealR_data[kk] = vecR_data[kk];                                                                                                                                                                                                    \
        vecImgR_data[kk] = 0.0;                                                                                                                                                                                                               \
        vecRealL_data[kk] = vecL_data[kk];                                                                                                                                                                                                    \
        vecImgL_data[kk] = 0.0;                                                                                                                                                                                                               \
        j++;                                                                                                                                                                                                                                  \
        }                                                                                                                                                                                                                                     \
      else                                                                                                                                                                                                                                    \
        {                                                                                                                                                                                                                                     \
        vecRealR_data[kk] = vecR_data[kk];                                                                                                                                                                                                    \
        vecImgR_data[kk] = vecR_data[kk+1];                                                                                                                                                                                                   \
                                                                                                                                                                                                                                              \
        vecRealR_data[kk+1] = vecR_data[kk];                                                                                                                                                                                                  \
        vecImgR_data[kk+1] = -vecR_data[kk+1];                                                                                                                                                                                                \
                                                                                                                                                                                                                                              \
        vecRealL_data[kk] = vecL_data[kk];                                                                                                                                                                                                    \
        vecImgL_data[kk] = vecL_data[kk+1];                                                                                                                                                                                                   \
                                                                                                                                                                                                                                              \
        vecRealL_data[kk+1] = vecL_data[kk];                                                                                                                                                                                                  \
        vecImgL_data[kk+1] = -vecL_data[kk+1];                                                                                                                                                                                                \
                                                                                                                                                                                                                                              \
        j += 2;                                                                                                                                                                                                                               \
        }                                                                                                                                                                                                                                     \
      }                                                                                                                                                                                                                                       \
    }                                                                                                                                                                                                                                         \
}                                                                                                                                                                                                                                             \
                                                                                                                                                                                                                                              \
template <class T> inline void                                                                                                                                                                                                                \
FuncName ( const RowMajorMatrixName& mat, VectorName& valReal, VectorName& valImg)                                                                                                                                                            \
{                                                                                                                                                                                                                                             \
  utlException(mat.GetRowsFuncName()!=mat.GetColsFuncName(), "The matrix should be square");                                                                                                                                                  \
  char JOBVL='N', JOBVR='N';                                                                                                                                                                                                                  \
  int N = mat.GetRowsFuncName(),INFO=0;                                                                                                                                                                                                       \
  int LDA=N, LDVL=N, LDVR=N;                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                              \
  RowMajorMatrixName matCopy(mat);                                                                                                                                                                                                            \
  valReal.ReSizeFuncName(N);                                                                                                                                                                                                                  \
  valImg.ReSizeFuncName(N);                                                                                                                                                                                                                   \
                                                                                                                                                                                                                                              \
  INFO = geev<T>( LAPACK_ROW_MAJOR, JOBVL, JOBVR, N, matCopy.MatrixGetDataFuncName(), LDA, valReal.VectorGetDataFuncName(), valImg.VectorGetDataFuncName(), NULL, LDVL, NULL, LDVR );                                                         \
}                                                                                                                                                                                                                                             \
                                                                                                                                                                                                                                              \
template <class T> inline void                                                                                                                                                                                                                \
FuncName ( const RowMajorMatrixName& mat, VectorName& valReal, VectorName& valImg, RowMajorMatrixName& vecRealR, RowMajorMatrixName& vecImgR)                                                                                                 \
{                                                                                                                                                                                                                                             \
  utlException(mat.GetRowsFuncName()!=mat.GetColsFuncName(), "The matrix should be square");                                                                                                                                                  \
  char JOBVL='N', JOBVR='V';                                                                                                                                                                                                                  \
  int N = mat.GetRowsFuncName(),INFO=0;                                                                                                                                                                                                       \
  int LDA=N, LDVL=N, LDVR=N;                                                                                                                                                                                                                  \
                                                                                                                                                                                                                                              \
  RowMajorMatrixName matCopy(mat), eigenVectorCopy(N,N);                                                                                                                                                                                      \
  valReal.ReSizeFuncName(N);                                                                                                                                                                                                                  \
  valImg.ReSizeFuncName(N);                                                                                                                                                                                                                   \
  vecRealR.ReSizeFuncName(N,N);                                                                                                                                                                                                               \
  vecImgR.ReSizeFuncName(N,N);                                                                                                                                                                                                                \
                                                                                                                                                                                                                                              \
  INFO = geev<T>( LAPACK_ROW_MAJOR, JOBVL, JOBVR, N, matCopy.MatrixGetDataFuncName(), LDA, valReal.VectorGetDataFuncName(), valImg.VectorGetDataFuncName(), NULL, LDVL, eigenVectorCopy.MatrixGetDataFuncName(), LDVR );                      \
                                                                                                                                                                                                                                              \
  T const * vecR_data=eigenVectorCopy.MatrixGetDataFuncName();                                                                                                                                                                                \
  T* vecRealR_data=vecRealR.MatrixGetDataFuncName();                                                                                                                                                                                          \
  T* vecImgR_data=vecImgR.MatrixGetDataFuncName();                                                                                                                                                                                            \
  for ( int i = 0; i < N; ++i )                                                                                                                                                                                                               \
    {                                                                                                                                                                                                                                         \
    int j=0;                                                                                                                                                                                                                                  \
    int k0 = i*LDVR;                                                                                                                                                                                                                          \
    while( j < N )                                                                                                                                                                                                                            \
      {                                                                                                                                                                                                                                       \
      int kk=k0+j;                                                                                                                                                                                                                            \
      if( valImg[j] == (double)0.0 )                                                                                                                                                                                                          \
        {                                                                                                                                                                                                                                     \
        vecRealR_data[kk] = vecR_data[kk];                                                                                                                                                                                                    \
        vecImgR_data[kk] = 0.0;                                                                                                                                                                                                               \
        j++;                                                                                                                                                                                                                                  \
        }                                                                                                                                                                                                                                     \
      else                                                                                                                                                                                                                                    \
        {                                                                                                                                                                                                                                     \
        vecRealR_data[kk] = vecR_data[kk];                                                                                                                                                                                                    \
        vecImgR_data[kk] = vecR_data[kk+1];                                                                                                                                                                                                   \
                                                                                                                                                                                                                                              \
        vecRealR_data[kk+1] = vecR_data[kk];                                                                                                                                                                                                  \
        vecImgR_data[kk+1] = -vecR_data[kk+1];                                                                                                                                                                                                \
        j += 2;                                                                                                                                                                                                                               \
        }                                                                                                                                                                                                                                     \
      }                                                                                                                                                                                                                                       \
    }                                                                                                                                                                                                                                         \
}                                                                                                                                                                                                                                             \


    /** @} */
}

#endif 
