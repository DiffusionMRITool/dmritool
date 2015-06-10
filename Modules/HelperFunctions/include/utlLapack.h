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

/** @addtogroup Math
@{ */

#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus

#ifdef UTL_USE_MKL
  #include <mkl_lapacke.h>
  #include <mkl_lapack.h>
  #include <mkl.h>
#else
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
gesvd(int matrix_order, char JOBU, char JOBVT, int M, int N, T *A, int LDA, T* S, T* U, int LDU, T* VT, int LDVT, T* superb);
template <class T> inline int 
gesdd(int matrix_order, char JOBZ, int M, int N, T *A, int LDA, T* S, T* U, int LDU, T* VT, int LDVT);

template <class T> inline int 
sytri(int matrix_order, char UPLO, int N, T *A, int LDA, const int* IPIV);
template <class T> inline int
sytrf(int matrix_order, char UPLO, int N, T* A, int LDA, int* IPIV); 	

template <class T> inline T
lange(int matrix_order, char norm, int m, int n, const T* A, int LDA);

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

template <> inline int 
sytri<double>(int matrix_order, char UPLO, int N, double *A, int LDA, const int* IPIV)
{ return LAPACKE_dsytri(matrix_order,UPLO,N,A,LDA,IPIV);}
template <> inline int 
sytri<float>(int matrix_order, char UPLO, int N, float *A, int LDA, const int* IPIV)
{ return LAPACKE_ssytri(matrix_order,UPLO,N,A,LDA,IPIV); }

template <> inline int
sytrf<double>(int matrix_order, char UPLO, int N, double* A, int LDA, int* IPIV)
{ return LAPACKE_dsytrf(matrix_order,UPLO,N,A,LDA,IPIV); }
template <> inline int
sytrf<float>(int matrix_order, char UPLO, int N, float* A, int LDA, int* IPIV)
{ return LAPACKE_ssytrf(matrix_order,UPLO,N,A,LDA,IPIV); }

template <> inline double
lange<double>(int matrix_order, char norm, int m, int n, const double* A, int LDA)
{ return LAPACKE_dlange(matrix_order,norm,m,n,A,LDA); }
template <> inline float
lange<float>(int matrix_order, char norm, int m, int n, const float* A, int LDA)
{ return LAPACKE_slange(matrix_order,norm,m,n,A,LDA); }


    /** @} */
}

#endif 
