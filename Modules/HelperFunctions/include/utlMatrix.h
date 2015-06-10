/**
 *       @file  utlMatrix.h
 *      @brief  utl::NDArray<T,2> class which uses blas mkl
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlMatrix_h
#define __utlMatrix_h

#include "utlVector.h"
#include "utlCore.h"
#include "utlLapack.h"
#include "utlMath.h"

namespace utl
{

/** \ingroup utlNDArray
* @{ */

/** http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ C = alpha * A * B + beta * C \f$
*  \note: C should be pre-allocated
*
*  define several functions. 
* 
*  template \<class T\> inline bool 
*  gemm_UtlMatrixTimesMatrix(const bool bATrans, const bool bBTrans, const T alpha, const utl::NDArray<T,2>& a, const utl::NDArray<T,2>& b, const T beta, utl::NDArray<T,2>& c);
*
*  \f$ C = alpha * A * B + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMM(const utl::NDArray<T,2>& A, const utl::NDArray<T,2>& B, utl::NDArray<T,2>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A * B^T + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMMt(const utl::NDArray<T,2>& A, const utl::NDArray<T,2>& B, utl::NDArray<T,2>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A^T * B + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMtM(const utl::NDArray<T,2>& A, const utl::NDArray<T,2>& B, utl::NDArray<T,2>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A^T * B^T + beta * C \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMtMt(const utl::NDArray<T,2>& A, const utl::NDArray<T,2>& B, utl::NDArray<T,2>& C, const double alpha=1.0, const double beta=0.0);
*
*  \f$ C = alpha * A1*A2*A3 \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMM(const utl::NDArray<T,2>& A1, const utl::NDArray<T,2>& A2, const utl::NDArray<T,2>& A3, utl::NDArray<T,2>& C);
*
*  \f$ C = alpha * A1*A2*A3*A4 \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMM(const utl::NDArray<T,2>& A1, const utl::NDArray<T,2>& A2, const utl::NDArray<T,2>& A3, const utl::NDArray<T,2>& A4, utl::NDArray<T,2>& C);
*
*  \f$ C = alpha * A1*A2*A3*A5 \f$ \n
*  template \<class T\> inline void 
*  ProductUtlMM(const utl::NDArray<T,2>& A1, const utl::NDArray<T,2>& A2, const utl::NDArray<T,2>& A3, const utl::NDArray<T,2>& A4, const utl::NDArray<T,2>& A5, utl::NDArray<T,2>& C); 
*
* */
__utl_gemm_MatrixTimesMatrix(T, gemm_UtlMatrixTimesMatrix, Utl, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, ReSize);


/**
 * \brief syrk_UtlMatrix 
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
 *  syrk_UtlMatrix( const bool trans, const T alpha, const utl::NDArray<T,2>& A, const T beta, utl::NDArray<T,2>& C ) 
 *
*  \f$ C = alpha *A*A^T  + beta*C \f$ \n
 *  template \<class T\> void  
 *  ProductUtlXXt ( const vnl_matrix<T>& A, utl::NDArray<T,2>& C, const double alpha=1.0, const double beta=0.0 ) 
 *
*  \f$ C = alpha *A^T*A  + beta*C \f$ \n
 *  template \<class T\> void  
 *  ProductUtlXtX ( const vnl_matrix<T>& A, utl::NDArray<T,2>& C, const double alpha=1.0, const double beta=0.0 ) 
 */
__utl_syrk_Matrix(T, syrk_UtlMatrix, Utl, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, ReSize); 


/**
 * \brief dsyev_VnlMatrix 
 *  eigen-decomposition for symmetric matrix.
 *  http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f.html
 *
 * \param mat  a symmetric matrix (only upper triangular matrix is used)
 * \param eigenValues Eigen-values are in \b increasing order. 
 * \param eigenVectors Eigen-vectors. each row is an eigen-vector
 */
template <class T> inline void 
syev_UtlMatrix ( const NDArray<T,2>& mat, NDArray<T,1>& eigenValues, NDArray<T,2>& eigenVectors)
{
  utlException(mat.Rows()!=mat.Cols(), "The matrix should be symmetric");
  char JOBZ='V', UPLO='U';
  int N = mat.Rows(),INFO=0, LWORK=-1;

  eigenVectors = mat;
  eigenValues.ReSize(N);

  // T query;
  // utl::syev<T>(&JOBZ,&UPLO,&N,eigenVectors.GetData(), &N, eigenValues.GetData(),&query,&LWORK,&INFO);
  // LWORK=static_cast<int>(query); 

  // T *const WORK = new T[LWORK];
  // utl::syev<T>(&JOBZ,&UPLO,&N,eigenVectors.GetData(), &N, eigenValues.GetData(),WORK,&LWORK,&INFO);
  // delete[] WORK;

  INFO = utl::syev<T>(LAPACK_COL_MAJOR, JOBZ,UPLO,N,eigenVectors.GetData(), N, eigenValues.GetData());
  utlGlobalException(INFO, "LAPACK library function dsyev_() returned error code INFO=" << INFO);
}

/**
 * \brief dsyevd_UtlMatrix 
 *  eigen-decomposition for symmetric matrix. dsyevd is faster than dsyev
 *  http://www.netlib.org/lapack/explore-html/d1/da2/dsyevd_8f.html
 *
 * \param mat  a symmetric matrix (only upper triangular matrix is used)
 * \param eigenValues Eigen-values are in \b increasing order. 
 * \param eigenVectors Eigen-vectors. each row is an eigen-vector
 */
template <class T> inline void 
syevd_UtlMatrix ( const NDArray<T,2>& mat, NDArray<T,1>& eigenValues, NDArray<T,2>& eigenVectors)
{
  utlException(mat.Rows()!=mat.Cols(), "The matrix should be symmetric");
  char JOBZ='V', UPLO='U';
  int N = mat.Rows(),INFO=0, LWORK=-1, LIWORK=-1;

  eigenVectors = mat;
  eigenValues.ReSize(N);

  // T query; 
  // int queryI;
  // utl::syevd<T>(&JOBZ,&UPLO,&N,eigenVectors.GetData(), &N, eigenValues.GetData(),&query,&LWORK,&queryI,&LIWORK,&INFO);
  // LWORK=static_cast<int>(query); 
  // LIWORK=static_cast<int>(queryI); 

  // T *const WORK = new T[LWORK];
  // int *const IWORK = new int[LIWORK];
  // utl::syevd<T>(&JOBZ,&UPLO,&N,eigenVectors.GetData(), &N, eigenValues.GetData(),WORK,&LWORK,IWORK,&LIWORK,&INFO);
  // delete[] WORK;
  // delete[] IWORK;

  INFO = utl::syevd<T>(LAPACK_COL_MAJOR, JOBZ,UPLO,N,eigenVectors.GetData(), N, eigenValues.GetData());
  utlGlobalException(INFO, "LAPACK library function dsyevd_() returned error code INFO=" << INFO);
}


/**
 * \brief dgesvd_UtlMatrix 
 *
 * \param mat matrix with size MxN.
 * \param U left singular vectors. If format is 'A', U is MxM matrix. If format is 'S', U size is M x min(M,N)
 * \param s singular values with size min(M,N). Sored in \b decreasing order.
 * \param V right singular vectors. If format is 'A', V size is NxN. If format is 'S', V size is N x min(M,N)
 * \param format 'S' or 'A'. 'A' means full size, 'S' means reduced size. 
 */
template <class T> inline void 
gesvd_UtlMatrix(const NDArray<T,2>& mat, NDArray<T,2>& U, NDArray<T,1>& s, NDArray<T,2>& V, char format='S' )
{
  utlException(format!='S' && format!='A', "wrong format. format should be 'A' or 'S' ");
  int M=mat.Rows(), N=mat.Cols(), INFO=0;
  int min_MN = utl::min(M,N), LDVT;
  s.ReSize(min_MN);
  if (format=='S')
    {
    U.ReSize(min_MN, M);
    V.ReSize(N, min_MN);
    LDVT = min_MN;
    }
  else
    {
    U.ReSize(M, M);
    V.ReSize(N, N);
    LDVT = N;
    }
  char formatU=format, formatV=format;
  int LDA=M, LDU=M, LWORK=-1;
  NDArray<T,2> matTmp;
  mat.GetTranspose(matTmp);

  // T query;
  // utl::gesvd<T>(&formatU, &formatV, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, &query, &LWORK, &INFO);
  // LWORK=static_cast<int>(query); 

  // T* WORK = new T[LWORK];
  // // dgesdd_(&format, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, WORK, &LWORK, IWORK, &INFO);
  // utl::gesvd<T>(&formatU, &formatV, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, WORK, &LWORK, &INFO);
  // delete[] WORK;

  T* superb = new T[min_MN];
  INFO=utl::gesvd<T>(LAPACK_COL_MAJOR, formatU, formatV, M, N, matTmp.GetData(), LDA, s.GetData(), U.GetData(), LDU, V.GetData(), LDVT, superb);
  delete[] superb;
  U = U.GetTranspose();
  // U.TransposeInplace();

  utlGlobalException(INFO, "LAPACK library function dgesvd_() returned error code INFO=" << INFO);
}

/**
 * \brief dgesdd_UtlMatrix 
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
gesdd_UtlMatrix(const NDArray<T,2>& mat, NDArray<T,2>& U, NDArray<T,1>& s, NDArray<T,2>& V, char format='S' )
{
  utlException(format!='S' && format!='A', "wrong format. format should be 'A' or 'S' ");
  int M=mat.Rows(), N=mat.Cols(), INFO=0;
  int min_MN = utl::min(M,N), LDVT;
  s.ReSize(min_MN);
  if (format=='S')
    {
    U.ReSize(min_MN, M);
    V.ReSize(N, min_MN);
    LDVT = min_MN;
    }
  else
    {
    U.ReSize(M, M);
    V.ReSize(N, N);
    LDVT = N;
    }
  char formatU=format, formatV=format;
  int LDA=M, LDU=M, LWORK;
  NDArray<T,2> matTmp;
  mat.GetTranspose(matTmp);

  // int* IWORK = new int[8*min_MN];
  // T query;
  // LWORK=-1;
  // utl::gesdd<T>(&format, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, &query, &LWORK, IWORK, &INFO);
  // // dgesvd_(&formatU, &formatV, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, &query, &LWORK, &INFO);
  // LWORK=static_cast<int>(query); 

  // // LWORK = min_MN*(6+4*min_MN)+utl::max(M,N);
  // T* WORK = new T[LWORK];
  // utl::gesdd<T>(&format, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, WORK, &LWORK, IWORK, &INFO);
  // // dgesvd_(&formatU, &formatV, &M, &N, matTmp.GetData(), &LDA, s.GetData(), U.GetData(), &LDU, V.GetData(), &LDVT, WORK, &LWORK, &INFO);
  // delete[] WORK;
  // delete[] IWORK;

  INFO=utl::gesdd<T>(LAPACK_COL_MAJOR, format, M, N, matTmp.GetData(), LDA, s.GetData(), U.GetData(), LDU, V.GetData(), LDVT);
  U = U.GetTranspose();
  // U.TransposeInplace();

  utlGlobalException(INFO, "LAPACK library function dgesdd_() returned error code INFO=" << INFO);
}

template <class T> inline void 
InverseMatrix( const NDArray<T,2>& mat, NDArray<T,2>& result, const T eps=1e-8 )
{
  utlSAException(mat.Rows()!=mat.Cols())(mat.Rows())(mat.Cols()).msg("the matrix should be symmetric");
  int n = mat.Rows();
  utlException(n==0, "matrix cannot be empty");
  if (n>=1 && n <=4)
    {
    result.ReSize(n,n); 
    utl::InverseSmallMatrix(mat, result, n);
    }
  else
    {
    utlException(true, "TODO: use LU factorization. See getrf getri CImg");
    }
}

/** inverse of a symmetric matrix (non-singular). If the matrix is singular, it stop with an error. 
 * It is fast than PInverseVnlMatrix, but only works for non-singular matrix. */
template <class T> inline void 
InverseSymmericMatrix( const NDArray<T,2>& mat, NDArray<T,2>& result, const T eps=1e-8 )
{
  utlSAException(mat.Rows()!=mat.Cols())(mat.Rows())(mat.Cols()).msg("the matrix should be symmetric");
  int n = mat.Rows();
  if (n<=4)
    InverseMatrix(mat, result, eps);
  else
    {
    result = mat;
    char uplo='U';
    int lwork=-1, INFO=0;
    int* ipiv= new int[n];

    // T *query, *work;
    // query = new T[1];
    // utl::sytrf<T>(&uplo,&n,result.GetData(),&n,ipiv,query,&lwork,&INFO);
    // lwork=static_cast<INTT>(*query); 
    // delete[] query;
    // work = new T[lwork];
    // utl::sytrf<T>(&uplo,&n,result.GetData(),&n,ipiv,work,&lwork,&INFO);
    // delete[] work;
    // work = new T[2*n];
    // utl::sytri<T>(&uplo,&n,result.GetData(),&n,ipiv,work,&INFO);
    // delete[] work;

    utl::sytrf<T>(LAPACK_COL_MAJOR, uplo,n,result.GetData(),n,ipiv);
    INFO=utl::sytri<T>(LAPACK_COL_MAJOR, uplo,n,result.GetData(),n,ipiv);
    delete[] ipiv;

    utlGlobalException(INFO>0, "The matrix is singular and its inverse could not be computed. \
      LAPACK library function sytrid_() returned error code INFO=" << INFO);
    utlGlobalException(INFO<0, "LAPACK library function sytrid_() returned error code INFO=" << INFO);

    T* data = result.GetData();
    for ( int i = 0; i < n; ++i ) 
      for ( int j = 0; j < i; ++j ) 
        data[j*n+i] = data[i*n+j];
    }
}

/** pseudo-inverse of a symmetric matrix which can be singular. If the matrix is not singular, it returns the inverse.  */
template <class T> inline void 
PInverseSymmericMatrix( const NDArray<T,2>& mat, NDArray<T,2>& result, const T eps=1e-8 )
{
  utlSAException(mat.Rows()!=mat.Cols())(mat.Rows())(mat.Cols()).msg("the matrix should be symmetric");
  int N = mat.Rows();
  NDArray<T,2> eigenVectors, tmp, S(N,N,0.0);
  NDArray<T,1> eigenValues, v;
  mat.EigenDecompositionSymmetricMatrix(eigenValues, eigenVectors);
  for (unsigned i=0; i<N; ++i)
    {
    if ( eigenValues[i]>eps || eigenValues[i]<-eps ) 
      S(i,i) = 1.0/eigenValues[i];
    else
      S(i,i)=0.0;
    }
  utl::ProductUtlMtM(eigenVectors, S, tmp);
  utl::ProductUtlMM(tmp, eigenVectors, result);
}

/** pseudo-inverse of a general matrix.  */
template <class T> inline void 
PInverseMatrix( const NDArray<T,2>& mat, NDArray<T,2>& result, const T eps=1e-8 )
{
  NDArray<T,2> U,V, tmp;
  NDArray<T,1> S, uVec, vVec;
  mat.SVD(U, S, V, 'S');

  NDArray<T,2> diag(S.Size(), S.Size(),0.0);
  for ( int i = 0; i < S.Size(); i += 1 )
    {
    if ( S[i]>eps || S[i]<-eps ) 
      diag(i,i) = 1.0/S[i];
    else
      diag(i,i) = 0.0;
    }
  utl::ProductUtlMM(V, diag, tmp);
  utl::ProductUtlMMt(tmp, U, result);
}


/**
 *   \class   NDArray
 *   \brief   NDArray<T,2> is a row-major matrix
 *   \author  Jian Cheng
 *   \date    08-23-2014
 *   \ingroup utlNDArray
 */
template < class T >
class NDArray<T,2> : public NDArrayBase<T,2>
{
public:
  typedef NDArray                  Self;
  typedef NDArrayBase<T,2>            Superclass;

  typedef typename Superclass::ValueType         ValueType;
  typedef typename Superclass::SizeType          SizeType;
  typedef typename Superclass::ShapeType         ShapeType;
  typedef typename Superclass::Pointer           Pointer;
  typedef typename Superclass::ConstPointer      ConstPointer;
  typedef typename Superclass::Reference         Reference;
  typedef typename Superclass::ConstReference    ConstReference;
  typedef typename Superclass::Iterator          Iterator;
  typedef typename Superclass::ConstIterator     ConstIterator;

  using Superclass::Dimension;
  
  using Superclass::SetData;
  using Superclass::CopyData;
  using Superclass::ReSize;
  using Superclass::operator();
  using Superclass::operator=;
  using Superclass::operator+=;
  using Superclass::operator-=;
  using Superclass::operator%=;
  using Superclass::operator/=;

  NDArray() : Superclass() 
    { }
  explicit NDArray(const SizeType rows, const SizeType columns) : Superclass() 
    { 
    SizeType shape[2];
    shape[0]=rows, shape[1]=columns;
    __utl_ndarray_alloc_blah(shape);
    }
  NDArray(const NDArray<T,2>& mat) : Superclass(mat)
    { }
  template<typename EType>
  NDArray(const Expr<EType>& expr) : Superclass(expr)
    { }
  NDArray(const T* vec, const SizeType rows, const SizeType columns) : Superclass() 
    { 
    SizeType shape[2];
    shape[0]=rows, shape[1]=columns;
    __utl_ndarray_alloc_blah(shape);
    utl::cblas_copy<T>(this->Size(), vec, 1, this->m_Data, 1);
    }
  NDArray(const SizeType rows, const SizeType columns, const T r) : Superclass()
    { 
    SizeType shape[2];
    shape[0]=rows, shape[1]=columns;
    __utl_ndarray_alloc_blah(shape);
    this->Fill(r);
    }

  template< typename TMatrixValueType >
  NDArray(const NDArray< TMatrixValueType,2> & r) : Superclass(r)
    {  }

  explicit NDArray(const ShapeType& shape) : Superclass(shape) { }

/**
 * Constructor assumes input points to array of correct size.
 * Values are copied individually instead of with a binary copy.  This
 * allows the T's assignment operator to be executed.
 */
  NDArray(const T* vec, const ShapeType& shape) : Superclass(vec,shape)   {  }

  /**
   * Constructor to initialize entire array to one value.
   */
  NDArray(const ShapeType& shape, const T r) : Superclass(shape, r)   {  }


  UTL_ALWAYS_INLINE SizeType Rows() const {return this->m_Shape[0];}
  UTL_ALWAYS_INLINE SizeType Columns() const {return this->m_Shape[1];}
  UTL_ALWAYS_INLINE SizeType Cols() const {return this->m_Shape[1];}
  
  inline bool ReSize(const SizeType rows, const SizeType cols)
    {
    SizeType shape[2];
    shape[0]=rows, shape[1]=cols;
    return Superclass::ReSize(shape);
    }
  
  // UTL_ALWAYS_INLINE Reference operator()(const ShapeType& shape)         { return this->m_Data[shape[0]*this->m_Shape[1]+shape[1]]; }
  // UTL_ALWAYS_INLINE ConstReference operator()(const ShapeType& shape) const { return this->m_Data[shape[0]*this->m_Shape[1]+shape[1]]; }

  UTL_ALWAYS_INLINE T & operator()(unsigned int row, unsigned int col)
  { return this->m_Data[row*this->m_Shape[1]+col];}
  UTL_ALWAYS_INLINE const T & operator()(unsigned int row, unsigned int col) const
  { return this->m_Data[row*this->m_Shape[1]+col];}
    // {
    // SizeType index[2];
    // index[0]=row, index[1]=col;
    // return Superclass::operator()(index);
    // }

// #define __Matrix_Saver_Expr(saver)                                 \
// template<typename EType>                                           \
// inline NDArray<T,2>& operator saver (const Expr<EType>& src){         \
//   Superclass::operator saver (src);                                \
//   return *this;                                                    \
// }                                                                  

//   // __Matrix_Saver_Expr(=)
//   __Matrix_Saver_Expr(+=)
//   __Matrix_Saver_Expr(-=)
//   [>* use %= for elementwse *=  <]
//   __Matrix_Saver_Expr(%=)
//   __Matrix_Saver_Expr(/=)

// #define __Matrix_Saver_Scalar(saver)                    \
// inline NDArray<T,2>& operator saver (const T val){         \
//   Superclass::operator saver(val);                      \
//   return *this;                                         \
// }                                                       

//   // __Matrix_Saver_Scalar(=)
//   __Matrix_Saver_Scalar(+=)
//   __Matrix_Saver_Scalar(-=)
//   __Matrix_Saver_Scalar(%=)
//   __Matrix_Saver_Scalar(/=)

  

  /** matrix product  */
  NDArray<T,2>& operator*=(const NDArray<T,2>& mat) 
    {
    utlSAException(this->m_Shape[1]!=mat.m_Shape[0])(this->m_Shape[1])(mat.m_Shape[0]).msg("wrong size");
    NDArray<T,2> tmp;
    utl::ProductUtlMM(*this, mat, tmp);
    Swap(tmp);
    return *this;
    }
  NDArray<T,2> operator*(const NDArray<T,2>& mat) const
    {
    utlSAException(this->m_Shape[1]!=mat.m_Shape[0])(this->m_Shape[1])(mat.m_Shape[0]).msg("wrong size");
    NDArray<T,2> tmp;
    utl::ProductUtlMM(*this, mat, tmp);
    return tmp;
    }
  
  // NDArray<T,2>& operator*=(const NDArray<T,1>& vec) 
  //   {
  //   utlSAException(vec.Size()!=this->m_Shape[1])(this->m_Shape[1])(vec.Size()).msg("wrong size");
  //   NDArray<T,2> tmp;
  //   tmp.SetData((T* const)vec.GetData(), vec.Size(), 1);
  //   operator*=(tmp);
  //   return *this;
  //   }
  // NDArray<T,2> operator*(const NDArray<T,1>& vec) const
  //   {
  //   return NDArray<T,2>(*this)*=vec;
  //   }
  


  /** set diagonal values to val */
  NDArray<T,2>& FillDiagonal(const T val)
    {
    int min_MN = utl::min(this->m_Shape[0], this->m_Shape[1]);
    for ( int i = 0; i < min_MN; ++i ) 
      (*this)(i,i) = val;
    return *this;
    }

  /** set diagonal values from vec  */
  NDArray<T,2>& SetDiagonal(const NDArray<T,1>& vec)
    {
    int min_MN = utl::min(this->m_Shape[0], this->m_Shape[1]);
    utlSAException(vec.Size()!=min_MN)(vec.Size())(min_MN).msg("wrong size of input vector");
    for ( int i = 0; i < min_MN; ++i ) 
      (*this)(i,i) = vec[i];
    return *this;
    }

  void GetDiagonal(NDArray<T,1>& vec) const
    {
    int min_MN = utl::min(this->m_Shape[0], this->m_Shape[1]);
    vec.ReSize(min_MN);
    for ( int i = 0; i < min_MN; ++i ) 
      vec[i] = (*this)(i,i);
    }

  NDArray<T,2> SetIdentity()
    {
    for ( int i = 0; i < this->m_Shape[0]; ++i ) 
      for ( int j = 0; j < this->m_Shape[1]; ++j ) 
        (*this)(i,j) = (i==j) ? (T)1.0 : (T)0.0;
    return *this;
    }
  
  /** set m_Data from data, the ownership is in data. */
  void SetData( T* const data, const unsigned rows, const unsigned cols )
    {
    SizeType shape[2];
    shape[0]=rows, shape[1]=cols;
    Superclass::SetData(data, shape);
    }

  /** copy m_Data from data  */
  void CopyData(T* const data, const unsigned rows, const unsigned cols )
    {
    SizeType shape[2];
    shape[0]=rows, shape[1]=cols;
    Superclass::CopyData(data, shape);
    }
  
  NDArray<T,2> GetTranspose(const T scale=1.0) const
    {
    NDArray<T,2> tmp;
    GetTranspose(tmp, scale);
    return tmp;
    }
  void GetTranspose(NDArray<T,2>& result, const T scale=1.0) const
    {
    result.ReSize(this->m_Shape[1], this->m_Shape[0]);
#ifdef UTL_USE_MKL
    utl::mkl_omatcopy<T>('R', 'T', this->m_Shape[0], this->m_Shape[1], scale, this->m_Data, this->m_Shape[1], result.m_Data, result.m_Shape[1]); 
#else
    for ( int i = 0; i < this->m_Shape[0]; ++i ) 
      for ( int j = 0; j < this->m_Shape[1]; ++j ) 
        result[j*this->m_Shape[0]+i] = (*this)[i*this->m_Shape[1]+j];
    result.Scale(scale);
#endif
    }

  NDArray<T,2>& TransposeInplace(const T scale=1.0)
    {
#ifdef UTL_USE_MKL
    utl::mkl_imatcopy<T>('R', 'T', this->m_Shape[0], this->m_Shape[1], scale, this->m_Data, this->m_Shape[1], this->m_Shape[0]); 
    std::swap(this->m_Shape[0], this->m_Shape[1]);
    this->ComputeOffSetTable();
#else
    // NOTE: force swap row and col, because when -O3 is used, it may not swap the shape. 
    SizeType row=this->m_Shape[0], col=this->m_Shape[1];
    this->GetTranspose(scale).Swap(*this);
    this->m_Shape[0]=col, this->m_Shape[1]=row;
#endif
    return *this;
    }

  NDArray<T,2>& SetRow(const int index, T const* vec)
    {
    utlSAException(index<0 || index>=this->m_Shape[0])(index)(this->m_Shape[0]).msg("wrong index");
    utl::cblas_copy<T>(this->m_Shape[1],vec,1,this->m_Data+this->m_Shape[1]*index,1);
    return *this;
    }
  NDArray<T,2>& SetRow(const int index, const NDArray<T,1>& vec)
    { 
    utlSAException(vec.Size()!=this->m_Shape[1])(vec.Size())(this->m_Shape[1]).msg("wrong size of vec");
    SetRow(index, vec.GetData());
    return *this;
    }
  NDArray<T,2>& SetColumn(const int index, T const* vec)
    {
    utlSAException(index<0 || index>=this->m_Shape[1])(index)(this->m_Shape[1]).msg("wrong index");
    utl::cblas_copy<T>(this->m_Shape[0],vec,1,this->m_Data+index,this->m_Shape[1]);
    return *this;
    }
  NDArray<T,2>& SetColumn(const int index, const NDArray<T,1>& vec)
    { 
    utlSAException(vec.Size()!=this->m_Shape[0])(vec.Size())(this->m_Shape[0]).msg("wrong size of vec");
    SetColumn(index, vec.GetData());
    return *this;
    }

  void GetRow(const int index, T* vec) const
    {
    utlSAException(index<0 || index>=this->m_Shape[0])(index)(this->m_Shape[0]).msg("wrong index");
    utl::cblas_copy<T>(this->m_Shape[1],this->m_Data+this->m_Shape[1]*index,1,vec,1);
    }
  void GetRow(const int index, NDArray<T,1>& vec) const
    {
    vec.ReSize(this->m_Shape[1]);
    GetRow(index, vec.GetData());
    }
  NDArray<T,1> GetRow(const int index) const
    {
    NDArray<T,1> vec(this->m_Shape[1]);
    GetRow(index, vec.GetData());
    return vec;
    }
  void GetColumn(const int index, T* vec) const
    {
    utlSAException(index<0 || index>=this->m_Shape[1])(index)(this->m_Shape[1]).msg("wrong index");
    utl::cblas_copy<T>(this->m_Shape[0],this->m_Data+index,this->m_Shape[1],vec,1);
    }
  void GetColumn(const int index, NDArray<T,1>& vec) const
    {
    vec.ReSize(this->m_Shape[0]);
    GetColumn(index, vec.GetData());
    }
  NDArray<T,1> GetColumn(const int index) const
    {
    NDArray<T,1> vec(this->m_Shape[0]);
    GetColumn(index, vec.GetData());
    return vec;
    }

  void GetRows(const std::vector<int>& indexVec, NDArray<T,2>& mat) const 
    {
    mat.ReSize(indexVec.size(), this->m_Shape[1]);
    for ( int i = 0; i < indexVec.size(); ++i ) 
      GetRow(indexVec[i], mat.m_Data+i*this->m_Shape[1]);
    }
  void GetNRows(const int index0, const int N, NDArray<T,2>& mat) const 
    {
    std::vector<int> indexVec;
    for ( int i = 0; i < N; ++i ) 
      indexVec.push_back(index0+i);
    GetRows(indexVec, mat);
    }
  void GetColumns(const std::vector<int>& indexVec, NDArray<T,2>& mat) const 
    {
    mat.ReSize(this->m_Shape[0], indexVec.size());
    NDArray<T,1> vec;
    for ( int i = 0; i < indexVec.size(); ++i ) 
      {
      GetColumn(indexVec[i], vec);
      mat.SetColumn(i, vec);
      }
    }
  void GetNColumns(const int index0, const int N, NDArray<T,2>& mat) const 
    {
    std::vector<int> indexVec;
    for ( int i = 0; i < N; ++i ) 
      indexVec.push_back(index0+i);
    GetColumns(indexVec, mat);
    }
  NDArray<T,2>& SetRows(const std::vector<int>& indexVec, const NDArray<T,2>& mat)
    {
    NDArray<T,1> vec;
    for ( int i = 0; i < indexVec.size(); ++i ) 
      {
      mat.GetRow(i,vec);
      SetRow(indexVec[i], vec);
      }
    return *this;
    }
  NDArray<T,2>& SetColumns(const std::vector<int>& indexVec, const NDArray<T,2>& mat)
    {
    NDArray<T,1> vec;
    for ( int i = 0; i < indexVec.size(); ++i ) 
      {
      mat.GetColumn(i,vec);
      SetColumn(indexVec[i], vec);
      }
    return *this;
    }

  /** get a submatrix 
   *    \param x0 starting row index
   *    \param x1 ending row index
   *    \param y0 starting column index
   *    \param y1 ending column index
   * */
  void GetCrop(const int x0, const int x1, const int y0, const int y1, NDArray<T,2>& mat) const
    {
    utlSAException(x0<0 || x0>=this->m_Shape[0] || x1<0 || x1>=this->m_Shape[0] || x0>x1)(x0)(x1)(this->m_Shape[0]).msg("wrong index");
    utlSAException(y0<0 || y0>=this->m_Shape[0] || y1<0 || y1>=this->m_Shape[1] || y0>y1)(y0)(y1)(this->m_Shape[1]).msg("wrong index");
    mat.ReSize(x1-x0+1, y1-y0+1);
    for ( int i = 0; i < mat.Rows(); ++i ) 
      for ( int j = 0; j < mat.Cols(); ++j ) 
        mat(i,j) = (*this)(i+x0, j+x1);
    }
  
  NDArray<T,2>& SetCrop(const int x0, const int x1, const int y0, const int y1, const NDArray<T,2>& mat) 
    {
    utlSAException(x0<0 || x0>=this->m_Shape[0] || x1<0 || x1>=this->m_Shape[0] || x0>x1)(x0)(x1)(this->m_Shape[0]).msg("wrong index");
    utlSAException(y0<0 || y0>=this->m_Shape[0] || y1<0 || y1>=this->m_Shape[1] || y0>y1)(y0)(y1)(this->m_Shape[1]).msg("wrong index");
    utlSAException(mat.Rows()!=x1-x0+1 || mat.Cols()!=y1-y0+1)(mat.Rows())(mat.Cols()).msg("wrong size");
    for ( int i = 0; i < mat.Rows(); ++i ) 
      for ( int j = 0; j < mat.Cols(); ++j ) 
        (*this)(i+x0, j+x1) = mat(i,j);
    return *this;
    }

  /**
   * \brief SVD
   *
   * \param U left singular vectors. If format is 'A', U is MxM matrix. If format is 'S', U size is M x min(M,N)
   * \param s singular values with size min(M,N). Sored in \b decreasing order.
   * \param V right singular vectors. If format is 'A', V size is NxN. If format is 'S', V size is N x min(M,N)
   * \param format 'S' or 'A'. 'A' means full size, 'S' means reduced size. 
   */
  void SVD(NDArray<T,2>& U, NDArray<T,1>& S, NDArray<T,2>& V, char format='S') const
    {
#ifdef UTL_USE_FASTLAPACK
  // NOTE: fast but may be wrong for big matrix and multi-thread if openblas is not correctly built
  utl::gesdd_UtlMatrix<T>(*this, U, S, V, format);
#else
  // slow but robust
  utl::gesvd_UtlMatrix<T>(*this, U, S, V, format);
#endif
    }


  void EigenDecompositionSymmetricMatrix (NDArray<T,1>& eigenValues, NDArray<T,2>& eigenVectors ) const
    {
    utlSAException(this->m_Shape[0]!=this->m_Shape[1])(this->m_Shape[0])(this->m_Shape[1]).msg("the matrix should be symmetric");
#ifdef UTL_USE_FASTLAPACK
  // NOTE: fast but may be wrong for big matrix and multi-thread if openblas is not correctly built
  syevd_UtlMatrix<T>(*this, eigenValues, eigenVectors);
#else
  // slow but robust
  syev_UtlMatrix<T>(*this, eigenValues, eigenVectors);
#endif
    }
  
  void InverseSymmericMatrix(NDArray<T,2>& result, const T eps=1e-8)
    {
    utl::InverseSymmericMatrix(*this, result, eps);
    }

  void PInverseSymmericMatrix(NDArray<T,2>& result, const T eps=1e-8)
    {
    utl::PInverseSymmericMatrix(*this, result, eps);
    }
  
  void PInverseMatrix(NDArray<T,2>& result, const T eps=1e-8)
    {
    utl::PInverseMatrix(*this, result, eps);
    }
  void InverseMatrix(NDArray<T,2>& result, const T eps=1e-8)
    {
    utl::InverseMatrix(*this, result, eps);
    }

  double Determinant() const
    {
    utlException(!IsSquareMatrix(), "should be a square matrix");
    utlException(this->m_Shape[0]==0, "should not be empty");
    const T* p = this->m_Data;
    switch ( this->m_Shape[0] )
      {
      case 1 :  return p[0];
      case 2 :  return p[0]* p[3] - p[1]*p[2];
      case 3 :
           {
          const double
            a = p[0], d = p[1], g = p[2],
            b = p[3], e = p[4], h = p[5],
            c = p[6], f = p[7], i = p[8];
          return i*a*e-a*h*f-i*b*d+b*g*f+c*d*h-c*g*e;
           }
      case 4 :
          return
            + p[0]*p[5]*p[10]*p[15]
            - p[0]*p[5]*p[14]*p[11]
            - p[0]*p[9]*p[6]*p[15]
            + p[0]*p[9]*p[14]*p[7]
            + p[0]*p[13]*p[6]*p[11]
            - p[0]*p[13]*p[10]*p[7]
            - p[4]*p[1]*p[10]*p[15]
            + p[4]*p[1]*p[14]*p[11]
            + p[4]*p[9]*p[2]*p[15]
            - p[4]*p[9]*p[14]*p[3]
            - p[4]*p[13]*p[2]*p[11]
            + p[4]*p[13]*p[10]*p[3]
            + p[8]*p[1]*p[6]*p[15]
            - p[8]*p[1]*p[14]*p[7]
            - p[8]*p[5]*p[2]*p[15]
            + p[8]*p[5]*p[14]*p[3]
            + p[8]*p[13]*p[2]*p[7]
            - p[8]*p[13]*p[6]*p[3]
            - p[12]*p[1]*p[6]*p[11]
            + p[12]*p[1]*p[10]*p[7]
            + p[12]*p[5]*p[2]*p[11]
            - p[12]*p[5]*p[10]*p[3]
            - p[12]*p[9]*p[2]*p[7]
            + p[12]*p[9]*p[6]*p[3];
      default :
           {
           utlException(true, "TODO use LU factorization for matrix");
           return 0;
           }
      }
    }

  bool IsSquareMatrix() const
    {
    return this->m_Shape[0] == this->m_Shape[1];
    }

  bool IsSymmetric(const double eps=1e-10) const
    {
    if (this->m_Shape[0]!=this->m_Shape[1])
      return false;
    for ( int i = 0; i < this->m_Shape[0]; ++i ) 
      for ( int j = 0; j < i; ++j ) 
        {
        double v1 = (*this)(i,j), v2 = (*this)(j,i);
        if (std::fabs(v1-v2)>eps*std::fabs(v1))
          return false;
        }
    return true;
    }

  T GetArrayOneNorm() const { return Superclass::GetOneNorm(); }
  T GetArrayTwoNorm() const { return Superclass::GetTwoNorm(); }
  T GetArrayInfNorm() const { return Superclass::GetInfNorm(); }
  int GetArrayZeroNorm() const { return Superclass::GetZeroNorm(); }

  T GetTwoNorm() const {return GetArrayTwoNorm();}
  T GetOneNorm() const
    { return utl::lange<T>(LAPACK_ROW_MAJOR, 'O', this->m_Shape[0], this->m_Shape[1], this->m_Data, this->m_Shape[1]); }
  T GetInfNorm() const
    { return utl::lange<T>(LAPACK_COL_MAJOR, 'O', this->m_Shape[1], this->m_Shape[0], this->m_Data, this->m_Shape[1]); }

  void Save(const std::string file) const
    {
    utl::SaveMatrix<Self>(*this, this->m_Shape[0], this->m_Shape[1], file);
    }

protected:

private:

  // not used for matrix
  // Reference operator()(unsigned int index) ;
  // ConstReference operator()(unsigned int index) const ;

}; // -----  end of template class Matrix  -----


template< typename T >
std::ostream & 
operator<<(std::ostream & os, const NDArray< T,2 > & arr)
{
  utl::PrintUtlMatrix(arr, "utl::NDArray<T,2>", " ", os);
  return os;
}

template <class T>
NDArray<T,2>
ConnectUtlMatrix ( const NDArray<T,2>& m1, const NDArray<T,2>& m2, const bool isConnectRow )
{
  if (isConnectRow)
    {
    utlException(m1.Columns()!=m2.Columns(), "wrong column size! m1.Columns()="<<m1.Columns()<<", m2.Columns()"<<m2.Columns());
    NDArray<T,2> result(m1.Rows()+m2.Rows(), m1.Columns());
    int m1Rows = m1.Rows();
    for ( int j = 0; j < m1.Columns(); j += 1 ) 
      {
      for ( int i = 0; i < m1.Rows(); i += 1 ) 
        result(i,j) = m1(i,j);
      for ( int i = 0; i < m2.Rows(); i += 1 ) 
        result(i+m1Rows,j) = m2(i,j);
      }
    return result;
    }
  else
    {
    utlException(m1.Rows()!=m2.Rows(), "wrong row size! m1.Rows()="<<m1.Rows()<<", m2.Rows()"<<m2.Rows());
    NDArray<T,2> result(m1.Rows(), m1.Columns()+m2.Columns());
    int m1Columns = m1.Columns();
    for ( int i = 0; i < m1.Rows(); i += 1 ) 
      {
      for ( int j = 0; j < m1.Columns(); j += 1 ) 
        result(i,j) = m1(i,j);
      for ( int j = 0; j < m2.Columns(); j += 1 ) 
        result(i,j+m1Columns) = m2(i,j);
      }
    return result;
    }
}

template <class T>
NDArray<T,2> 
SphericalToCartesian ( const NDArray<T,2>& in )
{
  utlAssert(in.Columns()==3 || in.Rows()==3, "wrong dimension");
  NDArray<T,2> out(in);
  if (in.Columns()==3)
    {
    for ( int i = 0; i < in.Rows(); i += 1 ) 
      utl::spherical2Cartesian(in(i,0), in(i,1), in(i,2), out(i,0), out(i,1), out(i,2));
    }
  else
    {
    for ( int i = 0; i < in.Columns(); i += 1 ) 
      utl::spherical2Cartesian(in(0,i), in(1,i), in(2,i), out(0,i), out(1,i), out(2,i));
    }
  return out;
}

template <class T>
NDArray<T,2> 
CartesianToSpherical ( const NDArray<T,2>& in )
{
  utlAssert(in.Columns()==3 || in.Rows()==3, "wrong dimension");
  NDArray<T,2> out(in);
  if (in.Columns()==3)
    {
    for ( int i = 0; i < in.Rows(); i += 1 ) 
      utl::cartesian2Spherical(in(i,0), in(i,1), in(i,2), out(i,0), out(i,1), out(i,2));
    }
  else
    {
    for ( int i = 0; i < in.Columns(); i += 1 ) 
      utl::cartesian2Spherical(in(0,i), in(1,i), in(2,i), out(0,i), out(1,i), out(2,i));
    }
  return out;
}

template <class T>
void 
GetEqualityConstraintProjection ( const NDArray<T,2>& Aeq, const NDArray<T,1>& beq, const NDArray<T,2>& QInverse, 
  NDArray<T,2>& projMatrix, NDArray<T,1>& projVector )
{
  int n = QInverse.Rows();
  utlException(Aeq.Rows()!=n, "wrong size! Aeq.rows()="<<Aeq.Rows()<<", n="<<n);
  NDArray<T,2> AeqT=Aeq.GetTranspose();
  projMatrix.ReSize(n,n);
  projMatrix.SetIdentity();
  NDArray<T,2> tmp, tmp2, tmp3;
  ProductUtlMM(AeqT, QInverse, Aeq, tmp2);
  PInverseSymmericMatrix(tmp2, tmp);
  ProductUtlMM( QInverse, Aeq, tmp, tmp2);
  ProductUtlMM( tmp2, AeqT, tmp3);
  projMatrix -= tmp3;

  ProductUtlMv( tmp2, beq, projVector);
}

/** @}  */

}

#endif 
