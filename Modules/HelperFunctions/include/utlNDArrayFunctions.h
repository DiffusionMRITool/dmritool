/**
 *       @file  utlNDArrayFunctions.h
 *      @brief  
 *     Created  "07-03-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlNDArrayFunctions_h
#define __utlNDArrayFunctions_h


#include "utlCoreMacro.h"
#include "utlSTDHeaders.h"
#include "utlExpression.h"
#include "utlCore.h"
#include "utlBlas.h"
#include "utlLapack.h"
#include "utlMath.h"
#include "utlTypeinfo.h"

// #include "utlNDArray.h"

/** \ingroup utlNDArray
* @{ */

namespace utl
{


template < class T, unsigned int Dim >
class NDArray;

template < class T, unsigned int Dim >
class NDArrayBase;

/** return true only if two arrays have the same dimension and same shape.  */
template <class T1, class T2, unsigned Dim1, unsigned Dim2>
inline bool 
IsSameShape( const NDArray<T1,Dim1>& arr1, const NDArray<T2,Dim2>& arr2 )
{
  if (Dim1!=Dim2)
    return false;
  else 
    return arr1.IsSameShape(arr2);
}

/**************************************************************************************/
/*  Conversion   */
/**************************************************************************************/

template <class T, unsigned int Dim, class EType>
NDArray<T,Dim> 
Eval ( const Expr<EType, typename EType::ValueType>& expr )
{
  NDArray<T,Dim> mat(expr);
  return mat;
}

/**
 *  Convert utl::Expression to utl::NDArray 
 */
template <class T, unsigned int Dim, class EType>
utl_shared_ptr< NDArray<T,Dim> >
ToNDArray ( const Expr<EType, typename EType::ValueType>& expr )
{
  utl_shared_ptr< NDArray<T,Dim> > mat(new NDArray<T,Dim>(expr) );
  return mat;
}

/** Convert utl::Expression to utl::Vector   */
template <class T, class EType>
utl_shared_ptr< NDArray<T,1> >
ToVector ( const Expr<EType, typename EType::ValueType>& expr )
{
  utl_shared_ptr< NDArray<T,1> > vec (new NDArray<T,1>(expr));
  return vec;
}

/** Convert utl::Expression to utl::Matrix   */
template <class T, class EType>
utl_shared_ptr< NDArray<T,2> >
ToMatrix ( const Expr<EType, typename EType::ValueType>& expr )
{
  utl_shared_ptr< NDArray<T,2> > mat(new NDArray<T,2>(expr));
  return mat;
}

template <class T>
std::vector<T> 
UtlVectorToStdVector ( const NDArray<T,1>& vec )
{
  std::vector<T> v(vec.Size());
  for ( int i = 0; i < vec.Size(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
NDArray<T,1> 
StdVectorToUtlVector ( const std::vector<T>& vec )
{
  utl::NDArray<T,1> v(vec.size());
  for ( int i = 0; i < vec.size(); ++i ) 
    v[i] = vec[i];
  return v;
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

/** The function only works for 4th order tensors which are in 3D space and have minor symmetry.  
 *  
 *  When \f$ A_{i,j,k,l} \f$ is minor symmetric, it can be mapped to a 2D  matrix 
 *  
 *  \f[ 
 *  \begin{bmatrix} 
 *  A_{1111} & A_{1122} & A_{1133} &  \sqrt{2} A_{1112} & \sqrt{2} A_{1113} & \sqrt{2} A_{1123} \\ 
 *  A_{2211} & A_{2222} & A_{2233} &  \sqrt{2} A_{2212} & \sqrt{2} A_{2213} & \sqrt{2} A_{2223} \\ 
 *  A_{3311} & A_{3322} & A_{3333} &  \sqrt{2} A_{3312} & \sqrt{2} A_{3313} & \sqrt{2} A_{3323} \\ 
 *  \sqrt{2} A_{1211} & \sqrt{2} A_{1222} &  \sqrt{2} A_{1233} &  2 A_{1212} & 2 A_{1213} & 2 A_{1223} \\ 
 *  \sqrt{2} A_{1311} & \sqrt{2} A_{1322} &  \sqrt{2} A_{1333} &  2 A_{1312} & 2 A_{1313} & 2 A_{1323} \\ 
 *  \sqrt{2} A_{2311} & \sqrt{2} A_{2322} &  \sqrt{2} A_{2333} &  2 A_{2312} & 2 A_{2313} & 2 A_{2323} 
 *  \end{bmatrix} 
 *  \f]
 *
 * */
template <class T>
inline void
Convert4To2Tensor(const utl::NDArray<T,4>& tensor, utl::NDArray<T,2>& mat)
{
  unsigned int const* shape = tensor.GetShape();

  utlException(shape[0]!=3 || shape[1]!=3 || shape[2]!=3 || shape[3]!=3, "It should be in 3D space");

  constexpr double s2 = utl::SQRT2;
  mat.ReSize(6,6);

  mat(0,0)=tensor(0,0,0,0), mat(0,1)=tensor(0,0,1,1), mat(0,2)=tensor(0,0,2,2), mat(0,3)=s2*tensor(0,0,0,1), mat(0,4)=s2*tensor(0,0,0,2), mat(0,5)=s2*tensor(0,0,1,2);
  mat(1,0)=tensor(1,1,0,0), mat(1,1)=tensor(1,1,1,1), mat(1,2)=tensor(1,1,2,2), mat(1,3)=s2*tensor(1,1,0,1), mat(1,4)=s2*tensor(1,1,0,2), mat(1,5)=s2*tensor(1,1,1,2);
  mat(2,0)=tensor(2,2,0,0), mat(2,1)=tensor(2,2,1,1), mat(2,2)=tensor(2,2,2,2), mat(2,3)=s2*tensor(2,2,0,1), mat(2,4)=s2*tensor(2,2,0,2), mat(2,5)=s2*tensor(2,2,1,2);
  mat(3,0)=s2*tensor(0,1,0,0), mat(3,1)=s2*tensor(0,1,1,1), mat(3,2)=s2*tensor(0,1,2,2), mat(3,3)=2.0*tensor(0,1,0,1), mat(3,4)=2.0*tensor(0,1,0,2), mat(3,5)=2.0*tensor(0,1,1,2);
  mat(4,0)=s2*tensor(0,2,0,0), mat(4,1)=s2*tensor(0,2,1,1), mat(4,2)=s2*tensor(0,2,2,2), mat(4,3)=2.0*tensor(0,2,0,1), mat(4,4)=2.0*tensor(0,2,0,2), mat(4,5)=2.0*tensor(0,2,1,2);
  mat(5,0)=s2*tensor(1,2,0,0), mat(5,1)=s2*tensor(1,2,1,1), mat(5,2)=s2*tensor(1,2,2,2), mat(5,3)=2.0*tensor(1,2,0,1), mat(5,4)=2.0*tensor(1,2,0,2), mat(5,5)=2.0*tensor(1,2,1,2);
}

template <class T>
inline void
Convert2To4Tensor(const utl::NDArray<T,2>& mat, utl::NDArray<T,4>& tensor)
{
  utlException(mat.Rows()!=6 || mat.Cols()!=6, "Mat should be in 6x6 matrix");
  tensor.ReSize(3,3,3,3);

  constexpr double si2 = utl::SQRT1_2;

  tensor(0,0,0,0)=mat(0,0);
  tensor(0,0,0,1)=si2*mat(0,3); tensor(0,0,1,0)=tensor(0,0,0,1);
  tensor(0,0,0,2)=si2*mat(0,4); tensor(0,0,2,0)=tensor(0,0,0,2);
  tensor(0,0,1,1)=mat(0,1); 
  tensor(0,0,1,2)=si2*mat(0,5); tensor(0,0,2,1)=tensor(0,0,1,2); 
  tensor(0,0,2,2)=mat(0,2); 
  tensor(0,1,0,0)=si2*mat(3,0); tensor(1,0,0,0)=tensor(0,1,0,0);
  tensor(0,1,0,1)=0.5*mat(3,3); tensor(1,0,0,1)=tensor(0,1,0,1); tensor(0,1,1,0)=tensor(0,1,0,1); tensor(1,0,1,0)=tensor(0,1,0,1);
  tensor(0,1,0,2)=0.5*mat(3,4); tensor(0,1,2,0)=tensor(0,1,0,2); tensor(1,0,0,2)=tensor(0,1,0,2); tensor(1,0,2,0)=tensor(0,1,0,2);
  tensor(0,1,1,1)=si2*mat(3,1); tensor(1,0,1,1)=tensor(0,1,1,1);
  tensor(0,1,1,2)=0.5*mat(3,5); tensor(1,0,1,2)=tensor(0,1,1,2); tensor(0,1,2,1)=tensor(0,1,1,2); tensor(1,0,2,1)=tensor(0,1,1,2);
  tensor(0,1,2,2)=si2*mat(3,2); tensor(1,0,2,2)=tensor(0,1,2,2);
  tensor(0,2,0,0)=si2*mat(4,0); tensor(2,0,0,0)=tensor(0,2,0,0);
  tensor(0,2,0,1)=0.5*mat(4,3); tensor(0,2,1,0)=tensor(0,2,0,1); tensor(2,0,0,1)=tensor(0,2,0,1); tensor(2,0,1,0)=tensor(0,2,0,1);
  tensor(0,2,0,2)=0.5*mat(4,4); tensor(2,0,0,2)=tensor(0,2,0,2); tensor(0,2,2,0)=tensor(0,2,0,2); tensor(2,0,2,0)=tensor(0,2,0,2);
  tensor(0,2,1,1)=si2*mat(4,1); tensor(2,0,1,1)=tensor(0,2,1,1);
  tensor(0,2,1,2)=0.5*mat(4,5); tensor(2,0,1,2)=tensor(0,2,1,2); tensor(0,2,2,1)=tensor(0,2,1,2); tensor(2,0,2,1)=tensor(0,2,1,2);
  tensor(0,2,2,2)=si2*mat(4,2); tensor(2,0,2,2)=tensor(0,2,2,2);
  tensor(1,1,0,0)=mat(1,0);
  tensor(1,1,0,1)=si2*mat(1,3); tensor(1,1,1,0)=tensor(1,1,0,1);
  tensor(1,1,0,2)=si2*mat(1,4); tensor(1,1,2,0)=tensor(1,1,0,2);
  tensor(1,1,1,1)=mat(1,1);
  tensor(1,1,1,2)=si2*mat(1,5); tensor(1,1,2,1)=tensor(1,1,1,2);
  tensor(1,1,2,2)=mat(1,2);
  tensor(1,2,0,0)=si2*mat(5,0); tensor(2,1,0,0)=tensor(1,2,0,0);
  tensor(1,2,0,1)=0.5*mat(5,3); tensor(2,1,0,1)=tensor(1,2,0,1); tensor(1,2,1,0)=tensor(1,2,0,1); tensor(2,1,1,0)=tensor(1,2,0,1);
  tensor(1,2,0,2)=0.5*mat(5,4); tensor(2,1,0,2)=tensor(1,2,0,2); tensor(1,2,2,0)=tensor(1,2,0,2); tensor(2,1,2,0)=tensor(1,2,0,2);
  tensor(1,2,1,1)=si2*mat(5,1); tensor(2,1,1,1)=tensor(1,2,1,1);
  tensor(1,2,1,2)=0.5*mat(5,5); tensor(2,1,1,2)=tensor(1,2,1,2); tensor(1,2,2,1)=tensor(1,2,1,2); tensor(2,1,2,1)=tensor(1,2,1,2);
  tensor(1,2,2,2)=si2*mat(5,2); tensor(2,1,2,2)=tensor(1,2,2,2);
  tensor(2,2,0,0)=mat(2,0);
  tensor(2,2,0,1)=si2*mat(2,3); tensor(2,2,1,0)=tensor(2,2,0,1);
  tensor(2,2,0,2)=si2*mat(2,4); tensor(2,2,2,0)=tensor(2,2,0,2);
  tensor(2,2,1,1)=mat(2,1);
  tensor(2,2,1,2)=si2*mat(2,5); tensor(2,2,2,1)=tensor(2,2,1,2);
  tensor(2,2,2,2)=mat(2,2);
}


/**************************************************************************************/
/* Print  */
/**************************************************************************************/

template <class T>
inline void
PrintUtlMatrix(const NDArray<T,2>& mat, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  utl::PrintMatrix<NDArray<T,2> >(mat, mat.Rows(), mat.Columns(), str, separate, os);
}

template <class T>
inline void
PrintUtlVector(const NDArray<T,1>& vec, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  PrintContainer(vec.Begin(), vec.End(), str, separate, os);
}

template <class T, unsigned int Dim>
inline void
PrintUtlNDArray(const NDArrayBase<T,Dim>& arr, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  arr.Print(os<< str, separate );
}

template <class T=double>
utl::NDArray<T,2>
Eye(const int n, const T val=1.0)
{
  utl::NDArray<T,2> mat(n,n);
  mat.Fill(0.0);
  mat.FillDiagonal(val);
  return mat;
}

template <class T=double>
utl::NDArray<T,1>
Ones(const int n)
{
  utl::NDArray<T,1> mat(n);
  mat.Fill(1.0);
  return mat;
}

template <class T=double>
utl::NDArray<T,2>
Ones(const int n, const int m)
{
  utl::NDArray<T,2> mat(n, m);
  mat.Fill(1.0);
  return mat;
}

template< typename T, unsigned int Dim >
std::ostream & 
operator<<(std::ostream & os, const NDArray<T,Dim> & arr)
{
  arr.Print(os<< "utl::NDArray<T,Dim>" );
  return os;
}

template< typename T >
std::ostream & 
operator<<(std::ostream & os, const NDArray<T,1> & arr)
{
  utl::PrintContainer(arr.Begin(), arr.End(), "utl::NDArray<T,1>", " ", os);
  return os;
}

template< typename T >
std::ostream & 
operator<<(std::ostream & os, const NDArray< T,2 > & arr)
{
  utl::PrintUtlMatrix(arr, "utl::NDArray<T,2>", " ", os);
  return os;
}

/**************************************************************************************/
/* Blas and Lapack functions  */
/**************************************************************************************/

/** http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
*  https://developer.apple.com/library/mac/documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html 
*
*  \f$ Y = alpha * A * X + beta * Y \f$
*  \note: Y should be pre-allocated
*
*  define several functions. 
*
*  template \<class T\>  inline bool  
*  gemv_UtlMatrixTimesVector(const bool bATrans, const T alpha, const utl::NDArray<T,2>& A, const utl::NDArray<T,1>& X, const T beta, utl::NDArray<T,1>& Y);
*
*  \f$ Y = alpha * A * X + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlMv(const utl::NDArray<T,2>& A, const utl::NDArray<T,1>& b, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
*  \f$ Y = alpha * A^T * X + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlMtv(const utl::NDArray<T,2>& A, const utl::NDArray<T,1>& b, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
* */
__utl_gemv_MatrixTimesVector(T, gemv_UtlMatrixTimesVector, Utl, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, utl::NDArray<T UTL_COMMA 1>, Size, GetData, ReSize);


/**
*  \f$ Y = alpha X^T * A + beta * Y \f$
*  \note: Y should be pre-allocated
*
*  define several functions. 
*
*  template \<class T\>  inline bool  
*  gemm_UtlVectorTimesMatrix(const bool bATrans, const T alpha, const utl::NDArray<T,1>& X, const utl::NDArray<T,2>& A, const T beta, utl::NDArray<T,1>& Y);
*
*  \f$ Y = alpha * b * A + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlvM(const utl::NDArray<T,1>& b, const utl::NDArray<T,2>& A, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
*  \f$ Y = alpha * b * A^T + beta * Y \f$ \n
*  template \<class T\>  inline void  
*  ProductUtlvMt(const utl::NDArray<T,1>& b, const utl::NDArray<T,2>& A, utl::NDArray<T,1>& c, const double alpha=1.0, const double beta=0.0);
*
* */
__utl_gevm_MatrixTimesVector(T, gevm_UtlVectorTimesMatrix, Utl, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, utl::NDArray<T UTL_COMMA 1>, Size, GetData, ReSize);
 
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
 *
 * \brief geev_UtlMatrix 
 *  calculate non-symmetric eigen-decomposition. 
 *
 * \param mat matrix with size NxN.
 * \param valReal real part of right eigen-values.
 * \param valImg imginary part of right eigen-values.
 * \param vecRealR real part of right eigen-vectors.
 * \param vecImgR part of right eigen-vectors.
 * \param vecRealL real part of left eigen-vectors.
 * \param vecImgL part of left eigen-vectors.
 *
 * template \<class T\> inline void 
 * geev_UtlMatrix ( const NDArray<T,2>& mat, NDArray<T,1>& valReal, NDArray<T,1>& valImg);
 *
 * template \<class T\> inline void 
 * geev_UtlMatrix ( const NDArray<T,2>& mat, NDArray<T,1>& valReal, NDArray<T,1>& valImg, NDArray<T,2>& vecRealR, NDArray<T,2>& vecImgR);
 * 
 * template \<class T\> inline void 
 * geev_UtlMatrix ( const NDArray<T,2>& mat, NDArray<T,1>& valReal, NDArray<T,1>& valImg, NDArray<T,2>& vecRealR, NDArray<T,2>& vecImgR, NDArray<T,2>& vecRealL, NDArray<T,2>& vecImgL);
 *
 * http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga8ec1625302675b981eb34ed024b27a47.html 
 * http://www.netlib.org/lapack/lug/node31.html
 *
 * */
__utl_geev_Matrix(T, geev_UtlMatrix, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, utl::NDArray<T UTL_COMMA 1>, GetData, ReSize);

/** 
 *  
 * \brief geev_UtlMatrix 
 *  calculate inverse of a general matrix via LU decomposition. 
 *
 * template \<class T\> inline void 
 * getri_UtlMatrix(const NDArray<T,2>& mat, NDArray<T,2>& result);
 *
 * */
__utl_getri_Matrix(T, getri_UtlMatrix, utl::NDArray<T UTL_COMMA 2>, Rows, Cols, GetData, ReSize); 

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
  matTmp = mat.GetTranspose();

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
gesdd_UtlMatrix(const NDArray<T,2>& mat, NDArray<T,2>& U, NDArray<utl::remove_complex_t<T>,1>& s, NDArray<T,2>& V, char format='S' )
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
  matTmp = mat.GetTranspose();

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

  INFO=utl::gesdd(LAPACK_COL_MAJOR, format, M, N, matTmp.GetData(), LDA, s.GetData(), U.GetData(), LDU, V.GetData(), LDVT);
  U = U.GetTranspose();
  // U.TransposeInplace();

  utlGlobalException(INFO, "LAPACK library function dgesdd_() returned error code INFO=" << INFO);
}

/**************************************************************************************/
/* Products  */
/**************************************************************************************/

/** 
 *  If value type is std::complex<double>, then utl::InnerProduct provides \f$ \sum_i \text{Conj}(v1_i) v2_i \f$. 
 *  For real values, the result is the same as utl::DotProduct
 * */
template <class T, unsigned int Dim>
inline T
InnerProduct ( const NDArrayBase<T,Dim>& v1, const NDArrayBase<T,Dim>& v2 )
{
  utlSAException(v1.Size() != v2.Size())(v1.Size())(v2.Size()).msg("vector sizes mismatch");
  return utl::cblas_dot<T>(v1.Size(), v1.GetData(), 1, v2.GetData(), 1);
}

/** \f$ \sum_i v1_i v2_i \f$. No conjugate for complex values. */
template <class T, unsigned int Dim>
inline T
DotProduct ( const NDArray<T,Dim>& v1, const NDArray<T,Dim>& v2 )
{
  utlSAException(v1.Size() != v2.Size())(v1.Size())(v2.Size()).msg("vector sizes mismatch");
  return InnerProduct(v1, v2);
}

template <unsigned int Dim>
inline std::complex<double>
DotProduct( const NDArray<std::complex<double>,Dim>& v1, const NDArray<std::complex<double>,Dim>& v2 )
{
  utlSAException(v1.Size() != v2.Size())(v1.Size())(v2.Size()).msg("vector sizes mismatch");
  std::complex<double> result;
  cblas_zdotu_sub(v1.Size(),v1.GetData(),1,v2.GetData(),1, &result);
  return result;
}

template <class T>
inline void
InnerProduct ( const NDArray<T,2>& mat, const NDArray<T,1>& vec, NDArray<T, 1>& result )
{
  utlSAException(vec.Size()!=mat.Cols())(mat.Cols())(vec.Size()).msg("wrong size");
  utl::ProductUtlMv(mat, vec, result);
}

template <class T>
inline void
InnerProduct ( const NDArray<T,1>& vec, const NDArray<T,2>& mat, NDArray<T, 1>& result )
{
  utlSAException(vec.Size()!=mat.Rows())(mat.Rows())(vec.Size()).msg("wrong size");
  utl::ProductUtlvM(vec, mat, result);
}

template <class T>
inline void
InnerProduct ( const NDArray<T,4>& tensor, const NDArray<T,2>& matrix, NDArray<T, 2>& result )
{
  unsigned int const* shapeT = tensor.GetShape();
  unsigned int const* shapeM = matrix.GetShape();
  utlSAException(shapeT[2]!=shapeM[0] || shapeT[3]!=shapeM[1])
    (shapeT[2])(shapeM[0])(shapeT[3])(shapeM[1]).msg("sizes do not match");

  result.ReSize(shapeT[0], shapeT[1]);

  NDArray<T,2> mat;
  int nn=0;
  for ( int i = 0; i < shapeT[0]; ++i ) 
    {
    for ( int j = 0; j < shapeT[1]; ++j ) 
      {
      mat = tensor.GetRefSubMatrix(i, j);
      result[nn] = utl::InnerProduct(mat, matrix);
      nn++;
      }
    }
}


template <class T>
inline void
InnerProduct ( const NDArray<T,2>& matrix, const NDArray<T,4>& tensor, NDArray<T, 2>& result )
{
  unsigned int const* shapeT = tensor.GetShape();
  unsigned int const* shapeM = matrix.GetShape();
  utlSAException(shapeT[0]!=shapeM[0] || shapeT[1]!=shapeM[1])
    (shapeT[0])(shapeM[0])(shapeT[1])(shapeM[1]).msg("sizes do not match");

  result.ReSize(shapeT[2], shapeT[3]);

  int nn=0;
  for ( int k = 0; k < shapeT[2]; ++k ) 
    {
    for ( int l = 0; l < shapeT[3]; ++l ) 
      {
      result[nn]=0;
      int off=0; 
      for ( int i = 0; i < shapeT[0]; ++i ) 
        {
        for ( int j = 0; j < shapeT[1]; ++j ) 
          {
          result[nn] += tensor(i,j,k,l) * matrix[off];
          off++;
          }
        }
      nn++;
      }
    }
}

template <class T>
inline void
OuterProduct ( const NDArrayBase<T,2>& mat1, const NDArrayBase<T,2>& mat2, NDArray<T, 4>& tensor )
{
  unsigned size[4];
  unsigned const* shape1 = mat1.GetShape();
  unsigned const* shape2 = mat2.GetShape();
  size[0]=shape1[0], size[1]=shape1[1];
  size[2]=shape2[0], size[3]=shape2[1];
  tensor.ReSize(size);

  int off1=0;
  int off0=0;
  for ( int i = 0; i < size[0]; ++i ) 
    {
    for ( int j = 0; j < size[1]; ++j ) 
      {
      int off2=0;
      for ( int k = 0; k < size[2]; ++k ) 
        {
        for ( int l = 0; l < size[3]; ++l ) 
          {
          tensor[off0] = mat1[off1] * mat2[off2];
          off2++;
          off0++;
          }
        }
      off1++;
      }
    }
}

template <class T>
void 
OuterProduct ( const NDArray<T,1>& v1, const NDArray<T,1>& v2, NDArray<T,2>& mat )
{
  bool isComplex = std::is_same<T, std::complex<double>>::value ||std::is_same<T, std::complex<float>>::value; 
  mat.ReSize(v1.Size(), v2.Size());
  for ( int i = 0; i < v1.Size(); ++i ) 
    for ( int j = 0; j < v2.Size(); ++j ) 
      {
      if (isComplex)
        mat(i,j) = v1[i]*std::conj(v2[j]);
      else
        mat(i,j) = v1[i]*v2[j];
      }
}

template <class T>
void 
CrossProduct ( const NDArray<T,1>& v1, const NDArray<T,1>& v2, NDArray<T,1>& v3 )
{
  utlException(v1.Size()!=3, "wrong size");
  utlException(v2.Size()!=3, "wrong size");
  v3.ReSize(3);
  v3[0]=v1[1]*v2[2] - v1[2]*v2[1];
  v3[1]=v1[2]*v2[0] - v1[0]*v2[2];
  v3[2]=v1[0]*v2[1] - v1[1]*v2[0];
}

template< typename T >
inline NDArray<T,2> 
operator*(const NDArray<T,2>& mat1, const NDArray<T,2>& mat2) 
{ 
  utlSAException(mat1.Cols()!=mat2.Rows())(mat1.Cols())(mat2.Rows()).msg("wrong size");
  NDArray<T,2> result;
  utl::ProductUtlMM(mat1, mat2, result);
  return result;
}

template<typename T>
inline NDArray<std::complex<T>,2> 
operator*(const NDArray<T,2>& mat1, const NDArray<std::complex<T>,2>& mat2) 
{ 
  utlSAException(mat1.Cols()!=mat2.Rows())(mat1.Cols())(mat2.Rows()).msg("wrong size");
  NDArray<std::complex<T>,2> result, tmp;
  tmp=mat1;
  utl::ProductUtlMM(tmp, mat2, result);
  return result;
}

template<typename T>
inline NDArray<std::complex<T>,2> 
operator*(const NDArray<std::complex<T>,2>& mat1, const NDArray<T,2>& mat2) 
{ 
  utlSAException(mat1.Cols()!=mat2.Rows())(mat1.Cols())(mat2.Rows()).msg("wrong size");
  NDArray<std::complex<T>,2> result, tmp;
  tmp=mat2;
  utl::ProductUtlMM(mat1, tmp, result);
  return result;
}

template< typename T >
inline NDArray<T,1> 
operator*(const NDArray<T,2>& mat, const NDArray<T,1>& vec) 
{ 
  utlSAException(vec.Size()!=mat.Cols())(mat.Cols())(vec.Size()).msg("wrong size");
  NDArray<T,1> result;
  utl::ProductUtlMv(mat, vec, result);
  return result;
}

template< typename T >
inline NDArray<T,1> 
operator*(const NDArray<T,1>& vec, const NDArray<T,2>& mat) 
{ 
  utlSAException(vec.Size()!=mat.Rows())(mat.Rows())(vec.Size()).msg("wrong size");
  NDArray<T,1> result;
  utl::ProductUtlvM(vec, mat, result);
  return result;
}

template<typename T, typename EType >
inline NDArray<T,1> 
operator*(const NDArray<T,1>& vec, const Expr<EType, typename EType::ValueType>& expr )
{
  utl_shared_ptr<NDArray<T,2> > mat = ToMatrix<T>(expr);
  return vec* (*mat);
}

/** inverse of a general matrix (non-singular), by using LU decomposition. If the matrix is singular, it stop with an error. 
 * It is fast than PInverseVnlMatrix, but only works for non-singular matrix. */
template <class T> 
inline NDArray<T,2>
InverseMatrix( const NDArray<T,2>& mat, const double eps=1e-10 )
{
  utlSAException(mat.Rows()!=mat.Cols())(mat.Rows())(mat.Cols()).msg("the matrix should be symmetric");
  int n = mat.Rows();
  NDArray<T,2> result(n,n);
  utlException(n==0, "matrix cannot be empty");
  if (n>=1 && n <=4)
    {
    utl::InverseSmallMatrix(mat, result, n);
    }
  else
    {
    getri_UtlMatrix(mat, result);
    }
  return result;
}

/** inverse of a symmetric matrix (non-singular), by using LU decomposition. If the matrix is singular, it stop with an error. 
 * It is fast than PInverseVnlMatrix, but only works for non-singular matrix. */
template <class T> 
inline NDArray<T,2>
InverseSymmericMatrix( const NDArray<T,2>& mat, const double eps=1e-10 )
{
  utlSAException(mat.Rows()!=mat.Cols())(mat.Rows())(mat.Cols()).msg("the matrix should be symmetric");
  int n = mat.Rows();
  NDArray<T,2> result;
  if (n<=4)
    result = InverseMatrix(mat, eps);
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
  return result;
}

/** pseudo-inverse of a symmetric matrix which can be singular, by using eigen-decomposition. 
 * If the matrix is not singular, it returns the inverse.  */
template <class T> 
inline NDArray<T,2>
PInverseSymmericMatrix( const NDArray<T,2>& mat, const double eps=1e-10 )
{
  utlSAException(mat.Rows()!=mat.Cols())(mat.Rows())(mat.Cols()).msg("the matrix should be symmetric");
  int N = mat.Rows();
  NDArray<T,2> result;
  NDArray<T,2> eigenVectors, tmp, S(N,N,0.0);
  NDArray<T,1> eigenValues, v;
  mat.EigenDecompositionSymmetricMatrix(eigenValues, eigenVectors);
  for (unsigned i=0; i<N; ++i)
    {
    if ( std::abs(eigenValues[i])>eps ) 
      S(i,i) = 1.0/eigenValues[i];
    else
      S(i,i)=0.0;
    }
  utl::ProductUtlMtM(eigenVectors, S, tmp);
  utl::ProductUtlMM(tmp, eigenVectors, result);
  return result;
}

/** pseudo-inverse of a general matrix, by using SVD  */
template <class T> 
inline NDArray<T,2>
PInverseMatrix( const NDArray<T,2>& mat, const double eps=1e-10 )
{
  NDArray<T,2> U,V, tmp, result;
  NDArray<T,1> uVec, vVec;
  NDArray<utl::remove_complex_t<T>,1> S;
  mat.SVD(U, S, V, 'S');

  NDArray<T,2> diag(S.Size(), S.Size(),0.0);
  for ( int i = 0; i < S.Size(); i += 1 )
    {
    if ( std::abs(S[i])>eps ) 
      diag(i,i) = 1.0/S[i];
    else
      diag(i,i) = 0.0;
    }
  utl::ProductUtlMM(V, diag, tmp);
  utl::ProductUtlMMt(tmp, U, result);
  return result;
}

template <class T>
NDArray<T,1>
ConnectUtlVector ( const NDArray<T,1>& m1, const NDArray<T,1>& m2 )
{
  NDArray<T,1> result(m1.Size()+m2.Size());
  int m1Size = m1.Size();
  for ( int i = 0; i < m1Size; i += 1 ) 
    result[i] = m1[i];
  for ( int i = 0; i < m2.Size(); i += 1 ) 
    result[i+m1Size] = m2[i];
  return result;
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

template <class T, unsigned Dim>
NDArray<T, Dim>
Real ( const NDArray<std::complex<T>, Dim >& mat )
{
  NDArray<T, Dim> result(mat.GetShape());
  for ( int i = 0; i < mat.Size(); ++i ) 
    result[i] = std::real(mat[i]);
  return result;
}

template <class T, unsigned Dim>
NDArray<T, Dim>
Imag ( const NDArray<std::complex<T>, Dim >& mat )
{
  NDArray<T, Dim> result(mat.GetShape());
  for ( int i = 0; i < mat.Size(); ++i ) 
    result[i] = std::imag(mat[i]);
  return result;
}


/** generate complex array from real and imaginary parts  */
template <class T, unsigned Dim>
NDArray<std::complex<T>, Dim>
ComplexCombine ( const NDArray<T,Dim>& arrReal, const NDArray<T,Dim>& arrImg )
{
  utlException(!IsSameShape(arrReal, arrImg), "should have the same shape");
  typedef std::complex<T> VT;
  NDArray<VT, Dim> result(arrReal.GetShape());
  for ( int i = 0; i < result.Size(); ++i ) 
    result[i] = VT(arrReal[i], arrImg[i]);
  return result;
}

/** generate complex array from imaginary part, if val==0  */
template <class T, unsigned Dim>
NDArray<std::complex<T>, Dim>
ComplexCombine ( const T val, const NDArray<T,Dim>& arrImg )
{
  typedef std::complex<T> VT;
  NDArray<VT, Dim> result(arrImg.GetShape());
  for ( int i = 0; i < result.Size(); ++i ) 
    result[i] = VT(val, arrImg[i]);
  return result;
}

/** generate complex array from real part, if val==0  */
template <class T, unsigned Dim>
NDArray<std::complex<T>, Dim>
ComplexCombine ( const NDArray<T,Dim>& arrReal, const T val )
{
  typedef std::complex<T> VT;
  NDArray<VT, Dim> result(arrReal.GetShape());
  for ( int i = 0; i < result.Size(); ++i ) 
    result[i] = VT(arrReal[i], val);
  return result;
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
  tmp = PInverseSymmericMatrix(tmp2);
  ProductUtlMM( QInverse, Aeq, tmp, tmp2);
  ProductUtlMM( tmp2, AeqT, tmp3);
  projMatrix -= tmp3;

  ProductUtlMv( tmp2, beq, projVector);
}

/** 
 * Mean of rotation matrices.  
 *
 * Reference: Rotation Averaging,  International journal of computer vision, 2013 
 * */
template <class T>
utl::NDArray<T,2>
GetMeanOfRotationMatrix(const std::vector<NDArray<T,2> >& matrixVec, const utl::NDArray<T,1>& weights)
{
  int N = matrixVec.size();
  utlSAException(N==0)(N).msg("wrong size");
  utlSAException(weights.Size()!=N)(N)(weights.Size()).msg("wrong size");

  NDArray<T,2> mat(3,3), tangentVec(3,3), matTmp(3,3);
  mat.Fill(0.0);
  mat.FillDiagonal(1.0);

  double eps=1e-4;
  double normVec=1.0;

  do 
    {
    tangentVec.Fill(0.0);
    for ( int i = 0; i < N; ++i ) 
      {
      utl::ProductUtlMtM(mat, matrixVec[i], matTmp);
      matTmp = matTmp.LogM() ;
      tangentVec += matTmp % weights[i];
      }

    mat = mat * tangentVec.ExpM();
    normVec = tangentVec.GetTwoNorm();
    } while ( normVec>eps );

  return mat;
}

template <class T>
utl::NDArray<T,2>
GetMeanOfRotationMatrix(const std::vector<NDArray<T,2> >& matrixVec)
{
  int N = matrixVec.size();
  utlSAException(N==0)(N).msg("wrong size");
  utl::NDArray<T,1>  weights(N, 1.0/N);
  GetMeanOfRotationMatrix(matrixVec, weights);
}

/** 
 * Convert a rotation matrix to axis and angle.
 * 
 * reference: rotm2axang in Matlab. 
 * https://en.wikipedia.org/wiki/Rotation_matrix
 * */
template <class T>
void
RotationMatrixToAxisAngle(const NDArray<T,2>& rotMat, NDArray<T,1>& axis, double& theta)
{
  utlSAException(rotMat.Rows()!=3 || rotMat.Cols()!=3)
    (rotMat.Rows())(rotMat.Cols()).msg("wrong size of rotMat");
  axis.ReSize(3);
  
  theta = std::acos((rotMat(0,0)+rotMat(1,1)+rotMat(2,2)-1.0)/2.0);
  // double thetaSin = std::sin(theta);

  axis[0]= rotMat(2,1)-rotMat(1,2);
  axis[1]= rotMat(0,2)-rotMat(2,0);
  axis[2]= rotMat(1,0)-rotMat(0,1);

  double normAxis = axis.GetTwoNorm();
  if (normAxis>1e-3)
    axis /= normAxis;
  else
    {
    NDArray<T,2> U, V, matTmp=utl::Eye<T>(3)-rotMat; 
    NDArray<T,1> S;
    matTmp.SVD(U, S, V);
    axis = V.GetColumn(2);
    }
}

/** 
 * Convert a set of axis and angle to a rotation matrix.
 * 
 * reference: axang2rotm in Matlab. 
 * https://en.wikipedia.org/wiki/Rotation_matrix
 * */
template <class T>
void
AxisAngleToRotationMatrix(const NDArray<T,1>& axis, const double theta, NDArray<T,2>& rotMat)
{
  NDArray<T,1> v(axis);
  double norm = v.GetTwoNorm();
  utlException(norm<1e-8, "norm is too small");
  v /= norm;

  double c=std::cos(theta), s=std::sin(theta);
  double t=1.0-c; 

  rotMat.ReSize(3,3);
  rotMat(0,0) = t*v[0]*v[0] + c;      rotMat(0,1) = t*v[0]*v[1] - v[2]*s;   rotMat(0,2) = t*v[0]*v[2] + v[1]*s;
  rotMat(1,0) = t*v[0]*v[1] + v[2]*s; rotMat(1,1) = t*v[1]*v[1] + c;        rotMat(1,2) = t*v[1]*v[2] - v[0]*s;
  rotMat(2,0) = t*v[0]*v[2] - v[1]*s; rotMat(2,1) = t*v[1]*v[2] + v[0]*s;   rotMat(2,2) = t*v[2]*v[2] + c;
}

template <class T>
utl::NDArray<T,1> 
MeanDirector ( const std::vector<utl::NDArray<T,1> >& dirVec, const bool isUnitNorm=true)
{
  int N = dirVec.size();
  utlException(N==0, "empty dirVec");

  int dim = dirVec[0].Size();
  utl::NDArray<T,2> mat(dim, dim), eigenVectors, matTmp;
  utl::NDArray<T,1> dir, eigenValues, meanV(dim);
  mat.Fill(0.0);
  for ( int i = 0; i < N; ++i ) 
    {
    dir = dirVec[i];
    utlSAException(dir.Size()!=dim)(dim)(dir.Size()).msg("wrong size of dir");
    utl::OuterProduct(dir, dir, matTmp);
    mat += matTmp;
    }

  if (mat.GetTwoNorm()>1e-10)
    {
    mat.EigenDecompositionSymmetricMatrix(eigenValues, eigenVectors);
    meanV = eigenVectors.GetRow(dim-1);
    if (isUnitNorm)
      return meanV;
    else
      return meanV % eigenValues[dim-1];
    }
  else
    {
    meanV.Fill(0.0);
    return meanV;
    }
}

}

/** @}  */

#endif 
