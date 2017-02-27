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

#include "utlNDArray.h"

namespace utl
{

template < class T, unsigned int Dim >
class NDArray;

template < class T, unsigned int Dim >
class NDArrayBase;


/** \ingroup utlNDArray
* @{ */

/**
 *   \class   NDArray<T,2>
 *   \brief   NDArray<T,2> is a row-major matrix
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \date    08-23-2014
 *   \ingroup utlNDArray Math
 */
template < class T >
class NDArray<T,2> : public NDArrayBase<T,2>
{
public:
  typedef NDArray                  Self;
  typedef NDArrayBase<T,2>            Superclass;

  typedef typename Superclass::ValueType         ValueType;
  typedef typename Superclass::ScalarValueType    ScalarValueType;

  typedef typename Superclass::SizeType          SizeType;
  typedef typename Superclass::ShapeType         ShapeType;
  typedef typename Superclass::Pointer           Pointer;
  typedef typename Superclass::ConstPointer      ConstPointer;
  typedef typename Superclass::Reference         Reference;
  typedef typename Superclass::ConstReference    ConstReference;
  typedef typename Superclass::Iterator          Iterator;
  typedef typename Superclass::ConstIterator     ConstIterator;

  using Superclass::Dimension;
  using Superclass::SubDimension;
  
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
    { 
    }
  NDArray(NDArray<T,2>&& mat) 
    { 
    operator=(std::move(mat)); 
    }

  template<typename EType>
  NDArray(const Expr<EType, typename EType::ValueType>& expr) : Superclass(expr)
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
  
  NDArray<T,2>& operator=(NDArray<T,2> & r)
    {
    Superclass::operator=(r);
    return *this;
    }

  NDArray<T,2>& operator=(NDArray<T,2> && r)
    {
    if ( this != &r ) 
      {
      this->Clear();
      this->Swap(r);
      }
    return *this;
    }

  UTL_ALWAYS_INLINE SizeType Rows() const {return this->m_Shape[0];}
  UTL_ALWAYS_INLINE SizeType Columns() const {return this->m_Shape[1];}
  UTL_ALWAYS_INLINE SizeType Cols() const {return this->m_Shape[1];}
  
  inline bool ReSize(const SizeType rows, const SizeType cols)
    {
    SizeType shape[2];
    shape[0]=rows, shape[1]=cols;
    return Superclass::ReSize(shape);
    }

  UTL_ALWAYS_INLINE T & operator()(unsigned int row, unsigned int col)
  { return this->m_Data[row*this->m_Shape[1]+col];}
  UTL_ALWAYS_INLINE const T & operator()(unsigned int row, unsigned int col) const
  { return this->m_Data[row*this->m_Shape[1]+col];}
  


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
    NDArray<T,2> result(this->m_Shape[1], this->m_Shape[0]);
#ifdef UTL_USE_MKL
    utl::mkl_omatcopy<T>('R', 'T', this->m_Shape[0], this->m_Shape[1], scale, this->m_Data, this->m_Shape[1], result.m_Data, result.m_Shape[1]); 
#else
    for ( int i = 0; i < this->m_Shape[0]; ++i ) 
      for ( int j = 0; j < this->m_Shape[1]; ++j ) 
        result[j*this->m_Shape[0]+i] = (*this)[i*this->m_Shape[1]+j];
    result.Scale(scale);
#endif
    return result;
    }

  NDArray<T,2> GetConjugateTranspose(const T scale=1.0) const 
    {
    NDArray<T, 2> result = GetTranspose(scale);
    result = utl::Conj(result);
    return result;
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
  NDArray<T,1> GetColumn(const int index) const
    {
    NDArray<T,1> vec(this->m_Shape[0]);
    GetColumn(index, vec.GetData());
    return vec;
    }

  NDArray<T,2> GetRows(const std::vector<int>& indexVec) const 
    {
    NDArray<T,2> mat(indexVec.size(), this->m_Shape[1]);
    for ( int i = 0; i < indexVec.size(); ++i ) 
      GetRow(indexVec[i], mat.m_Data+i*this->m_Shape[1]);
    return mat;
    }
  NDArray<T,2> GetNRows(const int index0, const int N) const 
    {
    std::vector<int> indexVec;
    for ( int i = 0; i < N; ++i ) 
      indexVec.push_back(index0+i);
    return GetRows(indexVec);
    }
  NDArray<T,2> GetColumns(const std::vector<int>& indexVec) const 
    {
    NDArray<T,2> mat(this->m_Shape[0], indexVec.size());
    NDArray<T,1> vec;
    for ( int i = 0; i < indexVec.size(); ++i ) 
      {
      vec = GetColumn(indexVec[i]);
      mat.SetColumn(i, vec);
      }
    return mat;
    }
  NDArray<T,2> GetNColumns(const int index0, const int N) const 
    {
    std::vector<int> indexVec;
    for ( int i = 0; i < N; ++i ) 
      indexVec.push_back(index0+i);
    return GetColumns(indexVec);
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
  NDArray<T,2> GetCrop(const int x0, const int x1, const int y0, const int y1) const
    {
    utlSAException(x0<0 || x0>=this->m_Shape[0] || x1<0 || x1>=this->m_Shape[0] || x0>x1)(x0)(x1)(this->m_Shape[0]).msg("wrong index");
    utlSAException(y0<0 || y0>=this->m_Shape[0] || y1<0 || y1>=this->m_Shape[1] || y0>y1)(y0)(y1)(this->m_Shape[1]).msg("wrong index");
    NDArray<T,2> mat(x1-x0+1, y1-y0+1);
    for ( int i = 0; i < mat.Rows(); ++i ) 
      for ( int j = 0; j < mat.Cols(); ++j ) 
        mat(i,j) = (*this)(i+x0, j+x1);
    return mat;
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
  void SVD(NDArray<T,2>& U, NDArray<ScalarValueType,1>& S, NDArray<T,2>& V, char format='S') const
    {
#ifdef UTL_USE_FASTLAPACK
  // NOTE: fast but may be wrong for big matrix and multi-thread if openblas is not correctly built
  utl::gesdd_UtlMatrix(*this, U, S, V, format);
#else
  // slow but robust
  utl::gesvd_UtlMatrix<T>(*this, U, S, V, format);
#endif
    }


  /** 
   * \brief Eigen-decomposition for symmetric matrix.  
   *
   * \param eigenValues Eigen-values are in \b increasing order. 
   * \param eigenVectors Eigen-vectors. each row is an eigen-vector
   **/
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
  
  void EigenDecompositionNonSymmetricMatrix (NDArray<T,1>& valReal, NDArray<T,1>& valImg) const
    {
    utl::geev_UtlMatrix(*this, valReal, valImg);
    }
  void EigenDecompositionNonSymmetricMatrix (NDArray<T,1>& valReal, NDArray<T,1>& valImg, NDArray<T,2>& vecRealR, NDArray<T,2>& vecImgR) const
    {
    utl::geev_UtlMatrix(*this, valReal, valImg, vecRealR, vecImgR);
    }
  void EigenDecompositionNonSymmetricMatrix (NDArray<T,1>& valReal, NDArray<T,1>& valImg, NDArray<T,2>& vecRealR, NDArray<T,2>& vecImgR, NDArray<T,2>& vecRealL, NDArray<T,2>& vecImgL) const
    {
    utl::geev_UtlMatrix(*this, valReal, valImg, vecRealR, vecImgR, vecRealL, vecImgL);
    }
  
  NDArray<T,2> InverseSymmericMatrix(const double eps=1e-10)
    {
    return utl::InverseSymmericMatrix(*this, eps);
    }

  NDArray<T,2> PInverseSymmericMatrix(const double eps=1e-10)
    {
    return utl::PInverseSymmericMatrix(*this, eps);
    }
  
  NDArray<T,2> PInverseMatrix(const double eps=1e-10)
    {
    return utl::PInverseMatrix(*this, eps);
    }

  NDArray<T,2> InverseMatrix(const double eps=1e-10)
    {
    return utl::InverseMatrix(*this, eps);
    }

  T Determinant() const
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
          const T
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
        T v1 = (*this)(i,j), v2 = (*this)(j,i);
        if (std::abs(v1-v2)>eps*std::abs(v1))
          return false;
        }
    return true;
    }

  void Symmetrize()
    {
    utlException(this->m_Shape[0]!=this->m_Shape[1], "rows!=cols");
    for ( int i = 0; i < this->m_Shape[0]; ++i ) 
      for ( int j = 0; j < i; ++j ) 
        {
        T v1 = (*this)(i,j), v2 = (*this)(j,i);
        (*this)(i,j) = (*this)(j,i) = (v1+v2)*0.5;
        }
    }

  /** Assume the matrix is symmetric  */
  NDArray<T,2> ExpM()
    {
    NDArray<T,1> eigenValues;
    NDArray<T,2> eigenVectors, result;
    
    EigenDecompositionSymmetricMatrix(eigenValues, eigenVectors);
    for ( int i = 0; i < eigenValues.Size(); ++i ) 
      eigenValues[i] = std::exp(eigenValues[i]);

    result = eigenVectors.GetTranspose() * eigenValues.GetDiagonalMatrix() * eigenVectors;
    return result;
    }
  
  /** Assume the matrix is symmetric  */
  NDArray<T,2> LogM()
    {
    NDArray<T,1> eigenValues;
    NDArray<T,2> eigenVectors, result;
    
    EigenDecompositionSymmetricMatrix(eigenValues, eigenVectors);
    for ( int i = 0; i < eigenValues.Size(); ++i ) 
      {
      utlSAException(eigenValues[i]<=0)(i)(eigenValues[i]).msg("negative eigenValue");
      eigenValues[i] = std::log(eigenValues[i]);
      }

    result = eigenVectors.GetTranspose() * eigenValues.GetDiagonalMatrix() * eigenVectors;
    return result;
    }

  double GetArrayOneNorm() const { return Superclass::GetOneNorm(); }
  double GetArrayTwoNorm() const { return Superclass::GetTwoNorm(); }
  double GetArrayInfNorm() const { return Superclass::GetInfNorm(); }
  int GetArrayZeroNorm() const { return Superclass::GetZeroNorm(); }

  double GetTwoNorm() const {return GetArrayTwoNorm();}
  double GetOneNorm() const
    { return utl::lange(LAPACK_ROW_MAJOR, 'O', this->m_Shape[0], this->m_Shape[1], this->m_Data, this->m_Shape[1]); }
  double GetInfNorm() const
    { return utl::lange(LAPACK_COL_MAJOR, 'O', this->m_Shape[1], this->m_Shape[0], this->m_Data, this->m_Shape[1]); }

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


/** @}  */

}

#endif 
