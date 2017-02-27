/**
 *       @file  utlMatrixGTest.cxx
 *      @brief  gtest for utl::Matrix
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-23-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlGTest.h"
#include "utlNDArray.h"
#include "utlVNL.h"
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_svd.h>
#include "utlMath.h"
#include "utlDMRI.h"


typedef utl::NDArray<double,1> UtlVectorType;
typedef utl::NDArray<double,2> UtlMatrixType;

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const double val1, const double val2 )
{
  vnl_vector<T> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<T>(val1,val2);
  return v1;
}

template <>
vnl_vector<std::complex<double> > 
__GenerateRandomVector<std::complex<double> > ( const int N, const double val1, const double val2 )
{
  vnl_vector<std::complex<double> > v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    {
    v1[i].real( utl::Random<double>(val1,val2) );
    v1[i].imag( utl::Random<double>(val1,val2) );
    }
  return v1;
}


template <class T>
UtlVectorType
__GenerateUtlVector(const int N, const T val1, const T val2)
{
  UtlVectorType v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

template <class T>
UtlMatrixType
__GenerateUtlMatrix(const int rows, const int cols, const T val1, const T val2)
{
  UtlMatrixType mat(rows, cols);
  for ( int i = 0; i < rows*cols; i += 1 ) 
    mat[i] = utl::Random<double>(val1,val2);
  return mat;
}

TEST(utlMatrix, Constructors)
{
  int Rows=3, Cols=4;
  int N=Rows*Cols;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_matrix<double> mat0(vec0.data_block(),Rows,Cols);
    {
    UtlMatrixType mat(vec0.data_block(),Rows,Cols);
    EXPECT_NEAR_VECTOR(vec0, mat, N, 1e-10);
    EXPECT_NEAR_MATRIX(mat0, mat, Rows, Cols, 1e-10);
    }
    {
    UtlMatrixType mat, mat2;
    mat.SetData(vec0.data_block(), Rows, Cols);
    EXPECT_NEAR_VECTOR(vec0, mat, N, 1e-10);
    EXPECT_NEAR_MATRIX(mat0, mat, Rows, Cols, 1e-10);
    EXPECT_EQ((void*)vec0.data_block(), (void*)mat.GetData());
    mat2.CopyData(vec0.data_block(), Rows, Cols);
    EXPECT_NEAR_VECTOR(vec0, mat2, N, 1e-10);
    EXPECT_NEAR_MATRIX(mat0, mat2, Rows, Cols, 1e-10);
    EXPECT_NE((void*)vec0.data_block(), (void*)mat2.GetData());
    }
    {
    UtlMatrixType mat(vec0.data_block(),Rows,Cols), mat2, mat3(mat);
    mat2 = mat;
    EXPECT_NEAR_MATRIX(mat0, mat2, Rows, Cols, 1e-10);
    EXPECT_NEAR_MATRIX(mat0, mat3, Rows, Cols, 1e-10);
    }
    {
    UtlMatrixType mat(vec0.data_block(),Rows,Cols);
    utl::NDArray<int,2> mat2(mat), mat3;
    mat3 = mat;
    for ( int i=0; i<mat.Size(); ++i ) 
      {
      EXPECT_NEAR((int)mat[i], mat2[i], 1e-10);
      EXPECT_NEAR((int)mat[i], mat3[i], 1e-10);
      }
    }
    {
    double val = utl::Random<double>(-2.0,2.0);
    UtlMatrixType mat(Rows,Cols, val);
    for ( UtlMatrixType::Iterator iter = mat.Begin();  iter!=mat.End(); ++iter ) 
      EXPECT_NEAR(val, *iter, 1e-10);
    mat.Fill(3.0);
    for ( UtlMatrixType::Iterator iter = mat.Begin();  iter!=mat.End(); ++iter ) 
      EXPECT_NEAR(3.0, *iter, 1e-10);
    }
}

TEST(utlMatrix, Operators)
{
  int Rows=3, Cols=4, K=5;
  int N=Rows*Cols, Nk=Cols*K;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> vecCol = __GenerateRandomVector<double>(Cols, -2.0, 2.0);
  vnl_vector<double> vecRow = __GenerateRandomVector<double>(Rows, -2.0, 2.0);
  vnl_matrix<double> mat0(vec0.data_block(),Rows,Cols);
  vnl_vector<double> vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_matrix<double> mat1(vec1.data_block(),Rows,Cols);
  vnl_vector<double> veck = __GenerateRandomVector<double>(Nk, -2.0, 2.0);
  vnl_matrix<double> matk(veck.data_block(),Cols,K);

  UtlMatrixType mat0Utl(vec0.data_block(),Rows,Cols);
  UtlMatrixType mat1Utl(vec1.data_block(),Rows,Cols);
  UtlMatrixType matkUtl(veck.data_block(),Cols,K);
  UtlVectorType vecColUtl(vecCol.data_block(),Cols);
  UtlVectorType vecRowUtl(vecRow.data_block(),Rows);

  vnl_matrix<double> mat2;
  UtlMatrixType mat2Utl;

    {
    // matrix + matrix
    mat2 = mat0+mat1;
    mat2Utl = mat0Utl+mat1Utl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix - matrix
    mat2 = mat0-mat1;
    mat2Utl = mat0Utl-mat1Utl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix / matrix
    mat2 = element_quotient(mat0, mat1);
    mat2Utl = mat0Utl/mat1Utl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix * matrix (element wize)
    mat2 = element_product(mat0, mat1);
    mat2Utl = mat0Utl%mat1Utl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix * matrix
    mat2 = mat0*matk;
    mat2Utl = mat0Utl*matkUtl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // vector * matrix
    vnl_vector<double> vecProd = vecRow*mat0;
    UtlVectorType vecUtlProd = vecRowUtl*mat0Utl;
    EXPECT_NEAR_VECTOR(vecProd, vecUtlProd, Cols, 1e-10);

    vecProd = vecRow*(mat0+mat0);
    vecUtlProd = vecRowUtl*(mat0Utl+mat0Utl);
    EXPECT_NEAR_VECTOR(vecProd, vecUtlProd, Cols, 1e-10);
    }

    {
    // matrix * vector
    vnl_vector<double> vecProd = mat0*vecCol;
    UtlVectorType vecUtlProd = mat0Utl*vecColUtl;
    // UtlMatrixType matUtlProd = mat0Utl*vecColUtl;
    EXPECT_NEAR_VECTOR(vecProd, vecUtlProd, Rows, 1e-10);
    // EXPECT_NEAR_VECTOR(vecProd, matUtlProd, Rows, 1e-10);
    // EXPECT_EQ(Rows, matUtlProd.Rows());
    // EXPECT_EQ(1, matUtlProd.Cols());
    }

    {
    // scales
    mat2 = mat0*3.0;
    mat2Utl= mat0Utl%3.0;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);

    mat2 = 3.0*mat0;
    mat2Utl= 3.0%mat0Utl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    
    mat2 = 3.0+mat0;
    mat2Utl= 3.0+mat0Utl;
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    
    mat2Utl= 3.0/(3.0+mat0Utl);
    for ( int i = 0; i < Rows; ++i ) 
      for ( int j = 0; j < Cols; ++j ) 
        EXPECT_NEAR(3.0/(3.0+mat0Utl(i,j)), mat2Utl(i,j), 1e-10);
    }

    {
    // in-place transpose
    mat2 = mat0;
    mat2Utl = mat0Utl;
    mat2.inplace_transpose();
    mat2Utl.TransposeInplace();
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Cols, Rows, 1e-10);

    // out-place transpose
    mat2 = 3.0*mat0.transpose();
    mat2Utl = mat0Utl.GetTranspose(3.0);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Cols, Rows, 1e-10);
    }

    {
    // outer product
    utl::OuterProduct(vecRowUtl, vecColUtl, mat2Utl);
    UtlVectorType vv, v1;
    utl::ProductUtlMv(mat2Utl, vecColUtl, vv);
    v1 = vecRowUtl%(vecColUtl.GetSquaredTwoNorm());
    EXPECT_NEAR_VECTOR(vv, v1, Rows, 1e-10);
    }
}

TEST(utlMatrix, OperatorsComplex)
{
  typedef std::complex<double> VT;
  int Rows=3, Cols=4, K=5;
  int N=Rows*Cols, Nk=Cols*K;
  vnl_vector<VT> vec0 = __GenerateRandomVector<VT>(N, -2.0, 2.0);
  vnl_vector<VT> vecCol = __GenerateRandomVector<VT>(Cols, -2.0, 2.0);
  vnl_vector<VT> vecRow = __GenerateRandomVector<VT>(Rows, -2.0, 2.0);
  vnl_matrix<VT> mat0(vec0.data_block(),Rows,Cols);
  vnl_vector<VT> vec1 = __GenerateRandomVector<VT>(N, -2.0, 2.0);
  vnl_matrix<VT> mat1(vec1.data_block(),Rows,Cols);
  vnl_vector<VT> veck = __GenerateRandomVector<VT>(Nk, -2.0, 2.0);
  vnl_matrix<VT> matk(veck.data_block(),Cols,K);

  utl::Matrix<VT> mat0Utl(vec0.data_block(),Rows,Cols);
  utl::Matrix<VT> mat1Utl(vec1.data_block(),Rows,Cols);
  utl::Matrix<VT> matkUtl(veck.data_block(),Cols,K);
  utl::Vector<VT> vecColUtl(vecCol.data_block(),Cols);
  utl::Vector<VT> vecRowUtl(vecRow.data_block(),Rows);

  vnl_matrix<VT> mat2;
  utl::Matrix<VT> mat2Utl;

  utl::Matrix<double> mat0UtlReal = utl::Real(mat0Utl);
  utl::Matrix<double> matkUtlReal = utl::Real(matkUtl);
  vnl_matrix<VT> mat0Real(mat0), matkReal(matk);
  std::transform(mat0.begin(),mat0.end(), mat0Real.begin(), 
        []( VT& val) {return VT(std::real(val), 0.0); } );
  std::transform(matk.begin(),matk.end(), matkReal.begin(), 
        []( VT& val) {return VT(std::real(val), 0.0); } );


    {
    // matrix + matrix
    mat2 = mat0+mat1;
    mat2Utl = mat0Utl+mat1Utl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix - matrix
    mat2 = mat0-mat1;
    mat2Utl = mat0Utl-mat1Utl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix / matrix
    mat2 = element_quotient(mat0, mat1);
    mat2Utl = mat0Utl/mat1Utl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix * matrix (element wize)
    mat2 = element_product(mat0, mat1);
    mat2Utl = mat0Utl%mat1Utl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // matrix * matrix
    mat2 = mat0*matk;
    mat2Utl = mat0Utl*matkUtl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);

    mat2 = mat0Real*matk;
    mat2Utl = mat0UtlReal*matkUtl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    
    mat2 = mat0*matkReal;
    mat2Utl = mat0Utl*matkUtlReal;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    }

    {
    // vector * matrix
    vnl_vector<VT> vecProd = vecRow*mat0;
    utl::Vector<VT> vecUtlProd = vecRowUtl*mat0Utl;
    EXPECT_NEAR_VECTOR_COMPLEX(vecProd, vecUtlProd, Cols, 1e-10);
    
    vecProd = vecRow*(mat0+mat0);
    vecUtlProd = vecRowUtl*(mat0Utl+mat0Utl);
    EXPECT_NEAR_VECTOR_COMPLEX(vecProd, vecUtlProd, Cols, 1e-10);
    }

    {
    // matrix * vector
    vnl_vector<VT> vecProd = mat0*vecCol;
    utl::Vector<VT> vecUtlProd = mat0Utl*vecColUtl;
    // UtlMatrixType matUtlProd = mat0Utl*vecColUtl;
    EXPECT_NEAR_VECTOR_COMPLEX(vecProd, vecUtlProd, Rows, 1e-10);
    // EXPECT_NEAR_VECTOR(vecProd, matUtlProd, Rows, 1e-10);
    // EXPECT_EQ(Rows, matUtlProd.Rows());
    // EXPECT_EQ(1, matUtlProd.Cols());
    }

    {
    // scales
    mat2 = mat0*3.0;
    mat2Utl= mat0Utl%3.0;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);

    mat2 = VT(3.0,0.0)*mat0;
    mat2Utl= 3.0%mat0Utl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    
    mat2 = VT(3.0,0.0)+mat0;
    mat2Utl= 3.0+mat0Utl;
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Rows, Cols, 1e-10);
    
    mat2Utl= 3.0/(3.0+mat0Utl);
    for ( int i = 0; i < Rows; ++i ) 
      for ( int j = 0; j < Cols; ++j ) 
        EXPECT_NEAR_COMPLEX(3.0/(3.0+mat0Utl(i,j)), mat2Utl(i,j), 1e-10);
    }

    {
    // in-place transpose
    mat2 = mat0;
    mat2Utl = mat0Utl;
    mat2.inplace_transpose();
    mat2Utl.TransposeInplace();
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Cols, Rows, 1e-10);

    // out-place transpose
    mat2 = VT(3.0,0.0)*mat0.transpose();
    mat2Utl = mat0Utl.GetTranspose(3.0);
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Cols, Rows, 1e-10);
    
    mat2 = VT(3.0,0.0)*mat0.conjugate_transpose();
    mat2Utl = mat0Utl.GetConjugateTranspose(VT(3.0,0.0));
    EXPECT_NEAR_MATRIX_COMPLEX(mat2, mat2Utl, Cols, Rows, 1e-10);
    }

    {
    // outer product
    utl::OuterProduct(vecRowUtl, vecColUtl, mat2Utl);
    utl::Vector<VT> vv, v1;
    utl::ProductUtlMv(mat2Utl, vecColUtl, vv);
    v1 = vecRowUtl%(vecColUtl.GetSquaredTwoNorm());
    EXPECT_NEAR_VECTOR_COMPLEX(vv, v1, Rows, 1e-10);
    }
}

TEST(utlMatrix, SetGetRowColumn)
{
  int Rows=3, Cols=4;
  int N=Rows*Cols;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> vecCol = __GenerateRandomVector<double>(Cols, -2.0, 2.0);
  vnl_vector<double> vecRow = __GenerateRandomVector<double>(Rows, -2.0, 2.0);
  vnl_matrix<double> mat0(vec0.data_block(),Rows,Cols);

  UtlMatrixType mat0Utl(vec0.data_block(),Rows,Cols);
  UtlVectorType vecColUtl(vecCol.data_block(),Cols);
  UtlVectorType vecRowUtl(vecRow.data_block(),Rows);

  vnl_matrix<double> mat2;
  UtlMatrixType mat2Utl;
  vnl_vector<double> vec2;
  UtlVectorType vec2Utl;
    {
    mat2 = mat0; 
    mat2.set_row(1,vecCol);
    mat2Utl = mat0Utl; 
    mat2Utl.SetRow(1,vecColUtl);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    vec2 = mat2.get_row(2);
    vec2Utl = mat2Utl.GetRow(2);
    EXPECT_NEAR_VECTOR(vec2, vec2Utl, Cols, 1e-10);

    mat2 = mat0; 
    mat2.set_column(1,vecRow);
    mat2Utl = mat0Utl; 
    mat2Utl.SetColumn(1,vecRowUtl);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    vec2 = mat2.get_column(2);
    vec2Utl = mat2Utl.GetColumn(2);
    EXPECT_NEAR_VECTOR(vec2, vec2Utl, Rows, 1e-10);
    }
    {
    std::vector<int> indexVec(2);
    indexVec[0]=0; indexVec[1]=2;
    mat2 = utl::GetRowsVnlMatrix(mat0, indexVec);
    mat2Utl = mat0Utl.GetRows(indexVec);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, indexVec.size(), Cols, 1e-10);
    mat2 = utl::GetColumnsVnlMatrix(mat0, indexVec);
    mat2Utl = mat0Utl.GetColumns(indexVec);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, indexVec.size(), 1e-10);
    }
}

TEST(utlMatrix, SVDInverse)
{
  int Rows=3, Cols=4;
  int N=Rows*Cols;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0), s1;
  vnl_matrix<double> mat0(vec0.data_block(),Rows,Cols);
  UtlMatrixType mat0Utl(vec0.data_block(),Rows,Cols);

  vnl_matrix<double> U, V, matEst; 
  vnl_svd<double> svd(mat0); 
  matEst = svd.U()*svd.W()*svd.V().transpose(); 
  s1=svd.W().get_diagonal();
  EXPECT_NEAR_VNLMATRIX(mat0, matEst, 1e-10); 

  UtlMatrixType UUtl, VUtl, matEstUtl; 
  UtlVectorType SUtl;
  utl::gesvd_UtlMatrix(mat0Utl, UUtl, SUtl, VUtl); 
  matEstUtl = UUtl*SUtl.GetDiagonalMatrix()*VUtl.GetTranspose(); 
  EXPECT_NEAR_MATRIX(mat0, matEstUtl, mat0.rows(), mat0.cols(), 1e-10); 
  EXPECT_NEAR_VECTOR(s1, SUtl, utl::min(SUtl.Size(),s1.size()),1e-10); 
}

TEST(utlMatrix, Norms)
{
  int Rows=3, Cols=4;
  int N=Rows*Cols;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0), s1;
  vnl_matrix<double> mat0(vec0.data_block(),Rows,Cols);
  UtlMatrixType mat0Utl(vec0.data_block(),Rows,Cols);

  double norm, normUtl;
    {
    norm = mat0.operator_one_norm();
    normUtl = mat0Utl.GetOneNorm();
    EXPECT_NEAR(norm, normUtl, std::fabs(norm)*1e-10);
    }
    {
    norm = mat0.operator_inf_norm();
    normUtl = mat0Utl.GetInfNorm();
    EXPECT_NEAR(norm, normUtl, std::fabs(norm)*1e-10);
    }
    {
    norm = mat0.frobenius_norm();
    normUtl = mat0Utl.GetTwoNorm();
    EXPECT_NEAR(norm, normUtl, std::fabs(norm)*1e-10);
    }
    {
    norm = (mat0+mat0*3.0).frobenius_norm();
    utl_shared_ptr<UtlMatrixType > pp=utl::ToMatrix<double>(mat0Utl+mat0Utl%3.0);
    normUtl = utl::ToMatrix<double>(mat0Utl+mat0Utl%3.0)->GetTwoNorm();
    EXPECT_NEAR(norm, normUtl, std::fabs(norm)*1e-10);
    }
}

TEST(utlMatrix, Det_Inv_smallMatrix)
{
  int N=4;
  for ( int n = 1; n <= N; ++n ) 
    {
    vnl_vector<double> vec0 = __GenerateRandomVector<double>(n*n, -2.0, 2.0);
    vnl_matrix<double> mat0(vec0.data_block(), n,n);
    UtlMatrixType mat0Utl(vec0.data_block(),n,n);

    double norm = vnl_determinant(mat0);
    double normUtl = mat0Utl.Determinant();
    EXPECT_NEAR(norm, normUtl, std::fabs(norm)*1e-8);
    normUtl = utl::DeterminantSmallMatrix(mat0Utl, n);
    EXPECT_NEAR(norm, normUtl, std::fabs(norm)*1e-8);

    vnl_matrix<double> mat0_inv = utl::GetVnlMatrixPInverse(mat0,1e-10);
    UtlMatrixType mat0Utl_inv(n,n);
    mat0Utl_inv = mat0Utl.InverseMatrix(1e-10);
    EXPECT_NEAR_MATRIX(mat0_inv, mat0Utl_inv, n, n, 1e-10);
    }
}

TEST(utlMatrix, ConvertDTI)
{
  UtlMatrixType mat1, mat2, mat3;
  UtlVectorType v1, v2, v3;

    {
    // convert back
    v1 = __GenerateUtlVector<double>(6,-2.0,2.0);
    mat1.ReSize(3,3);
    v2.ReSize(6);

    utl::ConvertTensor6DTo9D(v1, mat1, TENSOR_LOWER_TRIANGULAR);
    utl::ConvertTensor9DTo6D(mat1, v2, TENSOR_LOWER_TRIANGULAR);
    EXPECT_NEAR_VECTOR(v1, v2, 6, 1e-10);

    utl::ConvertTensor6DTo9D(v1, mat1, TENSOR_UPPER_TRIANGULAR);
    utl::ConvertTensor9DTo6D(mat1, v2, TENSOR_UPPER_TRIANGULAR);
    EXPECT_NEAR_VECTOR(v1, v2, 6, 1e-10);

    utl::ConvertTensor6DTo9D(v1, mat1, TENSOR_EMBED6D);
    utl::ConvertTensor9DTo6D(mat1, v2, TENSOR_EMBED6D);
    EXPECT_NEAR_VECTOR(v1, v2, 6, 1e-10);
    
    utl::Convert1To2Tensor(v1, mat1);
    utl::Convert2To1Tensor(mat1, v2);
    EXPECT_NEAR_VECTOR(v1, v2, 6, 1e-10);
    }

    {
    // keep inner product: mat [*] mat
    mat1 = __GenerateUtlMatrix<double>(3,3,-2.0,2.0);
    mat2 = __GenerateUtlMatrix<double>(3,3,-2.0,2.0);
    mat1.Symmetrize();
    mat2.Symmetrize();
    v1.ReSize(6);
    v2.ReSize(6);
    utl::Convert2To1Tensor(mat1, v1);
    utl::Convert2To1Tensor(mat2, v2);
    EXPECT_NEAR(utl::InnerProduct(mat1, mat2), utl::InnerProduct(v1,v2), 1e-10);
    }
    
    {
    // keep inner product: vec [*] mat [*] vec
    mat1 = __GenerateUtlMatrix<double>(3,3,-2.0,2.0);
    v2 = __GenerateUtlVector<double>(3,-2.0,2.0);
    v3 = __GenerateUtlVector<double>(3,-2.0,2.0);
    mat1.Symmetrize();
    v1.ReSize(6);

    utl::Convert2To1Tensor(mat1, v1);
    utl::Convert1To2Tensor(v1, mat2);
    EXPECT_NEAR_UTLMATRIX(mat1, mat2, 1e-10);

    utl::OuterProduct(v2, v3, mat2);
    EXPECT_NEAR(utl::InnerProduct(mat1, mat2), utl::InnerProduct(v2, mat1*v3), 1e-10);

    UtlVectorType v4;
    utl::OuterProduct(v2, v2, mat2);
    utl::Convert2To1Tensor(mat2, v4);
    utl::Convert1To2Tensor(v4, mat3);
    EXPECT_NEAR_UTLMATRIX(mat3, mat2, 1e-10);
    EXPECT_NEAR(utl::InnerProduct(mat1, mat2), utl::InnerProduct(v1, v4), 1e-10);
    }
}

TEST(utlMatrix, RotationMatrix)
{
  double theta=0, theta1;
  UtlVectorType axis(3), axis1(3);
  UtlMatrixType mat(3,3), mat1(3,3);
  
    {
    mat.Fill(0.0);
    mat(0,0)=1.0, mat(1,1)=-1.0, mat(2,2)=-1.0;
    utl::RotationMatrixToAxisAngle(mat, axis, theta);
    EXPECT_NEAR(1.0, axis[0], 1e-10 );
    EXPECT_NEAR(0.0, axis[1], 1e-10 );
    EXPECT_NEAR(0.0, axis[2], 1e-10 );
    EXPECT_NEAR(M_PI, theta, 1e-10 );
    }

    {
    mat.Fill(0.0);
    mat(0,0)=1.0, mat(1,1)=-1.0, mat(2,2)=-1.0;
    utl::RotationMatrixToAxisAngle(mat, axis, theta);
    EXPECT_NEAR(1.0, axis[0], 1e-10 );
    EXPECT_NEAR(0.0, axis[1], 1e-10 );
    EXPECT_NEAR(0.0, axis[2], 1e-10 );
    EXPECT_NEAR(M_PI, theta, 1e-10 );
    
    utl::AxisAngleToRotationMatrix(axis, theta, mat1);
    EXPECT_NEAR_UTLMATRIX(mat, mat1, 1e-10);
    }

  int N=10;
  for ( int i = 0; i < N; ++i ) 
    {
    for ( int j = 0; j < 3; ++j ) 
      axis[j] = utl::Random<double>(-1.0,1.0);
    axis /= axis.GetTwoNorm();
    if (axis[2]<0)
      axis %= -1.0;
    theta = utl::Random<double>(-M_PI, M_PI);
    utl::AxisAngleToRotationMatrix(axis, theta, mat);

    utl::RotationMatrixToAxisAngle(mat, axis1, theta1);
    if (axis1[2]<0)
      {
      axis1 %= -1.0;
      theta1 *= -1.0;
      }
    EXPECT_NEAR(theta, theta1, 1e-10 );
    for ( int j = 0; j < 3; ++j ) 
      EXPECT_NEAR(axis[j], axis1[j], 1e-10 ) << "axis1 = " << axis1 << ", axis=" << axis;
    
    utl::AxisAngleToRotationMatrix(axis1, theta1, mat1);
    EXPECT_NEAR_UTLMATRIX(mat, mat1, 1e-10);
    }

}

