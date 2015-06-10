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


typedef utl::NDArray<double,1> UtlVectorType;
typedef utl::NDArray<double,2> UtlMatrixType;

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const T val1, const T val2 )
{
  vnl_vector<double> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
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
    mat0Utl.GetTranspose(mat2Utl,3.0);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Cols, Rows, 1e-10);
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
    mat2Utl.GetRow(2,vec2Utl);
    EXPECT_NEAR_VECTOR(vec2, vec2Utl, Cols, 1e-10);

    mat2 = mat0; 
    mat2.set_column(1,vecRow);
    mat2Utl = mat0Utl; 
    mat2Utl.SetColumn(1,vecRowUtl);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, Rows, Cols, 1e-10);
    vec2 = mat2.get_column(2);
    mat2Utl.GetColumn(2,vec2Utl);
    EXPECT_NEAR_VECTOR(vec2, vec2Utl, Rows, 1e-10);
    }
    {
    std::vector<int> indexVec(2);
    indexVec[0]=0; indexVec[1]=2;
    mat2 = utl::GetRowsVnlMatrix(mat0, indexVec);
    mat0Utl.GetRows(indexVec, mat2Utl);
    EXPECT_NEAR_MATRIX(mat2, mat2Utl, indexVec.size(), Cols, 1e-10);
    mat2 = utl::GetColumnsVnlMatrix(mat0, indexVec);
    mat0Utl.GetColumns(indexVec, mat2Utl);
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
    mat0Utl.InverseMatrix(mat0Utl_inv, 1e-10);
    EXPECT_NEAR_MATRIX(mat0_inv, mat0Utl_inv, n, n, 1e-10);
    }
}



