/**
 *       @file  utlVNLLapackGTest.cxx
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

#include "utlGTest.h"
// #include "utl.h"
#include "utlVNLLapack.h"
#include "utlSTDHeaders.h"
#include "itkMultiThreader.h"
#include <vnl/algo/vnl_svd.h>
#include "utlNDArray.h"
#include "utlVNLIO.h"

typedef vnl_matrix<double>         MatrixType;
typedef utl::NDArray<double,2>     UtlMatrixType;
typedef utl_shared_ptr<MatrixType> MatrixPointer;

template <class T>
vnl_matrix<T> 
__GenerateRandomMatrix ( const int M, const int N, const T val1, const T val2 )
{
  vnl_matrix<double> mat(M, N);
  for ( int i = 0; i < mat.rows(); i += 1 ) 
    for ( int j = 0; j < mat.columns(); j += 1 ) 
      mat(i,j) = utl::Random<double>(val1,val2);
  return mat;
}

class utlVNLLapack_gtest : public testing::Test 
{
protected:

  virtual void SetUp() 
    {
    const double array [] =
      {
      16,   3.5,  6,    8.5,
      3.5,  11,   8.5,  11,
      6,    8.5,  6,    13.5,
      8.5,  11,   13.5, 1
      };
    const double array2 [] =
      {
      16,   2,  3,    15,
      5,  11,   10,  8,
      9,    7,  6,    12
      };
    matrixSym = MatrixType(array, 4, 4);
    matrixNonSym = MatrixType(array2, 3, 4);

    int m = utl::RandomInt(3,10);
    int n = utl::RandomInt(3,10);
    int min_mn = utl::min(m,n);
    matrixRandomNonSym = MatrixType(m,n,0.0);
    matrixRandomSym = MatrixType(min_mn,min_mn,0.0);
    for ( int i = 0; i < m; i += 1 ) 
      for ( int j = 0; j < n; j += 1 ) 
        {
        matrixRandomNonSym(i,j) = utl::Random<double>(-2.0,2.0);
        if (i<min_mn && j<min_mn && i<=j)
          {
          matrixRandomSym(i,j) = matrixRandomNonSym(i,j);
          matrixRandomSym(j,i) = matrixRandomNonSym(i,j);
          }
        }
    }
  

  MatrixType matrixSym;
  MatrixType matrixNonSym;
  MatrixType matrixRandomSym;
  MatrixType matrixRandomNonSym;
};


TEST_F(utlVNLLapack_gtest, dsyev_VnlMatrix)
{
  std::vector<MatrixType > matVec;
  matVec.push_back(matrixSym);
  matVec.push_back(matrixRandomSym);
  for ( int i = 0; i < matVec.size(); i += 1 ) 
    {
    typedef vnl_symmetric_eigensystem< double >  SymEigenSystemType; 
    SymEigenSystemType eig (matrixSym); 
    MatrixType eigenVectors, eigenVectors1; 
    vnl_vector<double> eigenValues, eigenValues1; 
    utl::syev_VnlMatrix(matrixSym, eigenValues, eigenVectors); 
    utl::syevd_VnlMatrix(matrixSym, eigenValues1, eigenVectors1); 
    for ( int j = 0; j < matrixSym.rows(); j += 1 )  
      { 
      EXPECT_NEAR(eig.D(j,j), eigenValues(j), 1e-10); 
      EXPECT_NEAR(eig.D(j,j), eigenValues1(j), 1e-10); 
      vnl_vector<double> v1, v2; 
      v1 = eig.V.get_column(j); 
      v2 = eigenVectors.get_row(j); 
      EXPECT_NEAR(1.0, std::abs(dot_product(v1,v2)), 1e-10); 
      v2 = eigenVectors1.get_row(j); 
      EXPECT_NEAR(1.0, std::abs(dot_product(v1,v2)), 1e-10); 
      } 
    MatrixType matEst = eigenVectors.transpose()*utl::GetDiagonalMatrix(eigenValues)*eigenVectors; 
    EXPECT_NEAR_VNLMATRIX(matrixSym, matEst, 1e-10); 
    matEst = eig.V*eig.D*eig.V.transpose(); 
    EXPECT_NEAR_VNLMATRIX(matrixSym, matEst, 1e-10); 
    }
}


TEST_F(utlVNLLapack_gtest, svd_VnlMatrix)
{
  std::vector<MatrixType > matVec;
  matVec.push_back(matrixNonSym);
  for ( int i = 0; i < matVec.size(); i += 1 ) 
    {
    MatrixType U, V, matEst, mm=matVec[i]; 
    UtlMatrixType UUtl, VUtl, matEstUtl, mmUtl=utl::VnlMatrixToUtlMatrix(mm);
    vnl_vector<double> S, s1; 
    utl::NDArray<double,1> SUtl, s1Utl; 

    vnl_svd<double> svd(mm); 
    matEst = svd.U()*svd.W()*svd.V().transpose(); 
    s1=svd.W().get_diagonal();
    EXPECT_NEAR_VNLMATRIX(mm, matEst, 1e-10); 

    utl::gesdd_VnlMatrix(mm, U, S, V); 
    matEst = U*utl::GetDiagonalMatrix(S)*V.transpose(); 
    EXPECT_NEAR_VNLMATRIX(mm, matEst, 1e-10); 
    EXPECT_NEAR_VECTOR(s1, S, utl::min(S.size(),s1.size()),1e-10); 

    utl::gesvd_VnlMatrix(mm, U, S, V); 
    matEst = U*utl::GetDiagonalMatrix(S)*V.transpose(); 
    EXPECT_NEAR_VNLMATRIX(mm, matEst, 1e-10); 

    utl::gesvd_UtlMatrix(mmUtl, UUtl, SUtl, VUtl); 
    matEstUtl = UUtl*SUtl.GetDiagonalMatrix()*VUtl.GetTranspose(); 
    EXPECT_NEAR_MATRIX(mm, matEstUtl, mm.rows(), mm.cols(), 1e-10); 
    EXPECT_NEAR_VECTOR(s1, SUtl, utl::min(SUtl.Size(),s1.size()),1e-10); 
    
    utl::gesdd_UtlMatrix(mmUtl, UUtl, SUtl, VUtl); 
    matEstUtl = UUtl*SUtl.GetDiagonalMatrix()*VUtl.GetTranspose(); 
    EXPECT_NEAR_MATRIX(mm, matEstUtl, mm.rows(), mm.cols(), 1e-10); 
    EXPECT_NEAR_VECTOR(s1, SUtl, utl::min(SUtl.Size(),s1.size()),1e-10); 
    }
}

TEST_F(utlVNLLapack_gtest, InverseSymmetricVnlMatrix)
{
  std::vector<MatrixType > matVec;

  // matVec.push_back(matrixRandomSym);
  for ( int i = 0; i < matVec.size(); i += 1 ) 
    {
      {
      MatrixType inv0, inv1, mat=matVec[i];
      inv0 = utl::GetVnlSymmericMatrixPInverse(mat,1e-10);
      utl::InverseSymmericVnlMatrix(mat, inv1,1e-10); 
      EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8); 

      UtlMatrixType inv1Utl, matUtl=utl::VnlMatrixToUtlMatrix(mat);
      matUtl.InverseSymmericMatrix(inv1Utl,1e-10);
      EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8); 
      }
    }
}


TEST_F(utlVNLLapack_gtest, PInverseSymmetricVnlMatrix)
{
  std::vector<MatrixType > matVec;
  matVec.push_back(matrixSym);
  matVec.push_back(matrixRandomSym);
  for ( int i = 0; i < matVec.size(); i += 1 ) 
    {
    MatrixType inv0, inv1, mat=matVec[i];
    inv0 = utl::GetVnlSymmericMatrixPInverse(mat,1e-10);

    utl::PInverseSymmericVnlMatrix(mat, inv1,1e-10); 
    EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8); 
    utl::PInverseVnlMatrix(mat, inv1); 
    EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8); 

    UtlMatrixType inv1Utl, matUtl=utl::VnlMatrixToUtlMatrix(mat);
    matUtl.PInverseSymmericMatrix(inv1Utl,1e-10);
    EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8); 
    matUtl.PInverseMatrix(inv1Utl,1e-10);
    EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8); 
    }
}

TEST_F(utlVNLLapack_gtest, PInverseVnlMatrix)
{
  std::vector<MatrixType > matVec;
  matVec.push_back(matrixNonSym);
  matVec.push_back(matrixRandomNonSym);
  for ( int i = 0; i < matVec.size(); i += 1 ) 
    {
    MatrixType inv0, inv1, mat=matVec[i];
    inv0 = utl::GetVnlMatrixPInverse(mat,1e-10);
    utl::PInverseVnlMatrix(mat, inv1,1e-10); 
    EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8); 

    UtlMatrixType inv1Utl, matUtl=utl::VnlMatrixToUtlMatrix(mat);
    matUtl.PInverseMatrix(inv1Utl,1e-10);
    EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8); 
    }
}

TEST(utlVNLLapack, VnlMatrix_TimeCost)
{
  int M = 400, N = 300;
  int min_mn = utl::min(M, N);
  MatrixType matrixRandomNonSym(M,N);
  MatrixType matrixRandomSym(min_mn,min_mn);
  for ( int i = 0; i < M; i += 1 ) 
    for ( int j = 0; j < N; j += 1 ) 
      {
      matrixRandomNonSym(i,j) = utl::Random<double>(-2.0,2.0);
      if (i<min_mn && j<min_mn && i<=j)
        {
        matrixRandomSym(i,j) = matrixRandomNonSym(i,j);
        matrixRandomSym(j,i) = matrixRandomNonSym(i,j);
        }
      }
  
  MatrixType eigenVectors, eigenVectors1, U,V; 
  vnl_vector<double> eigenValues, eigenValues1, S, s0; 

    {
    // eigen-decomposition for symmetric matrix
    typedef vnl_symmetric_eigensystem< double >  SymEigenSystemType; 
    utl::Tic(std::cout<< "VNL vnl_symmetric_eigensystem");
    SymEigenSystemType eig (matrixRandomSym); 
    utl::Toc();
    utl::Tic(std::cout<< "lapack dsyev");
    utl::syev_VnlMatrix(matrixRandomSym, eigenValues, eigenVectors); 
    utl::Toc();
    utl::Tic(std::cout<< "lapack dsyevd");
    utl::syevd_VnlMatrix(matrixRandomSym, eigenValues1, eigenVectors1); 
    utl::Toc();
    for ( int i = 0; i < matrixRandomSym.rows(); i += 1 )  
      { 
      EXPECT_NEAR(eig.D(i,i), eigenValues(i), 1e-10); 
      EXPECT_NEAR(eig.D(i,i), eigenValues1(i), 1e-10);
      vnl_vector<double> v1, v2; 
      v1 = eig.V.get_column(i); 
      v2 = eigenVectors.get_row(i); 
      EXPECT_NEAR(1.0, std::abs(dot_product(v1,v2)), 1e-10); 
      v2 = eigenVectors1.get_row(i); 
      EXPECT_NEAR(1.0, std::abs(dot_product(v1,v2)), 1e-10); 
      } 
    }

    {
    // svd
    utl::Tic(std::cout<< "VNL vnl_svd");
    vnl_svd<double> svd(matrixRandomNonSym); 
    utl::Toc();
    s0 = svd.W().get_diagonal();
    utl::Tic(std::cout<< "lapack dgesvd");
    utl::gesvd_VnlMatrix(matrixRandomNonSym, U, S, V); 
    utl::Toc();
    EXPECT_NEAR_VECTOR(s0, S, utl::min(s0.size(), S.size()), 1e-8);
    utl::Tic(std::cout<< "lapack dgesdd");
    utl::gesdd_VnlMatrix(matrixRandomNonSym, U, S, V); 
    utl::Toc();
    EXPECT_NEAR_VECTOR(s0, S, utl::min(s0.size(), S.size()), 1e-8);
    }

    {
    // pinverse
    // make it works for InverseSymmericVnlMatrix
    for ( int i = 0; i < matrixRandomSym.rows(); ++i ) 
      matrixRandomSym(i,i) += 2;

    MatrixType inv0, inv1;
    UtlMatrixType inv1Utl, matUtl;

    utl::Tic(std::cout<< "VNL symmetric pinverse");
    inv0 = utl::GetVnlSymmericMatrixPInverse(matrixRandomSym,1e-10);
    utl::Toc();
    utl::Tic(std::cout<< "lapack symmetric pinverse");
    utl::PInverseSymmericVnlMatrix(matrixRandomSym, inv1,1e-10); 
    utl::Toc();
    EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8);
    utl::Tic(std::cout<< "lapack symmetric inverse");
    utl::InverseSymmericVnlMatrix(matrixRandomSym, inv1,1e-10); 
    utl::Toc();
    EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8);

    matUtl=utl::VnlMatrixToUtlMatrix(matrixRandomSym);
    utl::Tic(std::cout<< "lapack utl::Matrix symmetric pinverse");
    matUtl.PInverseSymmericMatrix(inv1Utl, 1e-10);
    utl::Toc();
    EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8);
    utl::Tic(std::cout<< "lapack utl::Matrix symmetric inverse");
    matUtl.InverseSymmericMatrix(inv1Utl, 1e-10);
    utl::Toc();
    EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8);

    utl::Tic(std::cout<< "VNL nonsymmetric pinverse");
    inv0 = utl::GetVnlMatrixPInverse(matrixRandomNonSym,1e-10);
    utl::Toc();
    utl::Tic(std::cout<< "lapack nonsymmetric pinverse");
    utl::PInverseVnlMatrix(matrixRandomNonSym, inv1,1e-10); 
    utl::Toc();
    EXPECT_NEAR_VNLMATRIX(inv0, inv1, 1e-8);
    
    matUtl=utl::VnlMatrixToUtlMatrix(matrixRandomNonSym);
    utl::Tic(std::cout<< "lapack utl::Matrix nonsymmetric pinverse");
    matUtl.PInverseMatrix(inv1Utl, 1e-10);
    utl::Toc();
    EXPECT_NEAR_MATRIX(inv0, inv1Utl, inv0.rows(), inv0.cols(), 1e-8);
    }

}

ITK_THREAD_RETURN_TYPE __ThreadedMethod(void* arg)
{
  const itk::MultiThreader::ThreadInfoStruct* threadInfo = static_cast<itk::MultiThreader::ThreadInfoStruct*>(arg);
  if (threadInfo)
    {
    const unsigned int threadId = threadInfo->ThreadID;
    const int numberOfThreads = threadInfo->NumberOfThreads;
    std::vector<MatrixPointer>* dataVec = static_cast<std::vector<MatrixPointer>*>(threadInfo->UserData);
    if (dataVec)
      {
      utl::PInverseSymmericVnlMatrix(*(*dataVec)[2*numberOfThreads], *(*dataVec)[threadId],1e-10); 
      utl::PInverseVnlMatrix(*(*dataVec)[2*numberOfThreads+1], *(*dataVec)[threadId+numberOfThreads],1e-10); 
      } 
    else 
      {
      std::cerr << "ERROR: UserData was not of type ThreadDataVec*" << std::endl;
      return ITK_THREAD_RETURN_VALUE;
      }
    } 
  else 
    {
    std::cerr << "ERROR: arg was not of type itk::MultiThreader::ThreadInfoStruct*" << std::endl;
    return ITK_THREAD_RETURN_VALUE;
    }
  return ITK_THREAD_RETURN_VALUE;
}


TEST(utlVNLLapack, MultiThreads_OpenMP)
{
  int M = 400, N = 300;
  int min_mn = utl::min(M, N);
  MatrixType matrixRandomNonSym(M,N);
  MatrixType matrixRandomSym(min_mn,min_mn);
  for ( int i = 0; i < M; i += 1 ) 
    for ( int j = 0; j < N; j += 1 ) 
      {
      matrixRandomNonSym(i,j) = utl::Random<double>(-2.0,2.0);
      if (i<min_mn && j<min_mn && i<=j)
        {
        matrixRandomSym(i,j) = matrixRandomNonSym(i,j);
        matrixRandomSym(j,i) = matrixRandomNonSym(i,j);
        }
      }

  // NOTE: need to build openblas with "USE_OPENMP=1"
  #pragma omp parallel
    {
    typedef vnl_symmetric_eigensystem< double >  SymEigenSystemType; 
    SymEigenSystemType eig (matrixRandomSym); 
    MatrixType eigenVectors, eigenVectors1; 
    vnl_vector<double> eigenValues, eigenValues1; 
    utl::syev_VnlMatrix(matrixRandomSym, eigenValues, eigenVectors); 
    // std::ostringstream msg;
    // utl::PrintVnlVector(eig.D.get_diagonal(), "d0", " ", msg);
    // utl::PrintVnlMatrix(eig.V, "V", " ", msg);
    // utl::PrintVnlVector(eigenValues, "eigenValues", " ", msg);
    // utl::PrintVnlMatrix(eigenVectors, "eigenVectors", " ", msg);
    // std::cout << msg.str() << std::endl << std::flush;
    for ( int i = 0; i < matrixRandomSym.rows(); i += 1 )  
      { 
      EXPECT_NEAR(eig.D(i,i), eigenValues(i), 1e-8); 
      vnl_vector<double> v1, v2;
      v1 = eig.V.get_column(i); 
      v2 = eigenVectors.get_row(i); 
      EXPECT_NEAR(1.0, std::abs(dot_product(v1,v2)), 1e-8); 
      } 

#ifdef UTL_USE_FASTLAPACK
    // NOTE: for system openblas, this may have some problems for large matrix (300x300), it works for small matrix (30x30).
    //       If openblas is built with "USE_OPENMP=1", then it works well.
    utl::syevd_VnlMatrix(matrixRandomSym, eigenValues1, eigenVectors1); 
    for ( int i = 0; i < matrixRandomSym.rows(); i += 1 )  
      { 
      EXPECT_NEAR(eig.D(i,i), eigenValues(i), 1e-8); 
      vnl_vector<double> v1, v2;
      v1 = eig.V.get_column(i); 
      EXPECT_NEAR(eig.D(i,i), eigenValues1(i), 1e-8); 
      v2 = eigenVectors1.get_row(i); 
      EXPECT_NEAR(1.0, std::abs(dot_product(v1,v2)), 1e-8); 
      }
    // std::ostringstream msg1;
    // utl::PrintVnlVector(eigenValues1, "eigenValues1", " ", msg1);
    // utl::PrintVnlMatrix(eigenVectors1, "eigenVectors1", " ", msg1);
    // std::cout << msg1.str() << std::endl << std::flush;
#endif
    }

  // multi-thread openmp
  #pragma omp parallel
    {
    MatrixType inv0, inv1;
    inv0 = utl::GetVnlSymmericMatrixPInverse(matrixRandomSym,1e-10);
    utl::PInverseSymmericVnlMatrix(matrixRandomSym, inv1,1e-10); 
    EXPECT_NEAR(0, (inv0-inv1).absolute_value_max(), 1e-8);

    inv0 = utl::GetVnlMatrixPInverse(matrixRandomNonSym,1e-10);
    utl::PInverseVnlMatrix(matrixRandomNonSym, inv1,1e-10); 
    EXPECT_NEAR(0, (inv0-inv1).absolute_value_max(), 1e-8);
    }

}

TEST(utlVNLLapack, MultiThreads_ITK)
{
  int M = 400, N = 300;
  int min_mn = utl::min(M, N);
  MatrixType matrixRandomNonSym(M,N);
  MatrixType matrixRandomSym(min_mn,min_mn);
  for ( int i = 0; i < M; i += 1 ) 
    for ( int j = 0; j < N; j += 1 ) 
      {
      matrixRandomNonSym(i,j) = utl::Random<double>(-2.0,2.0);
      if (i<min_mn && j<min_mn && i<=j)
        {
        matrixRandomSym(i,j) = matrixRandomNonSym(i,j);
        matrixRandomSym(j,i) = matrixRandomNonSym(i,j);
        }
      }

  // multi-thread itk
  // NOTE: need to build openblas with "USE_THREAD=0" or set "export OMP_NUM_THREADS=1" and export OMP_NUM_THREADS=1"
    {
    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    int numberOfThreads = 3;
    std::vector<MatrixPointer> matVec(2*numberOfThreads+2);
    for ( int i = 0; i < matVec.size(); i += 1 ) 
      matVec[i] = MatrixPointer(new MatrixType());
    *matVec[2*numberOfThreads] = matrixRandomSym;
    *matVec[2*numberOfThreads+1] = matrixRandomNonSym;
    MatrixType inv0, inv1;
    inv0 = utl::GetVnlSymmericMatrixPInverse(matrixRandomSym,1e-10);
    inv1 = utl::GetVnlMatrixPInverse(matrixRandomNonSym,1e-10);

    threader->SetGlobalMaximumNumberOfThreads(numberOfThreads+1);
    threader->SetNumberOfThreads(numberOfThreads);
    threader->SetSingleMethod(__ThreadedMethod, &matVec);
    threader->SingleMethodExecute();
    for ( int ii = 0; ii < numberOfThreads; ii += 1 ) 
      {
      EXPECT_NEAR(0.0, (inv0-*(matVec[ii])).absolute_value_max(), 1e-8 );
      EXPECT_NEAR(0.0, (inv1-*(matVec[ii+numberOfThreads])).absolute_value_max(), 1e-8 );
      }
    }
}

