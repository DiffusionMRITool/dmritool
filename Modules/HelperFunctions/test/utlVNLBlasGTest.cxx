/**
 *       @file  utlVNLBlasGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "07-19-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

// #include "utl.h"
#include "utlGTest.h"
#include "utlVNLBlas.h"
#include "utlSTDHeaders.h"
#include "itkMultiThreader.h"
#include "utlMatrix.h"
#include "utlVNLIO.h"

typedef vnl_matrix<double>         MatrixType;
typedef utl::NDArray<double,2>     UtlMatrixType;
typedef utl_shared_ptr<MatrixType> MatrixPointer;
typedef utl_shared_ptr<UtlMatrixType> UtlMatrixPointer;

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const T val1, const T val2 )
{
  vnl_vector<double> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

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

class utlVNLBlas_gtest : public testing::Test 
{
protected:

  virtual void SetUp() 
    {
    const double array_1 [] =
      {
      1,  1, 
      -1, 2,
      2,  1
      };
    const double array_2 [] =
      {
      1,  1, 3, 4 ,
      -1, 2, 5, 6
      };
    vec1.set_size(2);
    vec1[0]=1; vec1[1]=4;

    mat1 = vnl_matrix<double>(array_1, 3, 2);
    mat2 = vnl_matrix<double>(array_2, 2, 4);
    }

  vnl_matrix<double> mat1;
  vnl_matrix<double> mat2;
  vnl_vector<double> vec1;
};

TEST_F(utlVNLBlas_gtest, MatrixTimesMatrix)
{
  vnl_matrix<double> mat_mul, mat_mulBlas;
  mat_mul = mat1*mat2;
  utl::ProductVnlMM(mat1, mat2, mat_mulBlas);
  EXPECT_NEAR_VNLMATRIX( mat_mul, mat_mulBlas, 1e-15 );

  vnl_matrix<double> mat2t = mat2.transpose();
  mat_mul = mat1*mat2t.transpose();
  utl::ProductVnlMMt(mat1, mat2t, mat_mulBlas);
  EXPECT_NEAR_VNLMATRIX( mat_mul, mat_mulBlas, 1e-15 );
}

TEST_F(utlVNLBlas_gtest, MatrixTimesVector)
{
  vnl_vector<double> vec_mul, vec_mulBlas;
  vec_mul = mat1*vec1;
  utl::ProductVnlMv(mat1, vec1, vec_mulBlas);
  EXPECT_NEAR_VNLVECTOR( vec_mul, vec_mulBlas, 1e-15 );

  vnl_matrix<double> mat1t = mat1.transpose();
  vec_mul = mat1t.transpose()*vec1;
  utl::ProductVnlMtv(mat1t, vec1, vec_mulBlas);
  EXPECT_NEAR_VNLVECTOR( vec_mul, vec_mulBlas, 1e-15 );
}

TEST_F(utlVNLBlas_gtest, VectorTimesMatrix)
{
  vnl_vector<double> vec_mul, vec_mulBlas;
  vec_mul = vec1*mat2;
  utl::ProductVnlvM(vec1, mat2, vec_mulBlas);
  EXPECT_NEAR_VNLVECTOR( vec_mul, vec_mulBlas, 1e-15 );

  vnl_matrix<double> mat2t = mat2.transpose();
  vec_mul = vec1*mat2t.transpose();
  utl::ProductVnlvMt(vec1, mat2t, vec_mulBlas);
  EXPECT_NEAR_VNLVECTOR( vec_mul, vec_mulBlas, 1e-15 );
}

TEST(utlVNLBlas, MatrixTimesMatrix_TimeCost)
{
  vnl_matrix<double> mat1 = __GenerateRandomMatrix(514, 180, -2.0,2.0);
  vnl_matrix<double> mat2 = __GenerateRandomMatrix(180, 254, -2.0,2.0);
  UtlMatrixType mat1Utl = utl::VnlMatrixToUtlMatrix(mat1);
  UtlMatrixType mat2Utl = utl::VnlMatrixToUtlMatrix(mat2);

  int N = 20;
    {
    utl::Tic(std::cout<<"vnl matrix multiplization time: \n");
    for ( int i = 0; i < N; i += 1 ) 
      {
      vnl_matrix<double> tmp;
      tmp = mat1*mat2;
      }
    utl::Toc();
    }
    {
    utl::Tic(std::cout<<"ProductVnlMM time:\n");
    for ( int i = 0; i < N; i += 1 ) 
      {
      vnl_matrix<double> tmp;
      utl::ProductVnlMM(mat1, mat2, tmp);
      }
    utl::Toc();
    utl::Tic(std::cout<<"ProductUtlMM time:\n");
    for ( int i = 0; i < N; i += 1 ) 
      {
      UtlMatrixType tmp;
      utl::ProductUtlMM(mat1Utl, mat2Utl, tmp);
      }
    utl::Toc();
    utl::Tic(std::cout<<"utl::Matrix product time:\n");
    for ( int i = 0; i < N; i += 1 ) 
      {
      UtlMatrixType tmp;
      tmp = mat1Utl*mat2Utl;
      }
    utl::Toc();
    }
}

TEST(utlVNLBlas, MultiThreads_OpenMP)
{
  vnl_matrix<double> mat1 = __GenerateRandomMatrix(514, 180, -2.0,2.0);
  vnl_matrix<double> mat2 = __GenerateRandomMatrix(180, 254, -2.0,2.0);
  UtlMatrixType mat1Utl = utl::VnlMatrixToUtlMatrix(mat1);
  UtlMatrixType mat2Utl = utl::VnlMatrixToUtlMatrix(mat2);
      
  vnl_matrix<double> product;
  utl::ProductVnlMM(mat1, mat2, product);
  UtlMatrixType productUtl = mat1Utl*mat2Utl;
  EXPECT_NEAR_MATRIX(product, productUtl, product.rows(), product.cols(), 1e-10);

#pragma omp parallel
    {
    vnl_matrix<double> tmp;
    utl::ProductVnlMM(mat1, mat2, tmp);
    EXPECT_NEAR(0.0, (tmp-product).absolute_value_max(), 1e-8 );
    UtlMatrixType tmpUtl = mat1Utl* mat2Utl;
    EXPECT_NEAR(0.0, utl::ToMatrix<double>(tmpUtl-productUtl)->AbsoluteMaxValue(), 1e-8 );
    }
}


ITK_THREAD_RETURN_TYPE __ThreadedMethod2(void* arg)
{
  const itk::MultiThreader::ThreadInfoStruct* threadInfo = static_cast<itk::MultiThreader::ThreadInfoStruct*>(arg);
  if (threadInfo)
    {
    const unsigned int threadId = threadInfo->ThreadID;
    const int numberOfThreads = threadInfo->NumberOfThreads;
    std::vector<UtlMatrixPointer>* dataVec = static_cast<std::vector<UtlMatrixPointer>*>(threadInfo->UserData);
    if (dataVec)
      {
      *(*dataVec)[threadId] =  *(*dataVec)[numberOfThreads] * (*(*dataVec)[numberOfThreads+1]);
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
      utl::ProductVnlMM(*(*dataVec)[numberOfThreads], *(*dataVec)[numberOfThreads+1], *(*dataVec)[threadId]);
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

TEST(utlVNLBlas, MultiThreads_ITK)
{
  vnl_matrix<double> mat1 = __GenerateRandomMatrix(514, 180, -2.0,2.0);
  vnl_matrix<double> mat2 = __GenerateRandomMatrix(180, 254, -2.0,2.0);
  UtlMatrixType mat1Utl = utl::VnlMatrixToUtlMatrix(mat1);
  UtlMatrixType mat2Utl = utl::VnlMatrixToUtlMatrix(mat2);

  vnl_matrix<double> product;
  utl::ProductVnlMM(mat1, mat2, product);
  UtlMatrixType productUtl = mat1Utl*mat2Utl;

  // multi-thread itk
    {
    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    int numberOfThreads = 3;
    std::vector<MatrixPointer> matVec(numberOfThreads+3);
    for ( int i = 0; i < matVec.size(); i += 1 ) 
      matVec[i] = MatrixPointer(new MatrixType());
    *matVec[numberOfThreads] = mat1;
    *matVec[numberOfThreads+1] = mat2;

    threader->SetGlobalMaximumNumberOfThreads(numberOfThreads+1);
    threader->SetNumberOfThreads(numberOfThreads);
    threader->SetSingleMethod(__ThreadedMethod, &matVec);
    threader->SingleMethodExecute();
    for ( int ii = 0; ii < numberOfThreads; ii += 1 ) 
      EXPECT_NEAR(0.0, (product-*(matVec[ii])).absolute_value_max(), 1e-8 );
    }
    {
    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    int numberOfThreads = 3;
    std::vector<UtlMatrixPointer> matVec(numberOfThreads+3);
    for ( int i = 0; i < matVec.size(); i += 1 ) 
      matVec[i] = UtlMatrixPointer(new UtlMatrixType());
    *matVec[numberOfThreads] = mat1Utl;
    *matVec[numberOfThreads+1] = mat2Utl;

    threader->SetGlobalMaximumNumberOfThreads(numberOfThreads+1);
    threader->SetNumberOfThreads(numberOfThreads);
    threader->SetSingleMethod(__ThreadedMethod2, &matVec);
    threader->SingleMethodExecute();
    for ( int ii = 0; ii < numberOfThreads; ii += 1 ) 
      EXPECT_NEAR(0.0, utl::ToMatrix<double>(productUtl-*(matVec[ii]))->AbsoluteMaxValue(), 1e-8 );
    }
}


TEST(utlVNLBlas, cblas_nrm2)
{
  int N = 100;
  vnl_vector<double> v1 = __GenerateRandomVector<double>(N, -2.0,2.0);
  EXPECT_NEAR(v1.magnitude(), utl::cblas_nrm2(v1.size(), v1.data_block(), 1), 1e-10*v1.magnitude());
}

TEST(utlVNLBlas, cblas_asum)
{
  int N = 100;
  vnl_vector<double> v1 = __GenerateRandomVector<double>(N, -2.0,2.0);
  double sum=0;
  for ( int i = 0; i < N; ++i ) 
    sum += std::fabs(v1[i]);
  EXPECT_NEAR(sum, utl::cblas_asum<double>(v1.size(), v1.data_block(), 1), sum*1e-10);
}


TEST(utlVNLBlas, InnerProductvv)
{
  int N = 100;
  vnl_vector<double> v1 = __GenerateRandomVector<double>(N, -2.0,2.0);
  vnl_vector<double> v2 = __GenerateRandomVector<double>(N, -2.0,2.0);
  EXPECT_NEAR(dot_product(v1,v2), utl::InnerProduct(v1,v2), 1e-10*std::fabs(dot_product(v1,v2)));
  
  int K = 1e5;
    {
    double sum=0;
    utl::Tic(std::cout<<"vnl dot_product time: \n");
    for ( int i = 0; i < K; ++i ) 
      sum += dot_product(v1,v2);
    utl::Toc();
    
    utl::Tic(std::cout<<"blas dot time: \n");
    for ( int i = 0; i < K; ++i ) 
      sum += utl::InnerProduct(v1,v2);
    utl::Toc();
    }
}

TEST(utlVNLBlas, OuterProduct)
{
  int N = 100;
  vnl_vector<double> v1 = __GenerateRandomVector<double>(N, -2.0,2.0);
  vnl_vector<double> v2 = __GenerateRandomVector<double>(N, -2.0,2.0);
    {
    vnl_matrix<double> mat0 = outer_product(v1,v2);
    vnl_matrix<double> mat1;
    utl::OuterProduct(v1,v2,mat1);
    EXPECT_NEAR((mat0-mat1).absolute_value_max(), 0.0, 1e-10);
    }
    {
    vnl_matrix<double> mat0 = outer_product(v1,v1);
    vnl_matrix<double> mat1;
    utl::OuterProduct(v1,mat1);
    EXPECT_NEAR((mat0-mat1).absolute_value_max(), 0.0, 1e-10);
    }

  int K = 1e4;
    {
    vnl_matrix<double> mat1, mat2;
    utl::Tic(std::cout<<"vnl outer_product time: \n");
    for ( int i = 0; i < K; ++i ) 
      mat1= outer_product(v1,v2);
    utl::Toc();
    
    utl::Tic(std::cout<<"blas outer product time: \n");
    for ( int i = 0; i < K; ++i ) 
      utl::OuterProduct(v1,v2, mat2);
    utl::Toc();
    }
}

TEST(utlVNLBlas, cblas_syrk)
{
  MatrixType mat1 = __GenerateRandomMatrix<double>(300, 400, -2.0,2.0);
    {
    MatrixType XXt_0 = mat1*mat1.transpose(), XXt_1;
    utl::ProductVnlXXt<double>(mat1, XXt_1);
    EXPECT_NEAR(0.0, (XXt_0-XXt_1).absolute_value_max(), 1e-8 );
    }
    {
    MatrixType XtX_0 = mat1.transpose()*mat1, XtX_1;
    utl::ProductVnlXtX<double>(mat1, XtX_1);
    EXPECT_NEAR(0.0, (XtX_0-XtX_1).absolute_value_max(), 1e-8 );
    }
    {
    int K=20;
    utl::Tic(std::cout<<"vnl matrix multiplization time: \n");
    vnl_matrix<double> tmp;
    for ( int i = 0; i < K; i += 1 ) 
      tmp = mat1*mat1.transpose();
    utl::Toc();
    utl::Tic(std::cout<<"cblas_syrk time: \n");
    for ( int i = 0; i < K; i += 1 ) 
      utl::ProductVnlXXt<double>(mat1, tmp);
    utl::Toc();
    utl::Tic(std::cout<<"ProductVnlMMt time: \n");
    for ( int i = 0; i < K; i += 1 ) 
      utl::ProductVnlMMt<double>(mat1, mat1, tmp);
    utl::Toc();
    }
}

TEST(utlVNLBlas, GetRow)
{
  vnl_matrix<double> mat1 = __GenerateRandomMatrix(300, 300, -2.0,2.0);
  vnl_vector<double> vec0 = mat1.get_row(100);
  vnl_vector<double> vec1;
  utl::GetRow(mat1, 100, vec1);
  EXPECT_NEAR((vec0-vec1).max_value(), 0.0, 1e-10);
  EXPECT_NEAR((vec0-vec1).min_value(), 0.0, 1e-10);

  int N = 1e5;
    {
    std::vector<int> indexVec;
    for ( int i = 0; i < N; ++i ) 
      indexVec.push_back(utl::RandomInt(0,200));

    utl::Tic(std::cout<<"vnl get_row time: \n");
    for ( int i = 0; i < N; ++i ) 
      {
      vec0 = mat1.get_row(indexVec[i]);
      }
    utl::Toc();
    utl::Tic(std::cout<<"blas GetRow time: \n");
    for ( int i = 0; i < N; ++i ) 
      {
      utl::GetRow(mat1, indexVec[i], vec1);
      }
    utl::Toc();
    }
}

TEST(utlVNLBlas, GetColumn)
{
  vnl_matrix<double> mat1 = __GenerateRandomMatrix(300, 300, -2.0,2.0);
  vnl_vector<double> vec0 = mat1.get_column(100);
  vnl_vector<double> vec1;
  utl::GetColumn(mat1, 100, vec1);
  EXPECT_NEAR((vec0-vec1).max_value(), 0.0, 1e-10);
  EXPECT_NEAR((vec0-vec1).min_value(), 0.0, 1e-10);

  int N = 1e5;
    {
    std::vector<int> indexVec;
    for ( int i = 0; i < N; ++i ) 
      indexVec.push_back(utl::RandomInt(0,200));

    utl::Tic(std::cout<<"vnl get_column time: \n");
    for ( int i = 0; i < N; ++i ) 
      {
      vec0 = mat1.get_column(indexVec[i]);
      }
    utl::Toc();
    utl::Tic(std::cout<<"blas GetColumn time: \n");
    for ( int i = 0; i < N; ++i ) 
      {
      utl::GetColumn(mat1, indexVec[i], vec1);
      }
    utl::Toc();
    }
}


TEST(utlVNLBlas, MatrixCopy)
{
  MatrixType mat0 = __GenerateRandomMatrix<double>(300,200,-2.0,2.0);
  MatrixType mat1, mat2;
  UtlMatrixType mat0Utl = utl::VnlMatrixToUtlMatrix(mat0), mat1Utl, mat2Utl;

    {
    // copy
    mat1=mat0;
    utl::MatrixCopy(mat0, mat2, 1.0, 'N');
    EXPECT_NEAR(0.0, (mat1-mat2).absolute_value_max(), 1e-8 );

    mat1Utl=mat0Utl;
    EXPECT_NEAR(0.0, utl::ToMatrix<double>(mat1Utl-mat0Utl)->AbsoluteMaxValue(), 1e-8 );

    int K = 1e3;
    // int K = 1;
    utl::Tic(std::cout<< "VNL matrix copy time:");
    for ( int i = 0; i < K; ++i ) 
      mat1=mat0;
    utl::Toc();
    utl::Tic(std::cout<< "utl::MatrixCopy matrix copy time:");
    for ( int i = 0; i < K; ++i ) 
      utl::MatrixCopy(mat0, mat2, 1.0, 'N');
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix matrix copy time:");
    // std::cout << "a1" << std::endl << std::flush;
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = mat0Utl;
    // std::cout << "a2" << std::endl << std::flush;
    utl::Toc();
    }
  
    {
    // scale
    double alpha = utl::Random<double>(-2.0,2.0);
    mat1 = alpha*mat0;
    utl::MatrixCopy(mat0, mat2, alpha, 'N');
    EXPECT_NEAR(0.0, (mat1-mat2).absolute_value_max(), 1e-8 );

    int K = 1e3;
    utl::Tic(std::cout<< "VNL matrix scale time:");
    for ( int i = 0; i < K; ++i ) 
      mat1=mat0*alpha;
    utl::Toc();
    utl::Tic(std::cout<< "utl::MatrixCopy matrix scale time:");
    for ( int i = 0; i < K; ++i ) 
      utl::MatrixCopy(mat0, mat2, alpha, 'N');
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix matrix x scale time:");
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = mat0Utl%alpha;
    utl::Toc();
    mat2Utl = mat0Utl;
    utl::Tic(std::cout<< "utl::Matrix matrix x scale time2:");
    for ( int i = 0; i < K; ++i ) 
      {
      UtlMatrixType tmp(mat0Utl);
      tmp.Scale(alpha);
      mat2Utl.Swap(tmp);
      }
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix scale x matrix time:");
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = alpha%mat0Utl;
    utl::Toc();
    }

    {
    // addtion
    double alpha = utl::Random<double>(-2.0,2.0);
    mat1 = alpha+mat0;

    int K = 1e3;
    utl::Tic(std::cout<< "VNL matrix + scale time:");
    for ( int i = 0; i < K; ++i ) 
      mat1=mat0+alpha;
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix matrix + scale time:");
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = mat0Utl+alpha;
    utl::Toc();
    mat2Utl = mat0Utl;
    utl::Tic(std::cout<< "utl::Matrix matrix + scale time2:");
    for ( int i = 0; i < K; ++i ) 
      {
      UtlMatrixType tmp(mat0Utl);
      tmp+=alpha;
      mat2Utl.Swap(tmp);
      }
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix scale + matrix time:");
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = alpha+mat0Utl;
    utl::Toc();
    }

    {
    // division
    double alpha = utl::Random<double>(-2.0,2.0);
    mat1 = mat0/alpha;

    int K = 1e3;
    utl::Tic(std::cout<< "VNL matrix / scale time:");
    for ( int i = 0; i < K; ++i ) 
      mat1=mat0/alpha;
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix matrix / scale time:");
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = mat0Utl/alpha;
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix matrix * 1.0/scale time:");
    for ( int i = 0; i < K; ++i ) 
      {
      double alpha_div=1.0/alpha;
      mat2Utl = mat0Utl%alpha_div;
      }
    utl::Toc();
    mat2Utl = mat0Utl;
    utl::Tic(std::cout<< "utl::Matrix matrix / scale time2:");
    for ( int i = 0; i < K; ++i ) 
      {
      UtlMatrixType tmp(mat0Utl);
      tmp/=alpha;
      mat2Utl.Swap(tmp);
      }
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix scale / matrix time:");
    for ( int i = 0; i < K; ++i ) 
      mat2Utl = alpha/mat0Utl;
    utl::Toc();
    }


    
    {
    // transpose
    mat1 = mat0.transpose();
    utl::MatrixCopy(mat0, mat2, 1.0, 'T');
    EXPECT_NEAR(0.0, (mat1-mat2).absolute_value_max(), 1e-8 );

    int K = 1e3;
    utl::Tic(std::cout<< "VNL matrix transpose time:");
    for ( int i = 0; i < K; ++i ) 
      mat1=mat0.transpose();
    utl::Toc();
    utl::Tic(std::cout<< "utl::MatrixCopy matrix transpose time:");
    for ( int i = 0; i < K; ++i ) 
      utl::MatrixCopy(mat0, mat2, 1.0, 'T');
    utl::Toc();
    utl::Tic(std::cout<< "utl::Matrix matrix transpose time:");
    for ( int i = 0; i < K; ++i ) 
       mat0Utl.GetTranspose(mat1Utl);
    utl::Toc();
    }
}
