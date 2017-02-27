/**
 *       @file  utlCoreMKLGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlGTest.h"

#include "utlCore.h"
#include "utlCoreMKL.h"
#include "utlVNL.h"


typedef vnl_vector<double>         VectorType;

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const T val1, const T val2 )
{
  vnl_vector<double> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

TEST(utlCoreMKL, vAdd)
{
  int N = 1000;
  VectorType vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  VectorType vec2 = __GenerateRandomVector<double>(N, -2.0, 2.0);

  VectorType result0 = vec1+vec2, result1(N);
  utl::vAdd(N, vec1.data_block(), vec2.data_block(), result1.data_block());
  EXPECT_NEAR_VNLVECTOR(result0, result1, 1e-10);

    {
    result0=vec2;
    result0 += vec1;
    result1 = vec2;
    utl::vAdd(N, result1.data_block(), vec1.data_block(), result1.data_block());
    EXPECT_NEAR_VNLVECTOR(result0, result1, 1e-10);
    }
  
  int K = 1e4;
    {
    utl::Tic(std::cout<<"vnl add");
    for ( int i = 0; i < K; ++i ) 
      result0 = vec1 + vec2;
    utl::Toc();
    utl::Tic(std::cout<<"mkl vadd");
    for ( int i = 0; i < K; ++i ) 
      utl::vAdd(N, vec1.data_block(), vec2.data_block(), result1.data_block());
    utl::Toc();
    }
}

TEST(utlCoreMKL, vMul)
{
  int N = 1000;
  VectorType vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  VectorType vec2 = __GenerateRandomVector<double>(N, -2.0, 2.0);

  VectorType result0 = element_product(vec1, vec2), result1(N);
  utl::vMul(N, vec1.data_block(), vec2.data_block(), result1.data_block());
  EXPECT_NEAR_VNLVECTOR(result0, result1, 1e-10);
  
  int K = 1e4;
    {
    utl::Tic(std::cout<<"vnl mul");
    for ( int i = 0; i < K; ++i ) 
      result0 = element_product(vec1, vec2);
    utl::Toc();
    utl::Tic(std::cout<<"mkl vMul");
    for ( int i = 0; i < K; ++i ) 
      utl::vMul(N, vec1.data_block(), vec2.data_block(), result1.data_block());
    utl::Toc();
    }
}

