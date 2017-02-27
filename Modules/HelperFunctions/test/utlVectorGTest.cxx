/**
 *       @file  utlVectorGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-22-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlGTest.h"
#include "utlNDArray.h"
#include "utlVNL.h"

typedef utl::NDArray<double,1> UtlVectorType;

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const T val1, const T val2 )
{
  vnl_vector<double> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

TEST(utlVector, Constructors)
{
  int N =4;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
    {
    // data constructor
    UtlVectorType vec1(vec0.data_block(), vec0.size()), vec2;
    EXPECT_NEAR_VECTOR(vec0, vec1, N, 1e-10);
    vec2 = utl::VnlVectorToStdVector(vec0);
    EXPECT_NEAR_VECTOR(vec0, vec2, N, 1e-10);

    UtlVectorType vec00(3);
    vec00[0]=1.0, vec00[1]=2.0, vec00[2]=3.0;
    vec2={1,2,3};
    UtlVectorType vec3({1,2,3});
    UtlVectorType vec4{1,2,3};
    EXPECT_NEAR_VECTOR(vec00, vec2, 3, 1e-10);
    EXPECT_NEAR_VECTOR(vec00, vec3, 3, 1e-10);
    EXPECT_NEAR_VECTOR(vec00, vec4, 3, 1e-10);
    }
    {
    // change size
    UtlVectorType vec1(N);
    EXPECT_EQ(N, vec1.Size());
    vec1.Clear();
    EXPECT_EQ(0, vec1.Size());
    vec1.ReSize(4);
    EXPECT_EQ(4, vec1.Size());
    }
    {
    // data copy
    UtlVectorType vec1(vec0.data_block(), vec0.size()), vec2;
    vec2 = vec1;
    EXPECT_NEAR_VECTOR(vec0, vec2, N, 1e-10);
    UtlVectorType vec3(vec1);
    EXPECT_NEAR_VECTOR(vec0, vec3, N, 1e-10);
    }
    {
    // data reference 
    UtlVectorType vec1(30);
    vec1.SetData(vec0.data_block(), vec0.size());
    EXPECT_NEAR_VECTOR(vec0, vec1, N, 1e-10);
    EXPECT_EQ((void*)vec0.data_block(), (void*)vec1.GetData());
    }
    {
    // copy in and out
    UtlVectorType vec1(N);
    vec1=10;
    for ( int i = 0; i < N; ++i ) 
      EXPECT_NEAR(10, vec1[i], 1e-10);
    double* data = vec0.data_block();
    vec1.CopyIn(data,N);
    EXPECT_NEAR_VECTOR(data, vec1, N, 1e-10);
    vec1[2] = 100; 
    vec1.CopyOut(data,N);
    EXPECT_NEAR_VECTOR(data, vec1, N, 1e-10);
    }
    {
    // different types 
    UtlVectorType vec1(vec0.data_block(), vec0.size());
    utl::NDArray<int,1> vec2(vec1), vec3;
    vec3 = vec1;
    UtlVectorType vec4(vec2), vec5;
    vec5 = vec3;
    for ( int i = 0; i < vec0.size(); ++i ) 
      {
      EXPECT_EQ(vec2[i], (int)vec1[i]);
      EXPECT_EQ(vec3[i], (int)vec1[i]);
      EXPECT_EQ(vec4[i], (int)vec1[i]);
      EXPECT_EQ(vec5[i], (int)vec1[i]);
      }
    }
}

TEST(utlVector_DeathTest, AssignNDArrayDifferentDimension)
{
  utl::NDArray<double, 1> vec(3);
  utl::NDArray<double, 2> mat(2,3);
  EXPECT_DEATH(vec=mat, "");
  EXPECT_DEATH(mat=vec, "");
}


TEST(utlVector, MaxMin)
{
  int N =10;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  UtlVectorType vec1(vec0.data_block(), vec0.size());
    {
    EXPECT_EQ(vec0.arg_max(), vec1.ArgMax());
    EXPECT_EQ(vec0.arg_min(), vec1.ArgMin());
    EXPECT_EQ(vec0.max_value(), vec1.MaxValue());
    EXPECT_EQ(vec0.min_value(), vec1.MinValue());
    }
}

TEST(utlVector, Rotate)
{
  UtlVectorType v(3), axis(3), vRotate(3);

    {
    v[0]=1, v[1]=1, v[2]=0;
    axis[0]=0, axis[1]=0, axis[2]=1;
    vRotate=v;
    vRotate.RotateAxis(M_PI/4, axis);
    EXPECT_NEAR(0, vRotate[0], 1e-10);
    EXPECT_NEAR(std::sqrt(2), vRotate[1], 1e-10);
    EXPECT_NEAR(0, vRotate[2], 1e-10);
    }
    
    {
    v[0]=1, v[1]=1, v[2]=1;
    axis[0]=0, axis[1]=0, axis[2]=1;
    vRotate=v;
    vRotate.RotateAxis(M_PI/4, axis);
    EXPECT_NEAR(0, vRotate[0], 1e-10);
    EXPECT_NEAR(std::sqrt(2), vRotate[1], 1e-10);
    EXPECT_NEAR(1, vRotate[2], 1e-10);
    }

}

TEST(utlVector, Norms)
{
  int N =10;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  UtlVectorType vec1(vec0.data_block(), vec0.size());
    {
    EXPECT_NEAR(vec0.two_norm(), vec1.GetTwoNorm(), vec0.two_norm()*1e-10);
    EXPECT_NEAR(vec0.squared_magnitude(), vec1.GetSquaredTwoNorm(), vec0.squared_magnitude()*1e-10);
    EXPECT_NEAR(vec0.inf_norm(), vec1.GetInfNorm(), vec0.inf_norm()*1e-10);
    EXPECT_NEAR(vec0.one_norm(), vec1.GetOneNorm(), vec0.one_norm()*1e-10);
    EXPECT_NEAR(vec0.rms(), vec1.GetRootMeanSquares(), vec0.one_norm()*1e-10);
    }
}

#define __utlVector_Operators(vec0, funcName, funcReal) \
  do \
    {\
    UtlVectorType vec1(vec0.data_block(), vec0.size()); \
    vec1.funcName(); \
    for ( int i = 0; i < N; ++i ) \
      EXPECT_NEAR(funcReal(vec0[i]), vec1[i], 1e-10);\
    } while ( 0 )


TEST(utlVector, Operators)
{
  int N =10;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0), vec2;
  __utlVector_Operators(vec0, ElementAbsolute, std::fabs);
  __utlVector_Operators(vec0, ElementInverse, 1.0/);
    {
    vnl_vector<double> vv = vec1 +5;
    __utlVector_Operators(vv, ElementSqrt, std::sqrt);
    }

  UtlVectorType d0(vec0.data_block(), vec0.size()); 
  UtlVectorType d1(vec1.data_block(), vec1.size()), d2; 
    {
    // +
    vec2 = vec0+vec1;
    d2 = d0+d1;
    EXPECT_NEAR_VECTOR(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // -
    vec2=vec0; vec2-=vec1;
    d2=d0; d2-=d1;
    EXPECT_NEAR_VECTOR(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // *
    vec2=vec0; vec2=element_product(vec2, vec1);
    d2=d0; d2%=d1;
    EXPECT_NEAR_VECTOR(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // /
    vnl_vector<double> vv = vec1 +5;
    vec2=vec0; vec2=element_quotient(vec2, vv);
    UtlVectorType dvv(vv.data_block(), vv.size());
    d2=d0; d2/=dvv; 
    EXPECT_NEAR_VECTOR(vec2.data_block(), d2.GetData(), N, 1e-10);
    }

    {
    // dot
    EXPECT_NEAR(dot_product(vec0, vec1), utl::InnerProduct(d0,d1), std::fabs(dot_product(vec0, vec1))*1e-10);
    }

    {
    // cross product
    UtlVectorType v1(3), v2(3), v3(3), v0;
    v1[0]=1, v1[1]=0, v1[2]=0;
    v2[0]=0, v2[1]=1, v2[2]=0;
    v3[0]=0, v3[1]=0, v3[2]=1;
    utl::CrossProduct(v1, v2, v0);
    EXPECT_NEAR_VECTOR(v3, v0, 3, 1e-10);
    }

    {
    // swap
    d0.Swap(d1);
    EXPECT_NEAR_VECTOR(vec0.data_block(), d1.GetData(), N, 1e-10);
    EXPECT_NEAR_VECTOR(vec1.data_block(), d0.GetData(), N, 1e-10);
    d0.Swap(d1);
    }

    {
    // flip
    d2 = d1;
    d2.Flip();
    for ( int i = 0; i < N/2; ++i ) 
      EXPECT_NEAR(d1[i], d2[N-1-i], 1e-10);
    }

    {
    // friend op
    d2 = d1%3.0;
    vec2 = vec1*3.0;
    EXPECT_NEAR_VECTOR(d2.GetData(), vec2.data_block(), N, 1e-10);

    d2 = 3.0%d1;
    vec2 = 3.0*vec1;
    EXPECT_NEAR_VECTOR(d2.GetData(), vec2.data_block(), N, 1e-10);

    d2 = 3.0+d1;
    vec2 = vnl_vector<double>(vec1.size(), 3)+vec1;
    EXPECT_NEAR_VECTOR(d2.GetData(), vec2.data_block(), N, 1e-10);
    
    d2 = 3.0-d1;
    vec2 = vnl_vector<double>(vec1.size(), 3)-vec1;
    EXPECT_NEAR_VECTOR(d2.GetData(), vec2.data_block(), N, 1e-10);
    
    d2 = 3.0/(3.0+d1);
    for ( int i = 0; i < N; ++i ) 
      EXPECT_NEAR(3.0/(3.0+vec1[i]), d2[i], 1e-10);
    
    d2 = (3.0-d1)%2.0;
    for ( int i = 0; i < N; ++i ) 
      EXPECT_NEAR((3.0-vec1[i])*2.0, d2[i], 1e-10);
    }

    {
    d2 = 3.0%d0-d1;
    UtlVectorType d3(d1);
    d3.ElementAxpby(d0.GetData(), 3.0, -1.0);
    EXPECT_NEAR_VECTOR(d2.GetData(), d3.GetData(), N, 1e-10);
    }
}

TEST(utlVector, Operators_TimeCost)
{
  int N = 200;
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0), vec2(N);
  vnl_vector<double> vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  
  UtlVectorType d0(vec0.data_block(), vec0.size()); 
  UtlVectorType d1(vec1.data_block(), vec1.size()), d2(N); 


  int Iter=30000;
    {
    utl::Tic(std::cout<<"expression raw time: \n");
    for ( int i = 0; i < Iter; i += 1 ) 
      {
      for ( int j = 0; j < N; ++j ) 
        {
        d2[j] = 3.0*d0[j]-i*d1[j];
        }
      }
    utl::Toc();
    }
    {
    utl::Tic(std::cout<<"utl vector expression time: \n");
    for ( int i = 0; i < Iter; i += 1 ) 
      {
      d2 = 3.0%d0-i*1.0%d1;
      }
    utl::Toc();
    }
    {
    utl::Tic(std::cout<<"vnl vector expression time: \n");
    for ( int i = 0; i < Iter; i += 1 ) 
      {
      vec2 = 3.0*vec0 -i*1.0*vec1;
      }
    utl::Toc();
    }
}


