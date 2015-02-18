/**
 *       @file  utlNDArrayGTest.cxx
 *      @brief  gtest for utl::NDArray
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "utlGTest.h"

#include "utlNDArray.h"
#include "utlVNL.h"
#include "utlFunctors.h"

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const T val1, const T val2 )
{
  vnl_vector<double> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
  return v1;
}

template <unsigned Dim>
void 
test_utlNDArray_Constructors ( unsigned shape[Dim] )
{
  int N = 1;
  for ( int i = 0; i < Dim; ++i ) 
    N *= shape[i];
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  typedef utl::NDArray<double,Dim> ArrayType;
    {
    // data constructor
    ArrayType vec1(vec0.data_block(), shape), vec2(shape);
    EXPECT_NEAR_VECTOR(vec0, vec1, N, 1e-10);
    vec2 = utl::VnlVectorToStdVector(vec0);;
    EXPECT_NEAR_VECTOR(vec0, vec2, N, 1e-10);
    }
    {
    // change size
    ArrayType vec1(shape);
    EXPECT_EQ(N, vec1.Size());
    vec1.Clear();
    EXPECT_EQ(0, vec1.Size());
    vec1.ReSize(shape);
    EXPECT_EQ(N, vec1.Size());
    }
    {
    // data copy
    ArrayType vec1(vec0.data_block(), shape), vec2;
    vec2 = vec1;
    EXPECT_NEAR_VECTOR(vec0, vec2, N, 1e-10);
    ArrayType vec3(vec1);
    EXPECT_NEAR_VECTOR(vec0, vec3, N, 1e-10);
    }
    {
    // data reference 
    ArrayType vec1(shape);
    vec1.SetData(vec0.data_block(), shape);
    EXPECT_NEAR_VECTOR(vec0, vec1, N, 1e-10);
    EXPECT_EQ((void*)vec0.data_block(), (void*)vec1.GetData());
    }
    {
    // copy in and out
    ArrayType vec1(shape);
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
    ArrayType vec1(vec0.data_block(), shape);
    utl::NDArray<int, Dim> vec2(vec1), vec3;
    vec3 = vec1;
    ArrayType vec4(vec2), vec5;
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

TEST(utlNDArray, Constructors)
{
    {
    unsigned shape[1];
    shape[0]=12;
    test_utlNDArray_Constructors<1>(shape);
    }

    {
    unsigned shape[2];
    shape[0]=3, shape[1]=4;
    test_utlNDArray_Constructors<2>(shape);
    }
    
    {
    unsigned shape[3];
    shape[0]=3, shape[1]=4, shape[2]=5;
    test_utlNDArray_Constructors<3>(shape);
    }
}

TEST(utlNDArray, MaxMin)
{
  unsigned shape[2];
  shape[0]=3, shape[1]=4;
  int N = shape[0]*shape[1];
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  utl::NDArray<double,2> vec1(vec0.data_block(), shape);
    {
    EXPECT_EQ(vec0.arg_max(), vec1.ArgMax());
    EXPECT_EQ(vec0.arg_min(), vec1.ArgMin());
    EXPECT_EQ(vec0.max_value(), vec1.MaxValue());
    EXPECT_EQ(vec0.min_value(), vec1.MinValue());
    }
}

TEST(utlNDArray, Norms)
{
  unsigned shape[2];
  shape[0]=3, shape[1]=4;
  int N = shape[0]*shape[1];
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  utl::NDArray<double,2> vec1(vec0.data_block(), shape);
    {
    EXPECT_NEAR(vec0.two_norm(), vec1.GetArrayTwoNorm(), vec0.two_norm()*1e-10);
    EXPECT_NEAR(vec0.squared_magnitude(), vec1.GetSquaredTwoNorm(), vec0.squared_magnitude()*1e-10);
    EXPECT_NEAR(vec0.inf_norm(), vec1.GetArrayInfNorm(), vec0.inf_norm()*1e-10);
    EXPECT_NEAR(vec0.one_norm(), vec1.GetArrayOneNorm(), vec0.one_norm()*1e-10);
    EXPECT_NEAR(vec0.rms(), vec1.GetRootMeanSquares(), vec0.one_norm()*1e-10);
    }
}

#define __utlNDArray_Operators(vec0, shape, funcName, funcReal) \
  do \
    {\
    utl::NDArray<double,2> vec1(vec0.data_block(), shape); \
    vec1.funcName(); \
    for ( int i = 0; i < N; ++i ) \
      EXPECT_NEAR(funcReal(vec0[i]), vec1[i], 1e-10);\
    } while ( 0 )


TEST(utlNDArray, Operators)
{
  unsigned shape[2];
  shape[0]=3, shape[1]=4;
  int N = shape[0]*shape[1];

  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0), vec2;

  __utlNDArray_Operators(vec0, shape, ElementAbsolute, std::fabs);
  __utlNDArray_Operators(vec0, shape, ElementInverse, 1.0/);
    {
    vnl_vector<double> vv = vec1 +5;
    __utlNDArray_Operators(vv, shape, ElementSqrt, std::sqrt);
    }

  utl::NDArray<double,2> d0(vec0.data_block(), shape); 
  utl::NDArray<double,2> d1(vec1.data_block(), shape), d2; 
    {
    // ()
    unsigned index[2];
    index[0]=1, index[1]=2;
    EXPECT_EQ(d0.GetOffset(index), index[0]*shape[1]+index[1]);
    EXPECT_NEAR(d0(index), d0[index[0]*shape[1]+index[1]], 1e-10);
    }
    {
    bool true_ = d0==d0;
    bool false_ = (d0==d0+1);
    EXPECT_EQ(true_,true);
    EXPECT_EQ(false_,false);
    d2 = d0;
    true_ = d2==d0;
    EXPECT_EQ(true_,true);
    }
    {
    // +
    vec2 = vec0+vec1;
    d2 = d0+d1;
    EXPECT_NEAR_VECTOR(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // -
    vec2=vec0; vec2=-vec1;
    d2=d0; d2=-d1;
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
    utl::NDArray<double,2> dvv(vv.data_block(), shape);
    d2=d0; d2/=dvv; 
    EXPECT_NEAR_VECTOR(vec2.data_block(), d2.GetData(), N, 1e-10);
    }

    {
    // dot
    EXPECT_NEAR(dot_product(vec0, vec1), utl::InnerProduct(d0,d1), std::fabs(dot_product(vec0, vec1))*1e-10);
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
    // friend op, expression
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
    utl::NDArray<double,2> d3(d1);
    d3.ElementAxpby(d0.GetData(), 3.0, -1.0);
    EXPECT_NEAR_VECTOR(d2.GetData(), d3.GetData(), N, 1e-10);
    }
}

TEST(utlNDArray, Operators_TimeCost)
{
  unsigned shape[2];
  shape[0]=40, shape[1]=20;
  int N = shape[0]*shape[1];
  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0), vec2(N);
  vnl_vector<double> vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  
  utl::NDArray<double,2> d0(vec0.data_block(), shape); 
  utl::NDArray<double,2> d1(vec1.data_block(), shape), d2(shape); 


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

TEST(utlNDArray, Functor)
{
  unsigned shape[2];
  shape[0]=3, shape[1]=4;
  int N = shape[0]*shape[1];
  vnl_vector<double> v1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> v2 = __GenerateRandomVector<double>(N, -2.0, 2.0), v3(N);
  utl::NDArray<double,2> vec1(v1.data_block(), shape), vec2(v2.data_block(), shape), vec3, vec4;

    {
    vec3 = utl::F<utl::Maximum<double> >(vec1,vec2);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = v1[i]>v2[i] ? v1[i] : v2[i];
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    }
    {
    vec3 = utl::F<utl::Absolute<double> >(vec1);
    vec1.Apply(utl::Absolute<double>(), vec4);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = std::abs(v1[i]);
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    EXPECT_NEAR_VECTOR(vec4, v3, N, 1e-10);
    }
    {
    vec3 = utl::F<utl::Absolute<double> >(vec1) + utl::F<utl::Maximum<double> >(vec1,vec2);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = (v1[i]>v2[i] ? v1[i] : v2[i])  + std::abs(v1[i]);
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    }
}

