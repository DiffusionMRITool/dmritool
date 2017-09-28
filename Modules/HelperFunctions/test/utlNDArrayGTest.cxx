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
void
__GenerateRandomUtlVector ( const int N, const double val1, const double val2,  utl::Vector<T>& vec)
{
  vec.ReSize(N);
  for ( int i = 0; i < N; i += 1 ) 
    vec[i] = utl::Random<double>(val1,val2);
}

template <class T>
vnl_vector<T> 
__GenerateRandomVector ( const int N, const double val1, const double val2 )
{
  vnl_vector<T> v1(N);
  for ( int i = 0; i < N; i += 1 ) 
    v1[i] = utl::Random<double>(val1,val2);
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

template <class T, unsigned Dim>
void 
test_utlNDArray_Constructors ( unsigned shape[Dim] )
{
  int N = 1;
  for ( int i = 0; i < Dim; ++i ) 
    N *= shape[i];
  vnl_vector<T> vec0 = __GenerateRandomVector<T>(N, -2.0, 2.0);
  // std::cout << "vec0 = " << vec0 << std::endl << std::flush;
  // utl::PrintVnlVector(vec0, "vec0");
  typedef utl::NDArray<T,Dim> ArrayType;
    {
    // data constructor
    ArrayType vec1(vec0.data_block(), shape), vec2(shape);
    EXPECT_NEAR_VECTOR_COMPLEX(vec0, vec1, N, 1e-10);
    vec2 = utl::VnlVectorToStdVector(vec0);;
    EXPECT_NEAR_VECTOR_COMPLEX(vec0, vec2, N, 1e-10);

     // subarray
      {
      int subInd = 1;
      utl::NDArray<T, ArrayType::SubDimension> vec1Sub = vec1.GetRefSubArray(subInd);
      // utl::PrintUtlNDArray(vec1, "vec1:\n");
      // utl::PrintUtlNDArray(vec1Sub, "vec1Sub:\n");

      const unsigned* offsettable = vec1.GetOffSetTable();
      int sizeSub = Dim==1? 1 : offsettable[1];
      EXPECT_EQ(sizeSub, vec1Sub.GetSize());
      EXPECT_FALSE(utl::IsSameShape(vec1,vec1Sub));

      unsigned ss[Dim];
      for ( int i = 0; i < Dim; ++i ) 
        ss[i]=0;
      ss[0]=subInd;
      int size=vec1.GetData()+vec1.Size()-vec1Sub.GetData();
      int off = vec1.GetOffset(ss);
      EXPECT_EQ((void*)(vec1.GetData()+ off), (void*)(vec1Sub.GetData()) );
      T* data1 = vec1.GetData()+ off;
      T* data2 = vec1Sub.GetData();
      EXPECT_NEAR_VECTOR_COMPLEX(data1, data2, size, 1e-8 );
      vec1Sub[1]*=100;
      EXPECT_NEAR_VECTOR_COMPLEX(data1, data2, size, 1e-8 );
      }

      {
      int subInd1 = 1, subInd2=2;
      utl::NDArray<T, Dim> vec1Sub = vec1.GetRefSubArray(subInd1, subInd2);
      // utl::PrintUtlNDArray(vec1, "vec1:\n");
      // utl::PrintUtlNDArray(vec1Sub, "vec1Sub:\n");
      
      const unsigned* offsettable = vec1.GetOffSetTable();
      int sizeSub = (Dim==1?1:offsettable[1]) *(subInd2-subInd1+1);
      EXPECT_EQ(sizeSub, vec1Sub.GetSize());

      unsigned ss[Dim];
      for ( int i = 0; i < Dim; ++i ) 
        ss[i]=0;
      ss[0]=subInd1;
      int size=vec1.GetData()+vec1.Size()-vec1Sub.GetData();
      int off = vec1.GetOffset(ss);
      EXPECT_EQ((void*)(vec1.GetData()+ off), (void*)(vec1Sub.GetData()) );
      T* data1 = vec1.GetData()+ off;
      T* data2 = vec1Sub.GetData();
      EXPECT_NEAR_VECTOR_COMPLEX(data1, data2, size, 1e-8 );
      vec1Sub[1]*=100;
      EXPECT_NEAR_VECTOR_COMPLEX(data1, data2, size, 1e-8 );
      }
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
    EXPECT_NEAR_VECTOR_COMPLEX(vec0, vec2, N, 1e-10);
    ArrayType vec3(vec1);
    EXPECT_NEAR_VECTOR_COMPLEX(vec0, vec3, N, 1e-10);
    }
    {
    // data reference 
    ArrayType vec1(shape);
    vec1.SetData(vec0.data_block(), shape);
    EXPECT_NEAR_VECTOR_COMPLEX(vec0, vec1, N, 1e-10);
    EXPECT_EQ((void*)vec0.data_block(), (void*)vec1.GetData());
    }
    {
    // copy in and out
    ArrayType vec1(shape);
    vec1=10;
    for ( int i = 0; i < N; ++i ) 
      EXPECT_NEAR(std::abs(10.0-vec1[i]), 0.0, 1e-10);
    T* data = vec0.data_block();
    vec1.CopyIn(data,N);
    EXPECT_NEAR_VECTOR_COMPLEX(data, vec1, N, 1e-10);
    vec1[2] = 100; 
    vec1.CopyOut(data,N);
    EXPECT_NEAR_VECTOR_COMPLEX(data, vec1, N, 1e-10);
    }
}

template <class T, unsigned Dim>
utl::NDArray<T,Dim>
__test_move_copy(const utl::NDArray<T,Dim>& arr, void*& address, utl::NDArray<T,Dim>* arr2)
{
  utl::NDArray<T,Dim> result;
  result = arr % 2.0;
  address = result.GetData();
  arr2 = &result;
  // utlPrintVar(true, &result, result.GetData(), address );
  return result;
}


template <class T, unsigned Dim>
void 
test_utlNDArray_Constructors_Convert ( unsigned shape[Dim] )
{
  int N = 1;
  for ( int i = 0; i < Dim; ++i ) 
    N *= shape[i];
  vnl_vector<T> vec0 = __GenerateRandomVector<T>(N, -2.0, 2.0);
  typedef utl::NDArray<T,Dim> ArrayType;
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

      {
      // move constructor, move assignment
      ArrayType vecTmp = vec1;
      T* data = vec1.GetData();
      ArrayType vec11(std::move(vec1)); 
      ArrayType vec12;
      vec12 = std::move(vec11);
      EXPECT_TRUE(vec1.IsEmpty());
      EXPECT_TRUE(vec11.IsEmpty());
      EXPECT_NEAR_UTLVECTOR(vec12, vecTmp, 1e-10);
      EXPECT_EQ(data, vec12.GetData());

      void *address;
      ArrayType* arr2=nullptr;
      vecTmp = __test_move_copy(vec12, address, arr2);
      // utlPrintVar(true, &vecTmp, vecTmp.GetData(),  address, arr2);
      EXPECT_EQ(address, vecTmp.GetData());
      EXPECT_EQ(arr2, nullptr);
      }
    utl::NDArray<std::complex<T>, Dim> vecTmp;
    vecTmp = vec5;
    EXPECT_NEAR_VECTOR_COMPLEX(vec5, vecTmp, vec5.Size(), 1e-10);
    }
}

TEST(utlNDArray, Constructors)
{
    {
    unsigned shape[1];
    shape[0]=12;
    test_utlNDArray_Constructors<double, 1>(shape);
    test_utlNDArray_Constructors<float, 1>(shape);
    test_utlNDArray_Constructors_Convert<double,1>(shape);
    test_utlNDArray_Constructors_Convert<float,1>(shape);
    test_utlNDArray_Constructors<std::complex<double>, 1>(shape);
    }

    {
    unsigned shape[2];
    shape[0]=3, shape[1]=4;
    test_utlNDArray_Constructors<double, 2>(shape);
    test_utlNDArray_Constructors<float, 2>(shape);
    test_utlNDArray_Constructors_Convert<double,2>(shape);
    test_utlNDArray_Constructors_Convert<float,2>(shape);
    test_utlNDArray_Constructors<std::complex<double>, 2>(shape);
    }
    
    {
    unsigned shape[3];
    shape[0]=3, shape[1]=4, shape[2]=5;
    test_utlNDArray_Constructors<double, 3>(shape);
    test_utlNDArray_Constructors<float, 3>(shape);
    test_utlNDArray_Constructors_Convert<double,3>(shape);
    test_utlNDArray_Constructors_Convert<float,3>(shape);
    test_utlNDArray_Constructors<std::complex<double>, 3>(shape);
    }
    
    {
    unsigned shape[4];
    shape[0]=3, shape[1]=4, shape[2]=5, shape[3]=4;
    test_utlNDArray_Constructors<double, 4>(shape);
    test_utlNDArray_Constructors<float, 4>(shape);
    test_utlNDArray_Constructors_Convert<double,4>(shape);
    test_utlNDArray_Constructors_Convert<float,4>(shape);
    test_utlNDArray_Constructors<std::complex<double>,4>(shape);
    }
}

TEST(utlNDArray, ConvertComplex)
{
  unsigned shape[2];
  shape[0]=3, shape[1]=4;
  int N = shape[0]*shape[1];

  vnl_vector<std::complex<double> > vec0 = __GenerateRandomVector<std::complex<double> >(N, -2.0, 2.0);
  utl::NDArray<std::complex<double>, 2> vec1(vec0.data_block(), shape);
  utl::NDArray<double, 2> vec2( reinterpret_cast<double*>(vec0.data_block()), shape);
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

    {
    vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
    utl::NDArray<double,2> vec1(vec0.data_block(), shape);

    EXPECT_NEAR(vec0.two_norm(), vec1.GetArrayTwoNorm(), vec0.two_norm()*1e-10);
    EXPECT_NEAR(vec0.squared_magnitude(), vec1.GetSquaredTwoNorm(), vec0.squared_magnitude()*1e-10);
    EXPECT_NEAR(vec0.inf_norm(), vec1.GetArrayInfNorm(), vec0.inf_norm()*1e-10);
    EXPECT_NEAR(vec0.one_norm(), vec1.GetArrayOneNorm(), vec0.one_norm()*1e-10);
    EXPECT_NEAR(vec0.rms(), vec1.GetRootMeanSquares(), vec0.one_norm()*1e-10);
    }
    
    {
    vnl_vector<std::complex<double> > vec0 = __GenerateRandomVector<std::complex<double> >(N, -2.0, 2.0);
    utl::NDArray<std::complex<double>, 2> vec1(vec0.data_block(), shape);

    EXPECT_NEAR(vec0.two_norm(), vec1.GetArrayTwoNorm(), vec0.two_norm()*1e-10);
    EXPECT_NEAR(vec0.squared_magnitude(), vec1.GetSquaredTwoNorm(), vec0.squared_magnitude()*1e-10);
    EXPECT_NEAR(vec0.rms(), vec1.GetRootMeanSquares(), vec0.one_norm()*1e-10);
    }
}

#define __utlNDArray_Operators(vec0, shape, funcName, funcReal, type) \
  do \
    {\
    utl::NDArray<type,2> vec1(vec0.data_block(), shape); \
    vec1.funcName(); \
    for ( int i = 0; i < N; ++i ) \
      EXPECT_NEAR(std::abs(funcReal(vec0[i])-vec1[i]), 0.0, 1e-10);\
    } while ( 0 )


TEST(utlNDArray, Operators)
{
  unsigned shape[2];
  shape[0]=3, shape[1]=4;
  int N = shape[0]*shape[1];

  vnl_vector<double> vec0 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> vec1 = __GenerateRandomVector<double>(N, -2.0, 2.0), vec2;

  __utlNDArray_Operators(vec0, shape, ElementAbsolute, std::fabs, double);
  __utlNDArray_Operators(vec0, shape, ElementInverse, 1.0/, double);
    {
    vnl_vector<double> vv = vec1 +5;
    __utlNDArray_Operators(vv, shape, ElementSqrt, std::sqrt, double);
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
    EXPECT_NEAR(dot_product(vec0, vec1), utl::DotProduct(d0,d1), std::fabs(dot_product(vec0, vec1))*1e-10);
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
    double aa=3.0;
    d2 = aa%d1;
    d2 = d1%aa;
    d2 = d1%3.0;
    vec2 = vec1*3.0;
    EXPECT_NEAR_VECTOR(d2.GetData(), vec2.data_block(), N, 1e-10);

    double bb=3.0;
    d2 = bb%d1;
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

TEST(utlNDArray, OperatorsComplex)
{
  typedef std::complex<double> VT;
  unsigned shape[2];
  shape[0]=1, shape[1]=2;
  int N = shape[0]*shape[1];

  vnl_vector<VT> vec0 = __GenerateRandomVector<VT>(N, -2.0, 2.0);
  vnl_vector<VT> vec1 = __GenerateRandomVector<VT>(N, -2.0, 2.0), vec2;
  // vec0[0]=std::complex<double>(1,0), vec0[1]=std::complex<double>(0,1);
  // vec1[0]=std::complex<double>(1,0), vec1[1]=std::complex<double>(0,1);
  // vec0[0]=std::complex<double>(-0.343329,-1.37744), vec0[1]=std::complex<double>(-1.31613,-0.0469455);
  // vec1[0]=std::complex<double>(-0.0648453,1.74731), vec1[1]=std::complex<double>(1.35157,0.781885);

  // __utlNDArray_Operators(vec0, shape, ElementAbsolute<double>, std::abs, VT);
  // __utlNDArray_Operators(vec0, shape, ElementInverse, 1.0/);
    {
    vnl_vector<VT> vv = vec1 +5;
    __utlNDArray_Operators(vv, shape, ElementSqrt, std::sqrt, VT);
    }

  utl::NDArray<VT,2> d0(vec0.data_block(), shape); 
  utl::NDArray<VT,2> d1(vec1.data_block(), shape), d2; 
    {
    // ()
    unsigned index[2];
    index[0]=1, index[1]=2;
    EXPECT_EQ(d0.GetOffset(index), index[0]*shape[1]+index[1]);
    EXPECT_NEAR_COMPLEX(d0(index), d0[index[0]*shape[1]+index[1]], 1e-10);
    }
    {
    bool true_ = d0==d0;
    bool false_ = (d0==d0+1);
    bool false2_ = (d0==d0+std::complex<double>(1.0,0.0));
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
    EXPECT_NEAR_VECTOR_COMPLEX(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // -
    vec2=vec0; vec2=-vec1;
    d2=d0; d2=-d1;
    EXPECT_NEAR_VECTOR_COMPLEX(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // -
    vec2=vec0; vec2-=vec1;
    d2=d0; d2-=d1;
    EXPECT_NEAR_VECTOR_COMPLEX(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // *
    vec2=vec0; vec2=element_product(vec2, vec1);
    d2=d0; d2%=d1;
    EXPECT_NEAR_VECTOR_COMPLEX(vec2.data_block(), d2.GetData(), N, 1e-10);
    }
    {
    // /
    vnl_vector<VT> vv = vec1 +5;
    vec2=vec0; vec2=element_quotient(vec2, vv);
    utl::NDArray<VT,2> dvv(vv.data_block(), shape);
    d2=d0; d2/=dvv; 
    EXPECT_NEAR_VECTOR_COMPLEX(vec2.data_block(), d2.GetData(), N, 1e-10);
    }

    {
    // dot
    std::complex<double> sumTmp(0,0), sumTmp1(0.0);
    for ( int i = 0; i < d0.Size(); ++i ) 
      {
      sumTmp += std::conj(d0[i]) * d1[i];
      sumTmp1 += d0[i] * d1[i];
      }
    EXPECT_NEAR_COMPLEX(sumTmp, utl::InnerProduct(d0,d1), std::abs(sumTmp)*1e-10);
    EXPECT_NEAR_COMPLEX(sumTmp1, utl::DotProduct(d0,d1), std::abs(sumTmp1)*1e-10);
    }

    {
    // swap
    d0.Swap(d1);
    EXPECT_NEAR_VECTOR_COMPLEX(vec0.data_block(), d1.GetData(), N, 1e-10);
    EXPECT_NEAR_VECTOR_COMPLEX(vec1.data_block(), d0.GetData(), N, 1e-10);
    d0.Swap(d1);
    }

    {
    // flip
    d2 = d1;
    d2.Flip();
    for ( int i = 0; i < N/2; ++i ) 
      EXPECT_NEAR_COMPLEX(d1[i], d2[N-1-i], 1e-10);
    }

    {
    // friend op, expression
    d2 = d1%3.0;
    vec2 = vec1*3.0;
    EXPECT_NEAR_VECTOR_COMPLEX(d2.GetData(), vec2.data_block(), N, 1e-10);

    d2 = 3.0%d1;
    vec2 = VT(3.0,0.0)*vec1;
    EXPECT_NEAR_VECTOR_COMPLEX(d2.GetData(), vec2.data_block(), N, 1e-10);

    d2 = 3.0+d1;
    vec2 = vnl_vector<VT>(vec1.size(), 3)+vec1;
    EXPECT_NEAR_VECTOR_COMPLEX(d2.GetData(), vec2.data_block(), N, 1e-10);
    
    d2 = 3.0-d1;
    vec2 = vnl_vector<VT>(vec1.size(), 3)-vec1;
    EXPECT_NEAR_VECTOR_COMPLEX(d2.GetData(), vec2.data_block(), N, 1e-10);
    
    d2 = 3.0/(3.0+d1);
    for ( int i = 0; i < N; ++i ) 
      EXPECT_NEAR_COMPLEX(3.0/(3.0+vec1[i]), d2[i], 1e-10);
    
    d2 = (3.0-d1)%2.0;
    for ( int i = 0; i < N; ++i ) 
      EXPECT_NEAR_COMPLEX((3.0-vec1[i])*2.0, d2[i], 1e-10);
    }

    {
    d2 = 3.0%d0-d1;
    utl::NDArray<VT,2> d3(d1);
    d3.ElementAxpby(d0.GetData(), 3.0, -1.0);
    EXPECT_NEAR_VECTOR_COMPLEX(d2.GetData(), d3.GetData(), N, 1e-10);
    }

    {
    vnl_vector<VT> v1 = __GenerateRandomVector<VT>(N, -2.0, 2.0);
    utl::NDArray<VT,1> vec1(v1.data_block(), shape);
    utl::NDArray<double,1> vec2;
    vec2= utl::Real(vec1);
    EXPECT_EQ(vec1.Size(), vec2.Size());
    for ( int i = 0; i < vec1.Size(); ++i ) 
      EXPECT_NEAR(std::real(vec1[i]), vec2[i], 1e-10);
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
    vec3 = utl::F<utl::Functor::Max<double> >(vec1,vec2);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = v1[i]>v2[i] ? v1[i] : v2[i];
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    }
    {
    vec3 = utl::F<utl::Functor::Abs<double> >(vec1);
    vec1.Apply(utl::Functor::Abs<double>(), vec4);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = std::abs(v1[i]);
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    EXPECT_NEAR_VECTOR(vec4, v3, N, 1e-10);
    vec4 = utl::Abs(vec1);
    EXPECT_NEAR_VECTOR(vec4, v3, N, 1e-10);
    }
    {
    vec3 = utl::F<utl::Functor::Log<double> >(utl::F<utl::Functor::Abs<double> >(vec1)) + utl::F<utl::Functor::Max<double> >(vec1,vec2);
    vec4 = utl::Log(utl::Abs(vec1)) + utl::Max(vec1,vec2);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = std::log(std::abs(v1[i])) + (v1[i]>v2[i] ? v1[i] : v2[i]);
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    EXPECT_NEAR_VECTOR(vec4, v3, N, 1e-10);
    }
    {
    vec3 = utl::Pow(vec1, 2.0)+3.0 + utl::Round(vec1) + utl::Sign(vec2) + utl::Neg(vec1);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = std::pow(vec1[i],2.0)+3.0 + std::round(vec1[i]) + utl::sign(vec2[i]) - vec1[i];
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    vec3 = 3.0+utl::Pow(2.0, vec1);
    for ( int i = 0; i < N; ++i ) 
      v3[i] = std::pow(2.0,vec1[i])+3.0;
    EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
    }
}

TEST(utlNDArray, Functor_TimeCost)
{
  unsigned shape[2];
  shape[0]=1e5, shape[1]=10;
  int N=shape[0]*shape[1];
  vnl_vector<double> v1 = __GenerateRandomVector<double>(N, -2.0, 2.0);
  vnl_vector<double> v2 = __GenerateRandomVector<double>(N, -2.0, 2.0), v3(N);
  utl::NDArray<double,2> vec1(v1.data_block(), shape), vec2(v2.data_block(), shape), vec3;

  utl::Tic(std::cout<<"elementwise raw: \n");
  for ( int i = 0; i < N; ++i ) 
    v3[i] = std::log(std::abs(v1[i])) + (v1[i]>v2[i] ? v1[i] : v2[i]) + (v1[i]>0?v1[i]:0) + std::pow(2.0,v2[i]);
  utl::Toc();
  utl::Tic(std::cout<<"utl functor : \n");
  vec3 = utl::Log(utl::Abs(vec1)) + utl::Max(vec1,vec2) + utl::Max(vec1, 0.0) + utl::Pow(2.0,vec2);
  utl::Toc();
  EXPECT_NEAR_VECTOR(vec3, v3, N, 1e-10);
  utl::NDArray<double,2> vec4 = utl::Functor::Square<utl::NDArray<double,2> >()(vec3);
  utl::NDArray<double,2> vec5(vec3);
  for ( int i = 0; i < N; ++i ) 
    vec5[i] = vec3[i]*vec3[i];
  EXPECT_NEAR_VECTOR(vec4, vec5, N, 1e-10);
}

TEST(utlNDArray, Functor_ScalarFunctorWrapper)
{
  utl::Functor::ScalarFunctorWrapper<utl::Functor::Square<double> > func;
  int N=10;
  utl::Vector<double> vec, vec1, vec2;
  __GenerateRandomUtlVector(N, -1.0, 1.0, vec);

  vec1 = func(vec);
  vec2 = utl::Square(vec);
  EXPECT_NEAR_VECTOR(vec1, vec2, N, 1e-10);
  
  utl::Functor::ScalarFunctorWrapper<utl::Functor::Square<double> > func2 = func;
}

TEST(utlNDArray, Functor_VectorUnaryFunctionWrapper)
{
  typedef utl::Functor::VectorUnaryFunctionWrapper<> functorWrapper;
  auto unaryFunc = 
    [](double xval)
      {
      return xval*xval +1;
      }; 
  functorWrapper func(unaryFunc);
  
  int N=10;
  utl::Vector<double> vec, vec1, vec2;
  __GenerateRandomUtlVector(N, -1.0, 1.0, vec);
  vec1 = vec%vec +1.0;
  vec2 = func(vec);
  EXPECT_NEAR_VECTOR(vec1, vec2, N, 1e-10);
}

TEST(utlNDArray, Functor_VectorMultiVariableFunctionWrapper)
{
  typedef utl::Functor::VectorMultiVariableFunctionWrapper<> functorWrapper;
  auto autoFunc = 
    [](const std::vector<double>& vec)
      {
      double sum=0;
      for ( int i = 0; i < vec.size(); ++i ) 
        sum += vec[i];
      return sum;
      }; 
  functorWrapper func(autoFunc);
  
  int N=10;
  utl::Vector<double> vec, vec1, vec2, vec3;
  __GenerateRandomUtlVector(N, -1.0, 1.0, vec);
  __GenerateRandomUtlVector(N, -1.0, 1.0, vec1);
  std::vector<utl::Vector<double> > vecs;
  vecs.push_back(vec);
  vecs.push_back(vec1);
  vec2 = func(vecs);
  vec3 = vec+vec1;
  EXPECT_NEAR_VECTOR(vec3, vec2, N, 1e-10);
}

