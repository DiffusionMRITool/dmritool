/**
 *       @file  utlCoreGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "07-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#include "utlGTest.h"
#include "utlCore.h"

TEST(utlCore, utlArraySize)
{
  int a1[3];  
  EXPECT_EQ(3, utlArraySize(a1));
  double a2[5];  
  EXPECT_EQ(5, utlArraySize(a2));

    {
    int nums []={2, 1, 2, 3, 4};
    EXPECT_EQ(5, utlArraySize(nums));
    }
}

TEST(utlCore, IsSame)
{
  EXPECT_TRUE(utl::IsSame<std::string>("abc", "abc"));
  EXPECT_TRUE(utl::IsSameArray("abc", "abc"));
  EXPECT_FALSE(utl::IsSame<std::string>("abc", " abc"));
  EXPECT_FALSE(utl::IsSameArray("abc", " abc"));
    {
    int a1[3], a2[3];
    for ( int i = 0; i < 3; ++i ) 
      a1[i]=a2[i]=i;
    EXPECT_TRUE(utl::IsSameArray(a1,a2));
    }
    {
    int a1[3], a2[4];
    for ( int i = 0; i < 3; ++i ) 
      a1[i]=a2[i]=i;
    EXPECT_FALSE(utl::IsSameArray(a1,a2));
    }
    {
    int a1[3], a2[3];
    std::vector<int> vec1, vec2;
    for ( int i = 0; i < 3; ++i ) 
      {
      a1[i]=i, a2[i]=2*i;
      vec1.push_back(a1[i]), vec2.push_back(a2[i]);
      }
    EXPECT_FALSE(utl::IsSameArray(a1,a2));
    EXPECT_FALSE(utl::IsSameVector(vec1,vec2));
    }
}

TEST(UtlCore, ZeroPad)
{
  EXPECT_STREQ("0005", utl::ZeroPad(5,4).c_str()); 
  EXPECT_STREQ("5", utl::ZeroPad(5,1).c_str()); 
  EXPECT_STREQ("dwi_000003.hdr", utl::GetSequentialFileName("dwi", 3, "hdr").c_str()); 
  EXPECT_STREQ("dwi_000003.nii.gz", utl::GetSequentialFileName("dwi", 3, "nii.gz").c_str()); 
}

TEST(UtlCore, ConvertStringToNumber)
{
  EXPECT_DOUBLE_EQ(3.1415,utl::ConvertStringToNumber<double>("3.1415"));
  EXPECT_DOUBLE_EQ(3.1415,utl::ConvertStringToNumber<double>(" 3.1415 "));
  EXPECT_FLOAT_EQ((float)3.1415, utl::ConvertStringToNumber<float>("3.1415"));
  EXPECT_EQ((int)3.1415, utl::ConvertStringToNumber<int>("3.1415"));
  EXPECT_DOUBLE_EQ(-3.1415, utl::ConvertStringToNumber<double>("-3.1415"));
  EXPECT_FLOAT_EQ((float)-3.1415, utl::ConvertStringToNumber<float>("-3.1415"));
  EXPECT_EQ((int)-3.1415, utl::ConvertStringToNumber<int>("-3.1415"));

  EXPECT_DOUBLE_EQ(1.8e3, utl::ConvertStringToNumber<double>("1.8e3"));
  EXPECT_DOUBLE_EQ(1.8e-3, utl::ConvertStringToNumber<double>("1.8e-3"));
}

TEST(UtlCore, RoundNumber)
{
  EXPECT_EQ(0, utl::RoundNumber(0));
  EXPECT_EQ(0, utl::RoundNumber(1e-4));
  EXPECT_EQ(1, utl::RoundNumber(1.3));
  EXPECT_EQ(2, utl::RoundNumber(1.5));
  EXPECT_EQ(2, utl::RoundNumber(1.8));
  EXPECT_EQ(-1, utl::RoundNumber(-1.3));
  EXPECT_EQ(-2, utl::RoundNumber(-1.5) );
  EXPECT_EQ(-2, utl::RoundNumber(-1.8) );
}

TEST(UtlCore, ConvertToUnixSlashes)
{
  std::string str("C:\\Windows\\System32");
  utl::ConvertToUnixSlashes(str);
  EXPECT_STREQ("C:/Windows/System32", str.c_str());
}

TEST(UtlCore, CreateExpandedPath)
{
  EXPECT_STREQ( (std::string(getenv("HOME"))+ "/dir").c_str(), utl::CreateExpandedPath("~/dir/").c_str());
}

TEST(utlCore, GetFileExtension)
{
  std::string ext, file;
  utl::GetPath("/home/dwi.hdr", ext, file);
  EXPECT_STREQ("/home/", ext.c_str());
  EXPECT_STREQ("dwi.hdr", file.c_str());
  utl::GetFileExtension("/home/dwi.hdr", ext, file);
  EXPECT_STREQ("hdr", ext.c_str());
  EXPECT_STREQ("/home/dwi", file.c_str());
  utl::GetFileExtension("/home/dwi.nii.gz", ext, file);
  EXPECT_STREQ("nii.gz", ext.c_str());
  EXPECT_STREQ("/home/dwi", file.c_str());
  utl::GetFileExtension("/home/dwi.gz", ext, file);
  EXPECT_STREQ("gz", ext.c_str());
  EXPECT_STREQ("/home/dwi", file.c_str());
}

TEST(utlCore, IsNumber)
{
  EXPECT_TRUE(utl::IsNumber("0.345"));
  EXPECT_TRUE(utl::IsNumber("1.7e-5"));

  EXPECT_FALSE(utl::IsNumber("1.7a-5"));
  EXPECT_FALSE(utl::IsNumber("3+1.7e-5"));
  EXPECT_FALSE(utl::IsNumber("0..3445"));
  EXPECT_FALSE(utl::IsNumber("31adadt465"));
}


TEST(utlCore, ComputeNDArrayOffset)
{
  std::vector<int> size(2);
  size[0]=3, size[1]=4;
  double off=0;

    {
    std::vector<int> index(2), index2(2);
    index[0]=1, index[1]=3;
    off=10;
    EXPECT_EQ(10, utl::ComputeNDArrayOffset(index, size, COLUMN_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, COLUMN_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 2, 1e-10);

    index[0]=0, index[1]=0;
    off=0;
    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, COLUMN_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, COLUMN_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 2, 1e-10);

    index[0]=2, index[1]=3;
    off=11;
    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, COLUMN_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, COLUMN_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 2, 1e-10);
    }
    
    {
    std::vector<int> index(2), index2(2);
    index[0]=1, index[1]=3;
    off=7;
    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, ROW_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, ROW_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 2, 1e-10);

    index[0]=0, index[1]=0;
    off=0;
    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, ROW_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, ROW_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 2, 1e-10);

    index[0]=2, index[1]=3;
    off=11;
    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, ROW_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, ROW_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 2, 1e-10);
    }

    {
    std::vector<int> size(1);
    size[0] = 10;

    std::vector<int> index(1), index2(1);
    index[0]=5;
    off=5;
    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, ROW_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, ROW_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 1, 1e-10);

    EXPECT_EQ(off, utl::ComputeNDArrayOffset(index, size, COLUMN_MAJOR));
    utl::ComputeNDArrayIndex(off, index2, size, COLUMN_MAJOR);
    EXPECT_NEAR_VECTOR(index, index2, 1, 1e-10);
    }
}

TEST(utlCore, print)
{
  double a_d=2.0;
  int a_i = 1;
  bool a_b = true;
  
    {
    std::ostringstream oss;
    utl::PrintOS(oss, a_d, a_i, a_b);
    EXPECT_EQ(oss.str(), "(2, 1, true)\n");
    }
    {
    std::ostringstream oss;
    PrintVar(true, oss, a_d, a_i, a_b);
    EXPECT_EQ(oss.str(), "(a_d, a_i, a_b) = (2, 1, true)\n");
    }
}

TEST(utlCore, InstanceOf)
{
    {
    std::vector<double> vec;
    EXPECT_EQ(utl::IsInstanceOf<std::vector<double> >(vec), true);
    EXPECT_EQ(utl::IsInstanceOf<std::vector<int> >(vec), false);

    int intNum;
    EXPECT_EQ(utl::IsInstanceOf<float >(intNum), false);
    EXPECT_EQ(utl::IsInstanceOf<int >(intNum), true);
    }
}
