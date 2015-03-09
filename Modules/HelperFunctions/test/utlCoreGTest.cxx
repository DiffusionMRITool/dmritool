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


#include "gtest/gtest.h"
#include "utlCore.h"


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


