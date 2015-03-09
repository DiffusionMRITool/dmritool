/**
 *       @file  utlMathGTest.cxx
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

#include "utlCore.h"
#include "utlMath.h"
#include "gtest/gtest.h"


TEST(UtlMath, Factorial)
{
  EXPECT_EQ(1, utl::Factorial(0));
  EXPECT_EQ(1, utl::Factorial(1));
  EXPECT_EQ(24, utl::Factorial(4));
  EXPECT_EQ(40320, utl::Factorial(8));
}

TEST(utlMath, GammaHalfInteger)
{
  EXPECT_NEAR( std::sqrt(M_PI), utl::GammaHalfInteger(0.5), 1e-15);
  EXPECT_NEAR( 1, utl::GammaHalfInteger(1), 1e-15);
  EXPECT_NEAR( 362880, utl::GammaHalfInteger(10), 1e-15);
  EXPECT_NEAR( 1133278.3889487855673346, utl::GammaHalfInteger(10.5), 1e-15);
}

TEST(utlMath, PowHalfInteger)
{
  double eps=1e-12;
  double val = utl::Random<double>(-2.0,2.0);
  for ( double p = -10; p < 10; p += 0.5 ) 
    {
    if (val<0 && utl::IsInt(p+0.5))
      EXPECT_NEAR( std::pow(-val,p), utl::PowHalfInteger(-val,p), eps*std::fabs(std::pow(-val,p))) 
        << "val = " << val << ", p = " << p << ", std::pow(-val,p) = " << std::pow(-val,p) <<", utl::PowHalfInteger(-val,p) =" << utl::PowHalfInteger(-val,p);
    else
      EXPECT_NEAR( std::pow(val,p), utl::PowHalfInteger(val,p), eps*std::fabs(std::pow(val,p))) 
        << "val = " << val << ", p = " << p << ", std::pow(val,p) = " << std::pow(val,p) <<", utl::PowHalfInteger(val,p) =" << utl::PowHalfInteger(val,p);
    }
}

TEST(utlMath, PowInteger)
{
  double eps=1e-12;
  double val = utl::Random<double>(-2.0,2.0);
  for ( double p = -10; p < 10; p += 1 ) 
    {
      EXPECT_NEAR( std::pow(val,p), utl::PowInteger(val,p), eps*std::fabs(std::pow(val,p))) 
        << "val = " << val << ", p = " << p << ", std::pow(val,p) = " << std::pow(val,p) <<", utl::PowInteger(val,p) =" << utl::PowHalfInteger(val,p);
    }
}

template <class T>
void 
test ( const std::vector<T>& data, const double val )
{
  int N = data.size();
  double sum1=0, sum2=0;
    {
    // PowHalfInteger
    utl::Tic(std::cout<<"PowHalfInteger for " << val << "\n");
    for ( int i = 0; i < N; i += 1 ) 
      {
      sum1+=utl::PowHalfInteger(val, data[i]);
      }
    utl::Toc();
    }

    {
    // std::pow
    utl::Tic(std::cout<<"std::pow for " << val << "\n");
    for ( int i = 0; i < N; i += 1 ) 
      {
      sum2+=std::pow(val, data[i]);
      }
    utl::Toc();
    }
  EXPECT_NEAR(sum1, sum2, std::fabs(sum1)*1e-10);
}


TEST(utlMath, PowHalfInteger_TimeCost)
{
  long N = 1e7;
  std::vector<double> data;
  std::vector<int> dataInt;
  for ( int i = 0; i < N; i += 1 ) 
    {
    data.push_back(utl::RandomInt(1,20)*0.5);
    dataInt.push_back(utl::RandomInt(1,10));
    }

  double val = utl::Random<double>(0,4);
  test(data, val);
  test(dataInt, val);
}

TEST(utlMath, ExpNegtiveLUT)
{
  double eps=1e-8;
  double distMax = 30.0, precision=1000;
  EXPECT_NEAR(std::exp(-3),utl::ExpNegtiveLUT(3,distMax,precision), eps);
  EXPECT_NEAR(std::exp(0), utl::ExpNegtiveLUT(0,distMax,precision), eps);
  EXPECT_NEAR(std::exp(-1.5), utl::ExpNegtiveLUT(1.5,distMax,precision), eps);
  EXPECT_NEAR(std::exp(-10.5), utl::ExpNegtiveLUT(10.5,distMax,precision), eps);
  EXPECT_NEAR(std::exp(-1e-5), utl::ExpNegtiveLUT(1e-5,distMax,precision), eps);
  EXPECT_NEAR(std::exp(-34), utl::ExpNegtiveLUT(34,distMax,precision), eps);
}

TEST(utlMath, PolynomialRoot)
{
  std::vector<std::complex<double> > root; 
  {
  // order 2
  std::vector<double> coef(3);
  for ( int i = 0; i < 3; ++i ) 
    coef[i] = utl::Random<double>(-3.0,3.0);
  root = utl::PolynomialRoot(coef);
  for ( int i = 0; i < root.size(); ++i ) 
    {
    std::complex<double> re, a=root[i];
    re = coef[0]+ coef[1]*a + coef[2]*a*a;
    EXPECT_NEAR(re.real(), 0.0, 1e-8);
    EXPECT_NEAR(re.imag(), 0.0, 1e-8);
    }
  }
  {
  // order 3
  std::vector<double> coef(4);
  for ( int i = 0; i < 4; ++i ) 
    coef[i] = utl::Random<double>(-3.0,3.0);
  root = utl::PolynomialRoot(coef);
  for ( int i = 0; i < root.size(); ++i ) 
    {
    std::complex<double> re, a=root[i];
    re = coef[0]+ coef[1]*a + coef[2]*a*a + coef[3]*a*a*a;
    EXPECT_NEAR(re.real(), 0.0, 1e-8);
    EXPECT_NEAR(re.imag(), 0.0, 1e-8);
    }
  }
}

