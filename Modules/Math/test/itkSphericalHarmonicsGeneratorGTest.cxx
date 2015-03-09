/**
 *       @file  itkSphericalHarmonicsGeneratorGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-08-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "itkSphericalHarmonicsGenerator.h"

typedef itk::SphericalHarmonicsGenerator<double> SHGenerator;

TEST(itkSphericalHarmonicsGenerator, ComplexSH)
{
  EXPECT_NEAR(-0.2891664529927667, std::real(SHGenerator::ComplexSH(2,1, 0.530274, 0.539671)), 1e-8); 
  EXPECT_NEAR(-0.17320564388961293, std::imag(SHGenerator::ComplexSH(2,1, 0.530274, 0.539671)), 1e-8); 
  
  EXPECT_NEAR(27.0/256.0*std::sqrt(35.0/(2.0*M_PI)), std::real(SHGenerator::ComplexSH(4,-4, M_PI/3.0, M_PI/2.0)), 1e-8); 
  EXPECT_NEAR(0.0, std::imag(SHGenerator::ComplexSH(4,-4, M_PI/3.0, M_PI/2.0)), 1e-8); 
  
  EXPECT_NEAR(0.125613008425, std::real(SHGenerator::ComplexSH(4,3, -M_PI/3.0, M_PI/5.0)), 1e-8); 
  EXPECT_NEAR(-0.386597088084, std::imag(SHGenerator::ComplexSH(4,3, -M_PI/3.0, M_PI/5.0)), 1e-8); 
  EXPECT_NEAR(0.125613008425, std::real(SHGenerator::ComplexSH(4,3, M_PI/3.0, M_PI/5.0)), 1e-8); 
  EXPECT_NEAR(-0.386597088084, std::imag(SHGenerator::ComplexSH(4,3, M_PI/3.0, M_PI/5.0)), 1e-8); 
}

TEST(itkSphericalHarmonicsGenerator, RealSH)
{
  EXPECT_NEAR(-0.2449497706682552, SHGenerator::RealSH(2,1, 0.530274, 0.539671), 1e-8); 
  EXPECT_NEAR(0.352032601190162, SHGenerator::RealSH(4,-4, M_PI/3.0, M_PI/2.0), 1e-8); 
  EXPECT_NEAR(0.282094791773878, SHGenerator::RealSH(0,0, M_PI/3.0, M_PI/2.0), 1e-8); 
  EXPECT_NEAR(0.282094791773878, SHGenerator::RealSH(0,0, 0, 0), 1e-8); 
  EXPECT_NEAR(-0.546730845142995, SHGenerator::RealSH(4,3, -M_PI/3.0, M_PI/5.0), 1e-8); 
  EXPECT_NEAR(-0.546730845142995, SHGenerator::RealSH(4,3, M_PI/3.0, M_PI/5.0), 1e-8); 
}

#define __itkSphericalHarmonicsGenerator_ComplexTripleIntegration(l1,m1,l2,m2,theta,phi)                                 \
do { \
  std::complex<double> y1, y2, y3;                                                                                             \
  y1 = SHGenerator::ComplexSH(l1,m1, theta, phi); \
  y2 = SHGenerator::ComplexSH(l2,m2, theta, phi); \
  y3 = 0; \
  int maxL = 2*std::max(l1,l2);        \
  for ( int l3 = 0; l3 <= maxL; l3 += 1 )  \
    for ( int m3 = -l3; m3 <= l3; m3 += 1 )  \
      y3 +=  SHGenerator::ComplexTripleIntegration(l1,m1,l2,m2,l3,m3) * std::conj(SHGenerator::ComplexSH(l3,m3, theta, phi)); \
  EXPECT_NEAR(std::real(y3), std::real(y1*y2), 1e-8);  \
  EXPECT_NEAR(std::imag(y3), std::imag(y1*y2), 1e-8); \
} while(0)


TEST(itkSphericalHarmonicsGenerator, ComplexTripleIntegration)
{
  __itkSphericalHarmonicsGenerator_ComplexTripleIntegration(4,3,2,1,M_PI/3.0,M_PI/5.0);
  __itkSphericalHarmonicsGenerator_ComplexTripleIntegration(2,1,5,-3,M_PI/3.0,M_PI/5.0);
  __itkSphericalHarmonicsGenerator_ComplexTripleIntegration(0,0,0,0,M_PI/3.0,M_PI/5.0);

  int N = 10;
  for ( int i = 0; i < N; ++i ) 
    {
    int l1 = utl::RandomInt(0,20), l2 = utl::RandomInt(0,20);
    int m1 = utl::RandomInt(-l1,l1), m2 = utl::RandomInt(-l2,l2);
    double theta = utl::Random<double>(0,M_PI), phi=utl::Random<double>(0, 2*M_PI);
    __itkSphericalHarmonicsGenerator_ComplexTripleIntegration(l1,m1,l2,m2,theta,phi);
    }
}

#define __itkSphericalHarmonicsGenerator_RealTripleIntegration(l1,m1,l2,m2,theta,phi)                                 \
do { \
  double y1, y2, y3;                                                                                             \
  y1 = SHGenerator::RealSH(l1,m1, theta, phi); \
  y2 = SHGenerator::RealSH(l2,m2, theta, phi); \
  y3 = 0; \
  int maxL = 2*std::max(l1,l2);        \
  for ( int l3 = 0; l3 <= maxL; l3 += 2 )  \
    for ( int m3 = -l3; m3 <= l3; m3 += 1 )  \
      y3 +=  SHGenerator::RealTripleIntegration(l1,m1,l2,m2,l3,m3) * SHGenerator::RealSH(l3,m3, theta, phi); \
  EXPECT_NEAR(y3, y1*y2, 1e-8);  \
} while(0)

TEST(itkSphericalHarmonicsGenerator, RealTripleIntegration)
{
  __itkSphericalHarmonicsGenerator_RealTripleIntegration(4,3,2,1,M_PI/3.0,M_PI/5.0);
  __itkSphericalHarmonicsGenerator_RealTripleIntegration(2,1,6,-3,M_PI/3.0,M_PI/5.0);
  __itkSphericalHarmonicsGenerator_RealTripleIntegration(6,-5,8,7,-M_PI/3.0,M_PI/5.0);
  __itkSphericalHarmonicsGenerator_RealTripleIntegration(0,0,0,0,M_PI/3.0,M_PI/5.0);

  int N = 10;
  for ( int i = 0; i < N; ++i ) 
    {
    int l1 = utl::RandomInt(0,10)*2, l2 = utl::RandomInt(0,10)*2;
    int m1 = utl::RandomInt(-l1,l1), m2 = utl::RandomInt(-l2,l2);
    double theta = utl::Random<double>(0,M_PI), phi=utl::Random<double>(0, 2*M_PI);
    __itkSphericalHarmonicsGenerator_RealTripleIntegration(l1,m1,l2,m2,theta,phi);
    }
}

