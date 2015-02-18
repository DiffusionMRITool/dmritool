/**
 *       @file  itkUnaryFunctorLookUpTableGTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-05-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "utlCore.h"
#include "itkUnaryFunctorLookUpTable.h"
#include "itkFunctors.h"


TEST(itkUnaryFunctorLookUpTable, EXP)
{
  double eps=1e-8;
  typedef itk::UnaryFunctorLookUpTable<itk::Functor::EXP<double> > LUTType;
  LUTType::Pointer lut = LUTType::New();
  lut->SetVariableMax(0);
  lut->SetVariableMin(-30);
  lut->SetNumberOfBins(30*1e3);
  lut->BuildTable();
  
  double val[]={-3,0,-1.5,-10.5,-1e-5,-34};
  for ( int i = 0; i < 5; i += 1 ) 
    {
    double real_func_val = std::exp(val[i]);
    EXPECT_NEAR(real_func_val, lut->GetFunctionValue(val[i]), eps*std::fabs(real_func_val));
    }
}

TEST(itkUnaryFunctorLookUpTable, EXP_TimeCost)
{
  typedef itk::UnaryFunctorLookUpTable<itk::Functor::EXP<double> > LUTType;
  LUTType::Pointer lut = LUTType::New();
  lut->SetVariableMax(0);
  lut->SetVariableMin(-30);
  lut->SetNumberOfBins(30*1e3);
  lut->BuildTable();

  long N = 1e7;
  std::vector<double> val;
  for ( int i = 0; i < N; i += 1 ) 
    val.push_back(utl::Random<double>(-30,0));

    {
    utl::Tic(std::cout<<"std::exp, " << N << "\n");
    double sum=0;
    for ( int i = 0; i < val.size(); i += 1 ) 
      sum += std::exp(val[i]);
    utl::Toc();
    if (sum<1)
      std::cout << "sum="<<sum << std::endl << std::flush;
    }

    {
    utl::Tic(std::cout<<"LookUpTable, " << N << "\n");
    double sum=0;
    for ( int i = 0; i < val.size(); i += 1 ) 
      sum += lut->GetFunctionValue(val[i]);
    utl::Toc();
    if (sum<1)
      std::cout << "sum="<<sum << std::endl << std::flush;
    }

    {
    utl::Tic(std::cout<<"utl::ExpNegtiveLUT, " << N << "\n");
    double sum=0;
    for ( int i = 0; i < val.size(); i += 1 ) 
      sum += utl::ExpNegtiveLUT(-val[i]);
    utl::Toc();
    if (sum<1)
      std::cout << "sum="<<sum << std::endl << std::flush;
    }
}

TEST(itkUnaryFunctorLookUpTable, COS)
{
  double eps=1e-8;
  typedef itk::UnaryFunctorLookUpTable<itk::Functor::COS<double> > LUTType;
  LUTType::Pointer lut = LUTType::New();
  lut->SetVariableMax(2*M_PI);
  lut->SetVariableMin(0);
  lut->SetNumberOfBins(2*M_PI*1e4);
  lut->BuildTable();

  double val[]={0, 2*M_PI, M_PI, utl::Random<double>(0,2*M_PI), utl::Random<double>(0,2*M_PI)};
  for ( int i = 0; i < 5; i += 1 ) 
    {
    double real_func_val = std::cos(val[i]);
    EXPECT_NEAR(real_func_val, lut->GetFunctionValue(val[i]), eps*std::fabs(real_func_val));
    }
}


TEST(itkUnaryFunctorLookUpTable, COS_TimeCost)
{
  typedef itk::UnaryFunctorLookUpTable<itk::Functor::COS<double> > LUTType;
  LUTType::Pointer lut = LUTType::New();
  lut->SetVariableMax(2*M_PI);
  lut->SetVariableMin(0);
  lut->SetNumberOfBins(2*M_PI*1e4);
  lut->BuildTable();

  long N = 1e7;
  std::vector<double> val;
  for ( int i = 0; i < N; i += 1 ) 
    val.push_back(utl::Random<double>(0,2*M_PI));

    {
    utl::Tic(std::cout<<"std::cos, " << N << "\n");
    double sum=0;
    for ( int i = 0; i < val.size(); i += 1 ) 
      sum += std::cos(val[i]);
    utl::Toc();
    if (sum==0)
      std::cout << "sum="<<sum << std::endl << std::flush;
    }

    {
    utl::Tic(std::cout<<"LookUpTable, " << N << "\n");
    double sum=0;
    for ( int i = 0; i < val.size(); i += 1 ) 
      sum += lut->GetFunctionValue(val[i]);
    utl::Toc();
    if (sum==0)
      std::cout << "sum="<<sum << std::endl << std::flush;
    }
}

