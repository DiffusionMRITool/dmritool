/**
 *       @file  itkFunctorHashTableTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */
#include "utl.h"
#include "itkFunctors.h"
#include "itkFunctorHashTable.h"
#include "itkTimeProbe.h"
#include "itkSpecialFunctionGenerator.h"
#include <gsl/gsl_sf_hyperg.h>

template< class TInput, class TArgument=TInput, class TOutput=TInput >
class GSL_1F1
{
public:
  GSL_1F1() {};
  ~GSL_1F1() {};
  bool operator!=(const GSL_1F1 &) const
  { return false; }
  bool operator==(const GSL_1F1 & other) const
  { return !( *this != other ); }
  inline TOutput operator()( const TInput & A )
  {
    return static_cast<TOutput>( gsl_sf_hyperg_1F1(2,0.5,A ) );
  }

};

template< class TInput, class TArgument=TInput, class TOutput=TInput >
class LAGURRE
{
public:
  LAGURRE() {};
  ~LAGURRE() {};
  bool operator!=(const LAGURRE &) const
  { return false; }
  bool operator==(const LAGURRE & other) const
  { return !( *this != other ); }
  inline TOutput operator()( const TInput & A )
  {
    return static_cast<TOutput>( utl::Lagurre(2,0.5,A ) );
  }

};


/**
 * \brief  test itk::FunctorHashTable
 */
int 
main (int argc, char const* argv[])
{
  itk::TimeProbe clock;
    int N=1e3, iter=1e5;
    // int N=1e7, iter=1;

    {
    typedef itk::FunctorHashTable<itk::Functor::EXP<double>, double, double> FunctorHashType;
    FunctorHashType::Pointer func = FunctorHashType::New();
    func->SetDebug(true);

    // int N=1e2, iter=1e4;
    std::vector<double> valVec;

    for ( int i = 0; i < N; i += 1 ) 
      {
      double val = utl::Random(-3.0,3.0);
      valVec.push_back(val);
      // double func_real = std::exp(val);
      double func_hash = func->GetFunctionValue(val);
      // utlPrintVar4(true, val, func_real, func_hash, std::fabs((func_real-func_hash)/func_real));
      }
    clock.Start();
    for ( int ii = 0; ii < iter; ii += 1 ) 
      {
      double aa=0;
      for ( int i = 0; i < N; i += 1 ) 
        {
        double val = valVec[i];
        double func_hash = func->GetFunctionValue(val);
        }
      }
    clock.Stop();
    std::cout << "Mean (hash exp): " << clock.GetMean() << std::endl;
    clock.Start();
    for ( int ii = 0; ii < iter; ii += 1 ) 
      {
      double aa=0;
      for ( int i = 0; i < N; i += 1 ) 
        {
        double val = valVec[i];
        double func_real = std::exp(val);
        }
      }
    clock.Stop();
    std::cout << "Mean (std exp): " << clock.GetMean() << std::endl;
    }

    {
    typedef itk::FunctorHashTable<GSL_1F1<double>, double, double> FunctorHashType;
    FunctorHashType::Pointer func = FunctorHashType::New();
    func->SetDebug(true);

    std::vector<double> valVec;

    for ( int i = 0; i < N; i += 1 ) 
      {
      double val = utl::Random(0.0,3.0);
      valVec.push_back(val);
      // double func_real = std::exp(val);
      double func_hash = func->GetFunctionValue(val);
      // utlPrintVar4(true, val, func_real, func_hash, std::fabs((func_real-func_hash)/func_real));
      }
    std::cout << "start" << std::endl << std::flush;
    clock.Start();
    for ( int ii = 0; ii < iter; ii += 1 ) 
      {
      double aa=0;
      for ( int i = 0; i < N; i += 1 ) 
        {
        double val = valVec[i];
        double func_hash = func->GetFunctionValue(val);
        }
      }
    clock.Stop();
    std::cout << "Mean (hash 1F1): " << clock.GetMean() << std::endl;
    clock.Start();
    for ( int ii = 0; ii < iter; ii += 1 ) 
      {
      double aa=0;
      for ( int i = 0; i < N; i += 1 ) 
        {
        double val = valVec[i];
        double func_real = gsl_sf_hyperg_1F1(2,0.5,val);
        }
      }
    clock.Stop();
    std::cout << "Mean (gsl 1F1): " << clock.GetMean() << std::endl;
    }

    {
    typedef itk::FunctorHashTable<LAGURRE<double>, double, double> FunctorHashType;
    FunctorHashType::Pointer func = FunctorHashType::New();
    func->SetDebug(true);

    std::vector<double> valVec;

    for ( int i = 0; i < N; i += 1 ) 
      {
      double val = utl::Random(0.0,3.0);
      valVec.push_back(val);
      // double func_real = std::exp(val);
      double func_hash = func->GetFunctionValue(val);
      // utlPrintVar4(true, val, func_real, func_hash, std::fabs((func_real-func_hash)/func_real));
      }
    std::cout << "start" << std::endl << std::flush;
    clock.Start();
    for ( int ii = 0; ii < iter; ii += 1 ) 
      {
      double aa=0;
      for ( int i = 0; i < N; i += 1 ) 
        {
        double val = valVec[i];
        double func_hash = func->GetFunctionValue(val);
        }
      }
    clock.Stop();
    std::cout << "Mean (hash Lagurre): " << clock.GetMean() << std::endl;
    clock.Start();
    for ( int ii = 0; ii < iter; ii += 1 ) 
      {
      double aa=0;
      for ( int i = 0; i < N; i += 1 ) 
        {
        double val = valVec[i];
        double func_real = utl::Lagurre(2,0.5,val);
        }
      }
    clock.Stop();
    std::cout << "Mean (gsl Lagurre): " << clock.GetMean() << std::endl;
    }

  return 0;
}
