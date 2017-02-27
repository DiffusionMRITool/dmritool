/**
 *       @file  itkUnaryFunctorLookUpTableTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-05-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkUnaryFunctorLookUpTable.h"
#include "utl.h"
#include "itkSpecialFunctionGenerator.h"
#include "itkFunctors.h"

/**
 * \brief  test itkUnaryFunctorLookUpTable
 */
int 
main (int argc, char const* argv[])
{
    {
    // exp 
    std::cout << "\nexp" << std::endl << std::flush;
    double valMax=0.0, valMin=-10.0;
    for ( int i = 0; i < 5; i += 1 ) 
      {
      double val = utl::Random<double>(valMin, valMax);
      if (i==0)
        val = -2.02e-05; 
      else if (i==1)
        val = -29.999;

      std::cout<< std::setprecision(15);
      double val_func = std::exp(val);

      double val_funcLut1 = utl::ExpNegtiveLUT(-val,30,1e3);
      utlPrintVar4(true, val, val_func, val_funcLut1, std::fabs(val_func-val_funcLut1)/val_func);

      typedef itk::UnaryFunctorLookUpTable<utl::Functor::Exp<double> > LUTType;
      LUTType::Pointer lut = LUTType::New();
      lut->SetVariableMax(0);
      lut->SetVariableMin(-30);
      lut->SetNumberOfBins(30*1e3);
      lut->Initialize();
      double val_funcLut2 = lut->GetFunctionValue(val);
      utlPrintVar4(true, val, val_func, val_funcLut2, std::fabs(val_func-val_funcLut2)/val_func);
      }
    }

    {
    // cos
    std::cout << "\ncos" << std::endl << std::flush;
    double valMax=2*M_PI, valMin=0;
    double val = utl::Random<double>(valMin, valMax);

    std::cout<< std::setprecision(15);
    double val_func = std::cos(val);

    typedef itk::UnaryFunctorLookUpTable<utl::Functor::Cos<double> > LUTType;
    LUTType::Pointer lut = LUTType::New();
    lut->SetVariableMax(2*M_PI);
    lut->SetVariableMin(0);
    lut->SetNumberOfBins(2*M_PI*1e3);
    lut->Initialize();
    double val_funcLut2 = lut->GetFunctionValue(val);
    utlPrintVar4(true, val, val_func, val_funcLut2, std::fabs(val_func-val_funcLut2)/std::fabs(val_func));
    }


  return 0;
}
