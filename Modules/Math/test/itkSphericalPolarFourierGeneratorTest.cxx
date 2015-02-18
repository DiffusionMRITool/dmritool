/**
 *       @file  itkSphericalPolarFourierGeneratorTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-05-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkSphericalPolarFourierGenerator.h"
#include "utl.h"


/**
 * \brief  test itkSphericalPolarFourierGenerator
 */
int 
main (int argc, char const* argv[])
{
  double md0 = 0.7e-3;
  double tau = ONE_OVER_4_PI_2;
  double scale = 1.0 / (8*M_PI*M_PI*tau*md0);
  int N=1, L=4;
  double radius = 0.015;

  typedef itk::SphericalPolarFourierRadialGenerator<double> SPFGenerator;
  SPFGenerator::Pointer spf = SPFGenerator::New();
  spf->SetScale(scale);
  spf->SetSPFType(SPFGenerator::SPF);
  spf->SetN(N);
  spf->SetL(L);
  double spfValue = spf->Evaluate(radius,false);

  utlPrintVar4(true, N, L, radius, spfValue );
  
  return 0;
}
