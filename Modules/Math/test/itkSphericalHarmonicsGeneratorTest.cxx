/**
 *       @file  itkSphericalHarmonicsGeneratorTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "03-25-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkSphericalHarmonicsGenerator.h"

/**
 * \brief  test itkSphericalHarmonicsGenerator
 */
int 
main (int argc, char const* argv[])
{
  std::cout << "RealTripleIntegration(6,-2,4,0,2,-2,false) = " << itk::SphericalHarmonicsGenerator<double>::RealTripleIntegration(6,-2,4,0,2,-2, false) << std::endl << std::flush;
  std::cout << "RealTripleIntegration(6,-2,4,0,2,-2,true) = " << itk::SphericalHarmonicsGenerator<double>::RealTripleIntegration(6,-2,4,0,2,-2, true) << std::endl << std::flush;
  return 0;
}
