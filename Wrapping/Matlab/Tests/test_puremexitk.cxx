/**
 *       @file  test_puremex.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-09-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include <iostream>
#include "mex.h" 
#include "tr1/memory" 
#include <vector>
#include "itkImage.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  // utlGlobalException(nrhs!=1, "Bad number of inputs arguments");
  // double num = mxGetScalar(prhs[0]);
  std::tr1::shared_ptr<std::vector<int> > pp(new std::vector<int>());
  std::cout << "num = " << 15 << std::endl << std::flush;

  typedef itk::Image<double,3> ImageType;
  ImageType::Pointer image = ImageType::New();

}
