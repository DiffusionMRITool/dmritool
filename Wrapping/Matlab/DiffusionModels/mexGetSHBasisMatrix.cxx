/**
 *       @file  GetSHBasisMatrix.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-04-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include <cstdio>
#include <iostream>
#include <cstdlib>

#include "mex.h" 
#include "utl.h"
#include "mexutils.h"
#include "utlMEX.h"


template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{
  utlGlobalException(!utl::mexCheckType<T>(prhs[1]),"type of argument 1 is not consistent");

  int order = mxGetScalar(prhs[0]);

  const mwSize* dimsOrientation = mxGetDimensions(prhs[1]);
  int row = static_cast<int>(dimsOrientation[0]);
  int column = static_cast<int>(dimsOrientation[1]);

  utlException(column!=3,"orientation matrix should have 3 columns (x,y,z)");

  utl::NDArray<T,2> orientationMatrix(row, column);
  utl::GetUtlMatrixFromMXArray(prhs[1], &orientationMatrix);

  std::string mode = "cartesian";
  if (nrhs==3)
    {
    utl::GetString(prhs[2], mode);
    }
    
  utlException(mode!="cartesian" && mode!="spherical","type of argument 3 is not consistent");

  utl_shared_ptr<utl::NDArray<T,2> > BMatrix (new utl::NDArray<T,2>());
  if (mode=="spherical")
    {
    BMatrix = utl::ComputeSHMatrix(order, orientationMatrix, SPHERICAL_TO_SPHERICAL);
    }
  else
    {
    BMatrix = utl::ComputeSHMatrix(order, orientationMatrix, CARTESIAN_TO_SPHERICAL); 
    }

  utl::GetMXArrayFromUtlMatrix(BMatrix.get(), plhs[0]);

  return;
}

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs != 2 && nrhs != 3,"Bad number of inputs arguments");
  utlGlobalException(nlhs != 1,"Bad number of outputs arguments");

  // if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS) 
  //   callFunction<float>(plhs,prhs,nlhs,nrhs);
  // else
    callFunction<double>(plhs,prhs,nlhs,nrhs);
} 

