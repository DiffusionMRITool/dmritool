/**
 *       @file  ReadGrad.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-05-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "utl.h"

#include "mex.h" 
#include "mexutils.h"
#include "utlMEX.h"

#include "itkSphericalPolarFourierGenerator.h"

template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{
  utlGlobalException(mxGetClassID(prhs[0])!=mxCHAR_CLASS, "the first input has to be a string");

  std::string homeStr = getenv ("HOME");
  std::string gradPath = homeStr + "/.dmritool/Data/";
  std::string gradFileStr;
  // std::cout << "gradPath = " << gradPath << std::endl << std::flush;

  std::string gradStr;
  utl::GetString(prhs[0], gradStr);

  typedef utl::NDArray<T,2> MatrixType;

  utl_shared_ptr<MatrixType> grad(new MatrixType());
  std::string ext, fileNoExt;
  utl::GetFileExtension(gradStr, ext, fileNoExt);

  if (ext=="txt")
    {
    gradFileStr = gradStr;
    }
  else if (ext=="")
    {
    utlGlobalException(nrhs!=2, "should have the second parameter");
    int num = mxGetScalar(prhs[1]);
    if (utl::StringCompareCaseIgnored(gradStr, "elec"))
      {
      if (num<10)
        gradFileStr = gradPath + "ElectricRepulsion/Elec00" + utl::ConvertNumberToString(num) + ".txt";
      else if (num<100)
        gradFileStr = gradPath + "ElectricRepulsion/Elec0" + utl::ConvertNumberToString(num) + ".txt";
      else
        gradFileStr = gradPath + "ElectricRepulsion/Elec" + utl::ConvertNumberToString(num) + ".txt";
      }
    else if (utl::StringCompareCaseIgnored(gradStr, "tess"))
      {
      gradFileStr = gradPath + "Tessellation/directions_t" + utl::ConvertNumberToString(num) + ".txt";
      }
    else
      utlGlobalException(true, "wrong mode");
    }
  else
    utlGlobalException(true, "wrong mode or gradient file name");

  grad = utl::ReadGrad<T>(gradFileStr,DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);

  utl::GetMXArrayFromUtlMatrix(grad.get(), plhs[0]);
}

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs!=1 && nrhs!=2, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=1, "Bad number of outputs arguments");

  // if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS) 
  //   callFunction<float>(plhs,prhs,nlhs,nrhs);
  // else
    callFunction<double>(plhs,prhs,nlhs,nrhs);
} 

