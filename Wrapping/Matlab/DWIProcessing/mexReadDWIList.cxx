/**
 *       @file  ReadDWIList.cxx
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
#include "utlMEX.h"

#include "itkDWIReader.h"

template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{
  utlGlobalException(mxGetClassID(prhs[0])!=mxCHAR_CLASS, "the first input has to be a string");
  utlGlobalException(nrhs==2 && mxGetClassID(prhs[1])!=mxSTRUCT_CLASS, "the second input has to be a struct");

  std::string dwiFile;
  utl::GetString(prhs[0], dwiFile);


  bool normalize = utl::GetScalarStructDef<bool>(nrhs==2?prhs[1]:NULL,"normalize",true);
  bool correctDWI = utl::GetScalarStructDef<bool>(nrhs==2?prhs[1]:NULL,"correctDWI",true);
  bool warn = utl::GetScalarStructDef<bool>(nrhs==2?prhs[1]:NULL,"warn",true);
  double bThreshold = utl::GetScalarStructDef<bool>(nrhs==2?prhs[1]:NULL,"bThreshold",-1);

  // utlPrintVar3(true, normalize, correctDWI, bThreshold);
  typedef itk::DWIReader<T,3>  DWIReaderType;
  typename DWIReaderType::Pointer dwiReader = DWIReaderType::New();

  dwiReader->SetConfigurationFile(dwiFile);
  dwiReader->SetNormalizeDWI(normalize);
  if (bThreshold>0)
    {
    dwiReader->GetSamplingSchemeQSpace()->SetBThresholdSingleShell(bThreshold);
    }
  dwiReader->SetCorrectDWIValues(correctDWI);
  dwiReader->SetShowWarnings(warn);

  typedef itk::Image<T,3> Image3DType;
  typename Image3DType::Pointer b0Image = Image3DType::New();
  mxArray* b0Array = utl::GetArrayStruct(prhs[1], "b0Image" );
  if (b0Array)
    {
    itk::GetITKImageFromMXArray(b0Array, b0Image);
    // itk::PrintImage(b0Image, "b0");
    dwiReader->SetB0Image(b0Image);
    }

  dwiReader->Update();

  typedef itk::VectorImage<T,3> VectorImageType;
  typename VectorImageType::Pointer dwiImage = dwiReader->GetOutput();
  // itk::PrintVectorImage(dwiImage, "dwiImage");
  itk::GetMXArrayFromITKVectorImage(dwiImage, plhs[0]);

  utl_shared_ptr<std::vector<double> > bVec = dwiReader->GetSamplingSchemeQSpace()->GetBVector();
  utl::GetMXArrayFromSTDVector<double>(bVec.get(), plhs[1]);

  utl_shared_ptr<utl::NDArray<double,2> >  grad = dwiReader->GetSamplingSchemeQSpace()->GetOrientationsCartesian();
  utl::GetMXArrayFromUtlMatrix<double>(grad.get(), plhs[2]);

  if (nlhs>3)
    {
    typename Image3DType::Pointer b0Image = dwiReader->GetB0Image();
    // itk::PrintImage(b0Image, "b02");
    itk::GetMXArrayFromITKImage(b0Image, plhs[3]);
    }
}
  
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utl_shared_ptr<std::vector<double> > grad(new std::vector<double>());
  utlGlobalException(nrhs!=1 && nrhs!=2, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=3 && nlhs!=4, "Bad number of outputs arguments");

  // utlPrintVar1(true, mxGetClassID(prhs[0]));
  // utlPrintVar2(true, mxDOUBLE_CLASS, mxSINGLE_CLASS);

  // if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS) 
  //   callFunction<float>(plhs,prhs,nlhs,nrhs);
  // else
    callFunction<double>(plhs,prhs,nlhs,nrhs);

}
