/**
 *       @file  GetSamplesFromSPFCoefficients.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-11-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "mex.h" 
#include "utl.h"
#include "utlMEX.h"

#include "itkProfileFromSPFImageFilter.h"
#include "itkSPFScaleFromMeanDiffusivityImageFilter.h"


template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{
  utlGlobalException(!utl::mexCheckType<T>(prhs[0]),"type of argument 1 is not consistent");
  utlGlobalException(!utl::mexCheckType<T>(prhs[1]),"type of argument 2 is not consistent");

  typedef itk::VectorImage<T, 3>  VectorImageType;
  typedef itk::Image<double, 3>  ImageType;
  typedef utl::NDArray<T,2> MatrixType;

  typedef itk::ProfileFromSPFImageFilter<VectorImageType, VectorImageType> ProfileFromSPFFilterType;
  typename ProfileFromSPFFilterType::Pointer profEstimator = ProfileFromSPFFilterType::New();

  utlGlobalException(mxGetNumberOfDimensions(prhs[0])!=4, "the input should have 4 dimension");

  const mwSize* dimsDWIs = mxGetDimensions(prhs[0]);
  int Nx = static_cast<int>(dimsDWIs[0]);
  int Ny = static_cast<int>(dimsDWIs[1]);
  int Nz = static_cast<int>(dimsDWIs[2]);
  int numberOfDWIs = static_cast<int>(dimsDWIs[3]);

  typename VectorImageType::Pointer dwiImage = VectorImageType::New();
  itk::GetITKVectorImageFromMXArray(prhs[0], dwiImage);
  profEstimator->SetInput(dwiImage);

  // utlGlobalException(nrhs==3 && mxGetNumberOfDimensions(prhs[1])!=1, "radius should be a scalar value. mxGetNumberOfDimensions(prhs[1])="<< mxGetNumberOfDimensions(prhs[1]));
  if (nrhs==3)
    {
    double radius = mxGetScalar(prhs[1]);
    profEstimator->SetRadius(radius);
    }

  const mxArray* params = nrhs==3?prhs[2]:prhs[3];

  double MD0 = utl::GetScalarStructDef<double>(params,"MD0",-1.0);
  double tau = utl::GetScalarStructDef<double>(params,"tau",ONE_OVER_4_PI_2);
  double scale = utl::GetScalarStructDef<double>(params,"scale",-1.0);
  int sh = utl::GetScalarStructDef<int>(params,"sh",-1);
  utlGlobalException(sh<=0, "need to set sh");
  int ra = utl::GetScalarStructDef<int>(params,"ra",-1);
  utlGlobalException(ra<=0, "need to set ra");
  std::string basisType = utl::GetScalarStructDef<std::string>(params,"basisType","SPF");
  double radius = utl::GetScalarStructDef<double>(params,"radius",0.015);
  mxArray* scaleImageArray = utl::GetArrayStruct(params, "scaleImage" );
  mxArray* maskArray = utl::GetArrayStruct(params, "mask" );
  bool fourier = utl::GetScalarStructDef<bool>(params,"fourier",false);
  bool inqspace = utl::GetScalarStructDef<bool>(params,"inqspace",true);

  bool debug = utl::GetScalarStructDef<bool>(params,"debug",false);
  int thread = utl::GetScalarStructDef<double>(params,"thread",-1);

  if (nrhs==4)
    {
    typename ProfileFromSPFFilterType::STDVectorPointer radiusVec(new typename ProfileFromSPFFilterType::STDVectorType());
    utl::GetSTDVectorFromMXArray(prhs[1], radiusVec.get() );
    profEstimator->SetRadiusVector(radiusVec);

    utl_shared_ptr<MatrixType> grad(new MatrixType());
    utl::GetUtlMatrixFromMXArray(prhs[2], grad.get() );
    *grad = utl::CartesianToSpherical(*grad); // convert to spherical format
    profEstimator->SetOrientations(grad);
    }
  profEstimator->SetSHRank(sh);
  profEstimator->SetRadialRank(ra);
  if (MD0>0)
    profEstimator->SetMD0(MD0);
  if (tau>0)
    profEstimator->SetTau(tau);
  profEstimator->SetBasisScale(scale);
  profEstimator->SetIsFourier(fourier);
  profEstimator->SetIsInQSpace(inqspace);
  
  if (basisType=="SPF")
    {
    profEstimator->SetBasisType(ProfileFromSPFFilterType::SPF);
    }
  else if (basisType=="DSPF")
    profEstimator->SetBasisType(ProfileFromSPFFilterType::DSPF);

  if (scaleImageArray)
    {
    typename ProfileFromSPFFilterType::ScalarImagePointer scaleImage = ProfileFromSPFFilterType::ScalarImageType::New();
    itk::GetITKImageFromMXArray(scaleImageArray, scaleImage);
    profEstimator->SetScaleImage(scaleImage);
    }

  if (maskArray)
    {
    typename ProfileFromSPFFilterType::ScalarImagePointer maskImage = ProfileFromSPFFilterType::ScalarImageType::New();
    itk::GetITKImageFromMXArray(maskArray, maskImage);
    profEstimator->SetMaskImage(maskImage);
    }

  profEstimator->SetDebug(debug);
  utl::InitializeThreadedLibraries(thread);
  if (thread>0)
    profEstimator->SetNumberOfThreads(thread);

  // profEstimator->Print(std::cout<<"profEstimator=");
  profEstimator->Update();

  typename VectorImageType::Pointer samples = profEstimator->GetOutput();
  itk::GetMXArrayFromITKVectorImage(samples, plhs[0]);
}

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs!=3 && nrhs!=4, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=1, "Bad number of outputs arguments");

  // if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS) 
  //   callFunction<float>(plhs,prhs,nlhs,nrhs);
  // else
    callFunction<double>(plhs,prhs,nlhs,nrhs);
} 
