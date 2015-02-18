/**
 *       @file  SphericalPolarFourierImaging.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-08-2013
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

#include "itkSphericalPolarFourierImageFilter.h"
#include "itkDWIReader.h"
#include "itkGeneralizedHighOrderTensorImageFilter.h"
#include "itkFeaturesFromSPFImageFilter.h"
#include "itkProfileFromSPFImageFilter.h"
#include "itkODFFromSPFImageFilter.h"
#include "itkL1RegularizedLeastSquaresFISTASolver.h"

template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{
  utlGlobalException(!utl::mexCheckType<T>(prhs[0]),"type of argument 1 is not consistent");
  utlGlobalException(!utl::mexCheckType<T>(prhs[1]),"type of argument 2 is not consistent");
  utlGlobalException(!utl::mexCheckType<T>(prhs[2]),"type of argument 3 is not consistent");

  typedef itk::VectorImage<T, 3>  VectorImageType;
  typedef itk::Image<double, 3>  ImageType;
  typedef itk::SphericalPolarFourierEstimationImageFilter<VectorImageType, VectorImageType> SPFIFilterBaseType;

  utlException(mxGetNumberOfDimensions(prhs[0])!=4, "the input should have 4 dimension");

  const mwSize* dimsDWIs = mxGetDimensions(prhs[0]);
  int Nx = static_cast<int>(dimsDWIs[0]);
  int Ny = static_cast<int>(dimsDWIs[1]);
  int Nz = static_cast<int>(dimsDWIs[2]);
  int numberOfDWIs = static_cast<int>(dimsDWIs[3]);
  

  const mwSize* dimsOrientation = mxGetDimensions(prhs[1]);
  utlException(dimsOrientation[0]!=numberOfDWIs, "wrong number of gradients");
  utlException(dimsOrientation[1]!=3, "the column of gradient should be 3");

  const mwSize* dimsBVec = mxGetDimensions(prhs[2]);
  utlException(dimsBVec[0]!=numberOfDWIs, "wrong number of bVec");
  utlException(dimsBVec[1]!=1, "the column of bVec should be 1");

  double MD0 = utl::GetScalarStructDef<double>(prhs[3],"MD0",-1.0);
  double tau = utl::GetScalarStructDef<double>(prhs[3],"tau",ONE_OVER_4_PI_2);
  double scale = utl::GetScalarStructDef<double>(prhs[3],"scale",-1.0);
  int sh = utl::GetScalarStructDef<int>(prhs[3],"sh",-1);
  utlGlobalException(sh<=0, "need to set sh");
  int ra = utl::GetScalarStructDef<int>(prhs[3],"ra",-1);
  utlGlobalException(ra<=0, "need to set ra");
  // std::string basisType = utl::GetScalarStructDef<std::string>(prhs[3],"basisType","SPF");
  std::string estimation = utl::GetScalarStructDef<std::string>(prhs[3],"estimation","LS");
  std::string solver = utl::GetScalarStructDef<std::string>(prhs[3],"solver","SPAMS");
  double lambdaSH = utl::GetScalarStructDef<double>(prhs[3],"lambdaSH",0.0);
  double lambdaRA = utl::GetScalarStructDef<double>(prhs[3],"lambdaRA",0.0);
  double lambdaL1 = utl::GetScalarStructDef<double>(prhs[3],"lambdaL1",0.0);
  // int odfOrder = utl::GetScalarStructDef<double>(prhs[3],"odfOrder",2);
  // double radius = utl::GetScalarStructDef<double>(prhs[3],"radius",0.015);
  mxArray* mdImageArray = utl::GetArrayStruct(prhs[3], "mdImage" );
  double numericalB0Weight = utl::GetScalarStructDef<double>(prhs[3],"numericalB0Weight",-1.0);
  mxArray* dictionaryArray = utl::GetArrayStruct(prhs[3], "dictionary" );
  mxArray* energyArray = utl::GetArrayStruct(prhs[3], "energy" );
  double energyPower = utl::GetScalarStructDef<double>(prhs[3],"energyPower",1.0);
  int maxIter = utl::GetScalarStructDef<int>(prhs[3],"maxIter",1000);
  double minChange = utl::GetScalarStructDef<double>(prhs[3],"minChange",0.0001);
  mxArray* maskArray = utl::GetArrayStruct(prhs[3], "mask" );
  bool debug = utl::GetScalarStructDef<bool>(prhs[3],"debug",false);
  int thread = utl::GetScalarStructDef<double>(prhs[3],"thread",-1.0);


  typename SPFIFilterBaseType::Pointer spfiFilter=NULL;
  typedef itk::SphericalPolarFourierImageFilter<VectorImageType, VectorImageType> SPFIFilterType;
  spfiFilter = SPFIFilterType::New();
  std::cout << "Use SPF basis" << std::endl << std::flush;
  spfiFilter->SetIsOriginalBasis(true);

  if (maskArray)
    {
    typename ImageType::Pointer maskImage = ImageType::New();
    itk::GetITKImageFromMXArray(maskArray, maskImage);
    spfiFilter->SetMaskImage(maskImage);
    }
  typename VectorImageType::Pointer dwiImage = VectorImageType::New();
  itk::GetITKVectorImageFromMXArray(prhs[0], dwiImage);
  spfiFilter->SetInput(dwiImage);

  if (MD0>0)
    spfiFilter->SetMD0(MD0);
  //NOTE: set tau before spfiFilter->GetSamplingSchemeQSpace()->SetBVector(bVec) because it uses tau to convert b values to q values
  if (tau>0)
    spfiFilter->GetSamplingSchemeQSpace()->SetTau(tau);

  typename SPFIFilterBaseType::MatrixPointer grad(new typename SPFIFilterBaseType::MatrixType());
  utl::GetUtlMatrixFromMXArray(prhs[1], grad.get() );
  *grad = utl::CartesianToSpherical(*grad); // convert to spherical format
  spfiFilter->GetSamplingSchemeQSpace()->SetOrientationsSpherical(grad);

  typename SPFIFilterBaseType::STDVectorPointer bVec(new typename SPFIFilterBaseType::STDVectorType());
  utl::GetSTDVectorFromMXArray(prhs[2], bVec.get() );
  spfiFilter->GetSamplingSchemeQSpace()->SetBVector(bVec);
  spfiFilter->SetSHRank(sh);
  spfiFilter->SetRadialRank(ra);
  spfiFilter->SetBasisScale(scale);
    spfiFilter->SetIsAnalyticalB0(true);
  
  spfiFilter->SetLambdaSpherical(lambdaSH);
  spfiFilter->SetLambdaRadial(lambdaRA);
  spfiFilter->SetLambdaL1(lambdaL1);

  if (estimation=="LS")
    {
    spfiFilter->SetEstimationType(SPFIFilterBaseType::LS);
    }
  else if (estimation=="L1_2" || estimation=="L1_DL")
    {

    if (estimation=="L1_2")
      {
      utlGlobalException(lambdaSH<=0 && lambdaRA<=0, "need to set lambdaSH and lambdaRA when estimation=\"L1_2\".");
      spfiFilter->SetEstimationType(SPFIFilterBaseType::L1_2);
      }
    else if (estimation=="L1_DL")
      {
      utlGlobalException(lambdaL1<=0, "need to set lambdaL1 when estimation=\"L1_DL\".");
      spfiFilter->SetEstimationType(SPFIFilterBaseType::L1_DL);
      }

    if (solver=="FISTA_LS")
      {
      typedef itk::L1RegularizedLeastSquaresFISTASolver<double> L1SolverType;
      L1SolverType::Pointer l1Sol = L1SolverType::New();
      l1Sol->SetUseL2SolverForInitialization(solver=="FISTA_LS");
      l1Sol->SetMaxNumberOfIterations(maxIter);
      l1Sol->SetMinRelativeChangeOfCostFunction(minChange);
      l1Sol->SetMinRelativeChangeOfPrimalResidual(minChange);
      spfiFilter->SetL1FISTASolver(l1Sol);
      spfiFilter->SetL1SolverType(SPFIFilterBaseType::FISTA_LS);
      }
    if (solver=="SPAMS")
      {
      typedef itk::SpamsWeightedLassoSolver<double> L1SolverType;
      L1SolverType::Pointer l1Sol = L1SolverType::New();
      spfiFilter->SetL1SpamsSolver(l1Sol);
      spfiFilter->SetL1SolverType(SPFIFilterBaseType::SPAMS);
      }
    }
  else
    utlGlobalException(true, "wrong estimation type");

  if (mdImageArray)
    {
    typename ImageType::Pointer mdImage = ImageType::New();
    itk::GetITKImageFromMXArray(mdImageArray, mdImage);
    spfiFilter->SetMDImage(mdImage);
    }

  // spfiFilter->SetBasisEnergyPowerDL(energyPower);
  spfiFilter->SetBasisEnergyPowerDL(1.0);

  utl::InitializeThreadedLibraries(thread);
  if (thread>0)
    spfiFilter->SetNumberOfThreads(thread);
  spfiFilter->SetDebug(debug);

  std::cout << "SPF estimation starts" << std::endl << std::flush;
  spfiFilter->Update();
  std::cout << "SPF estimation ends" << std::endl << std::flush;
  
  typename VectorImageType::Pointer spf = spfiFilter->GetOutput();
  itk::GetMXArrayFromITKVectorImage(spf, plhs[0]);
  
  if (nlhs==2)
    {
    typename ImageType::Pointer scaleImage = spfiFilter->GetScaleImage();
    itk::GetMXArrayFromITKImage(scaleImage, plhs[1]);
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs!=4, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=1 && nlhs!=2, "Bad number of outputs arguments");

  if (mxGetClassID(prhs[0]) == mxSINGLE_CLASS) 
    callFunction<float>(plhs,prhs,nlhs,nrhs);
  else
    callFunction<double>(plhs,prhs,nlhs,nrhs);
} 
