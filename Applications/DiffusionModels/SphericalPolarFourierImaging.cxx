/**
 *       @file  SphericalPolarFourierImaging.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-04-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "SphericalPolarFourierImagingCLP.h"

#include "itkL1RegularizedLeastSquaresFISTASolver.h"

#include "itkSphericalPolarFourierImageFilter.h"
#include "itkDWIReader.h"
#include "itkGeneralizedHighOrderTensorImageFilter.h"
#include "itkFeaturesFromSPFImageFilter.h"
#include "itkProfileFromSPFImageFilter.h"
#include "itkODFFromSPFImageFilter.h"
#include "itkSpamsWeightedLassoSolver.h"
#include "itkSamplingSchemeQSpace.h"
#include "itkSamplingScheme3D.h"

#include "itkCommandProgressUpdate.h"

#include "utl.h"

/**
 * \brief  general SPFI using SPF basis
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;
  
  utl::LogLevel = _Verbose;
  bool _Debug = _Verbose>=LOG_DEBUG;

  utlGlobalException(!_OutputEAPProfileFileArg.isSet() && !_OutputODFFileArg.isSet() && !_OutputSPFFileArg.isSet(), "no output");

  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3>  VectorImageType;
  typedef itk::Image<PrecisionType, 3>  ImageType;

  typedef itk::DWIReader<PrecisionType>  DWIReaderType;
  
  ImageType::Pointer mdImage=NULL, maskImage=NULL;
  if (_MDImageFileArg.isSet())
    itk::ReadImage<ImageType>(_MDImageFile, mdImage);
  if (_MaskFileArg.isSet())
    itk::ReadImage<ImageType>(_MaskFile, maskImage);

  // read DWI, mdImage, mask
  DWIReaderType::Pointer reader = DWIReaderType::New();
  reader->GetSamplingSchemeQSpace()->SetTau(_Tau);
  reader->SetConfigurationFile(_InputFile);
  if (_MaskFileArg.isSet())
    reader->SetMaskImage(maskImage);
  reader->Update();



  // SPFI
  VectorImageType::Pointer spf;
  typedef itk::SphericalPolarFourierEstimationImageFilter<VectorImageType, VectorImageType> SPFIFilterBaseType;
  SPFIFilterBaseType::Pointer spfiFilter=NULL;
    typedef itk::SphericalPolarFourierImageFilter<VectorImageType, VectorImageType> SPFIFilterType;
    spfiFilter = SPFIFilterType::New();
    std::cout << "Use SPF basis" << std::endl << std::flush;

  utl::InitializeThreadedLibraries(_NumberOfThreads);
  if (_NumberOfThreads>0)
    spfiFilter->SetNumberOfThreads(_NumberOfThreads);
  if (_MaskFileArg.isSet())
    spfiFilter->SetMaskImage(maskImage);
  spfiFilter->SetInput(reader->GetOutput());
  spfiFilter->SetSamplingSchemeQSpace(reader->GetSamplingSchemeQSpace());
  spfiFilter->SetSHRank(_SHRank);
  spfiFilter->SetRadialRank(_RadialRank);
  spfiFilter->SetMD0(_MD0);
  spfiFilter->SetBasisScale(_Scale);

  spfiFilter->SetIsAnalyticalB0(true);

  if (_EstimationType=="LS")
    spfiFilter->SetEstimationType(SPFIFilterBaseType::LS);
  else if (_EstimationType=="L1_2")
    {
    utlGlobalException(_LambdaSpherical<=0 && _LambdaRadial<=0, "need to set _LambdaSpherical and _LambdaRadial when _EstimationType=\"L1_2\".");
    spfiFilter->SetEstimationType(SPFIFilterBaseType::L1_2);
    }
  else if (_EstimationType=="L1_DL")
    {
    utlGlobalException(_LambdaL1<0, "need to set _LambdaL1 when _EstimationType=\"L1_DL\".");
    spfiFilter->SetEstimationType(SPFIFilterBaseType::L1_DL);
    }
  spfiFilter->SetLambdaSpherical(_LambdaSpherical);
  spfiFilter->SetLambdaRadial(_LambdaRadial);
  spfiFilter->SetLambdaL1(_LambdaL1);

  if (_SolverType=="FISTA_LS")
    {
    typedef itk::L1RegularizedLeastSquaresFISTASolver<double> L1SolverType;
    L1SolverType::Pointer l1Solver = L1SolverType::New();
    l1Solver->SetUseL2SolverForInitialization(_SolverType=="FISTA_LS");
    l1Solver->SetMaxNumberOfIterations(_MaxIter);
    l1Solver->SetMinRelativeChangeOfCostFunction(_MinChange);
    l1Solver->SetMinRelativeChangeOfPrimalResidual(_MinChange);
    spfiFilter->SetL1FISTASolver(l1Solver);
    spfiFilter->SetL1SolverType(SPFIFilterBaseType::FISTA_LS);
    }
  if (_SolverType=="SPAMS")
    {
    typedef itk::SpamsWeightedLassoSolver<double> L1SolverType;
    L1SolverType::Pointer l1Solver = L1SolverType::New();
    spfiFilter->SetL1SpamsSolver(l1Solver);
    spfiFilter->SetL1SolverType(SPFIFilterBaseType::SPAMS);
    }


  typedef SPFIFilterBaseType::SamplingSchemeQSpaceType  SamplingSchemeQSpaceType;
  typedef SPFIFilterBaseType::SamplingSchemeRSpaceType  SamplingSchemeRSpaceType;
  SamplingSchemeQSpaceType::Pointer samplingQSpace = SamplingSchemeQSpaceType::New();
  SamplingSchemeRSpaceType::Pointer samplingRSpace = SamplingSchemeRSpaceType::New();

  if (_MDImageFileArg.isSet())
    spfiFilter->SetMDImage(mdImage);

  spfiFilter->SetBasisEnergyPowerDL(1.0);

  itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  if (_ShowProgressArg.isSet())
    spfiFilter->AddObserver( itk::ProgressEvent(), observer );
  spfiFilter->SetDebug(_Verbose>=LOG_DEBUG);
  spfiFilter->SetLogLevel(utl::LogLevel);

  std::cout << "SPF estimation starts" << std::endl << std::flush;
  spfiFilter->Update();
  std::cout << "SPF estimation ends" << std::endl << std::flush;
  spf = spfiFilter->GetOutput();
  // save SPF coefficients if needed 
  if (_OutputSPFFileArg.isSet())
    {
    itk::SaveImage<VectorImageType>(spf, _OutputSPFFile);
    }

  // output features
  typedef itk::FeaturesFromSPFImageFilter<VectorImageType, VectorImageType> FeaturesFromSPFFilterType;
  // FeaturesFromSPFFilterType::Pointer featureFromSPFFilter=0;

  // output EAP 
  if (_OutputEAPProfileFileArg.isSet())
    {
    typedef itk::ProfileFromSPFImageFilter<VectorImageType, VectorImageType> ProfileFromSPFFilterType;
    ProfileFromSPFFilterType::Pointer featureFromSPFFilter = ProfileFromSPFFilterType::New();
    if (_MaskFileArg.isSet())
      featureFromSPFFilter->SetMaskImage(maskImage);
    if (_NumberOfThreads>0)
      featureFromSPFFilter->SetNumberOfThreads(_NumberOfThreads);
    featureFromSPFFilter->SetSHRank(_SHRank);
    featureFromSPFFilter->SetRadialRank(_RadialRank);
    featureFromSPFFilter->SetMD0(spfiFilter->GetMD0());
    featureFromSPFFilter->SetTau(spfiFilter->GetSamplingSchemeQSpace()->GetTau());
    featureFromSPFFilter->SetBasisScale(spfiFilter->GetBasisScale());
      featureFromSPFFilter->SetBasisType(FeaturesFromSPFFilterType::SPF);
    if (_MDImageFileArg.isSet())
      {
      featureFromSPFFilter->SetScaleImage(spfiFilter->GetScaleImage());
      }
    if (_Debug)
      featureFromSPFFilter->DebugOn();
    featureFromSPFFilter->SetInput(spf);
    featureFromSPFFilter->SetRadius(_Radius);
    featureFromSPFFilter->SetIsFourier(true);
    featureFromSPFFilter->SetIsInQSpace(true);
    std::cout << "EAP profile estimation starts" << std::endl << std::flush;
    featureFromSPFFilter->Update();
    std::cout << "EAP profile estimation ends" << std::endl << std::flush;
    VectorImageType::Pointer eap = featureFromSPFFilter->GetOutput();
    itk::SaveImage<VectorImageType>(eap, _OutputEAPProfileFile);
    }

  // output ODF
  if (_OutputODFFileArg.isSet())
    {
    typedef itk::ODFFromSPFImageFilter<VectorImageType, VectorImageType> ODFFromSPFFilterType;
    ODFFromSPFFilterType::Pointer featureFromSPFFilter = ODFFromSPFFilterType::New();
    if (_MaskFileArg.isSet())
      featureFromSPFFilter->SetMaskImage(maskImage);
    if (_NumberOfThreads>0)
      featureFromSPFFilter->SetNumberOfThreads(_NumberOfThreads);
    featureFromSPFFilter->SetSHRank(_SHRank);
    featureFromSPFFilter->SetRadialRank(_RadialRank);
    featureFromSPFFilter->SetMD0(spfiFilter->GetMD0());
    featureFromSPFFilter->SetTau(spfiFilter->GetSamplingSchemeQSpace()->GetTau());
    featureFromSPFFilter->SetBasisScale(spfiFilter->GetBasisScale());
      featureFromSPFFilter->SetBasisType(FeaturesFromSPFFilterType::SPF);
    if (_MDImageFileArg.isSet())
      {
      featureFromSPFFilter->SetScaleImage(spfiFilter->GetScaleImage());
      }
    if (_Debug)
      featureFromSPFFilter->DebugOn();
    featureFromSPFFilter->SetInput(spf);
    featureFromSPFFilter->SetODFOrder(_ODFOrder);
    featureFromSPFFilter->SetBMax(-1.0);
    featureFromSPFFilter->SetIsFourier(true);
    featureFromSPFFilter->SetIsInQSpace(true);
    std::cout << "ODF estimation starts" << std::endl << std::flush;
    featureFromSPFFilter->Update();
    std::cout << "ODF estimation ends" << std::endl << std::flush;
    VectorImageType::Pointer odf = featureFromSPFFilter->GetOutput();
    itk::SaveImage<VectorImageType>(odf, _OutputODFFile);
    }

  return 0;
}

