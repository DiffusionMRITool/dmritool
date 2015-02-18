/**
 *       @file  MeanDiffusivityEstimator.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-26-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "MeanDiffusivityEstimatorCLP.h"
#include "utl.h"
#include "itkGeneralizedHighOrderTensorImageFilter.h"
#include "itkDWIReader.h"
#include "itkMeanDiffusivityFromGHOTImageFilter.h"
// #include "itkCommandProgressUpdate.h"
// #include "itkSPFScaleFromMeanDiffusivityImageFilter.h"




/**
 * \brief  estimate mean diffusivity in each voxel
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3>  VectorImageType;
  typedef itk::Image<PrecisionType, 3>  ImageType;
  
  typedef itk::DWIReader<PrecisionType>  DWIReaderType;
  typedef itk::GeneralizedHighOrderTensorImageFilter<VectorImageType> GHOTFilterType;
  
  DWIReaderType::Pointer reader = DWIReaderType::New();
  // reader->GetSamplingSchemeQSpace()->SetTau(_Tau);
  reader->SetConfigurationFile(_InputFile);
  reader->Update();

  GHOTFilterType::Pointer ghot = GHOTFilterType::New();
  typedef  GHOTFilterType::MaskImageType MaskImageType;
  MaskImageType::Pointer maskImage;
  if (_MaskFileArg.isSet())
    {
    itk::ReadImage(_MaskFile, maskImage);
    ghot->SetMaskImage(maskImage);
    }
  ghot->SetInput(reader->GetOutput());
  ghot->SetSamplingSchemeQSpace(reader->GetSamplingSchemeQSpace());
  ghot->SetSHRank(_SHRank);
  ghot->SetRadialRank(_RadialRank);
  if (_NumberOfThreads>0)
    ghot->SetNumberOfThreads( _NumberOfThreads );
  // itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  // ghot->AddObserver( itk::ProgressEvent(), observer );
  if (_DebugArg.isSet())
    ghot->DebugOn();
  std::cout << "Estimate MD using GHOT model" << std::endl << std::flush;
  ghot->Update();
  std::cout << "DONE" << std::endl << std::flush;

  // VectorImageType::Pointer outputImage = ghot->GetOutput();
  // outputImage->Print(std::cout<<"outputImage");
  // itk::PrintVectorImage<PrecisionType,3>(outputImage, "ghot");

  typedef itk::MeanDiffusivityFromGHOTImageFilter<VectorImageType, ImageType> MDFromGHOTFilterType;
  MDFromGHOTFilterType::Pointer mdFromGHOT = MDFromGHOTFilterType::New();
  mdFromGHOT->SetInput(ghot->GetOutput());
  mdFromGHOT->SetTau(ghot->GetSamplingSchemeQSpace()->GetTau());
  mdFromGHOT->SetScale(ghot->GetBasisScale());
  mdFromGHOT->Update();
  ImageType::Pointer md = mdFromGHOT->GetOutput();

  // md->Print(std::cout<<"md");
  // itk::PrintImage<PrecisionType,3>(md, "md");
  itk::SaveImage(md, _OutputFile);


  return 0;
}
