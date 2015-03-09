/**
 *       @file  DWIAddNoise.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-19-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */


#include "utl.h"
#include "itkAddNoiseToDWIImageFilter.h"

#include "DWINoiseGeneratorCLP.h"

/**
 * \brief  add Gaussian or Rician noise to DWI data
 */
int 
main (int argc, char const* argv[])
{
  
  // GenerateCLP
  PARSE_ARGS;

  // Time Probe
  itk::TimeProbe clock;
  
  // Input and Output Image Types
  typedef itk::VectorImage<double, 3> ImageType;
  
  typedef itk::Image<double, 3> MaskImageType;
  typedef itk::Image<double, 3> B0ImageType;

  ImageType::Pointer dwi;
  ImageType::Pointer dwiNoise;
  B0ImageType::Pointer b0;
  MaskImageType::Pointer mask;

  typedef itk::AddNoiseToDWIImageFilter<ImageType, B0ImageType, MaskImageType> AddNoiseFilterType;
  AddNoiseFilterType::Pointer addNoiseFilter = AddNoiseFilterType::New();

  itk::ReadVectorImage<double>(_InputFile, dwi);
  // utlPrintVar1(true, dwi->GetSpacing());
  addNoiseFilter->SetInput(dwi);
  if (_b0FileArg.isSet())
    {
    itk::ReadImage<B0ImageType>(_b0File, b0);
    // utlPrintVar1(true, b0->GetSpacing());
    addNoiseFilter->SetB0Image(b0);
    }
  if (_maskFileArg.isSet())
    {
    itk::ReadImage<MaskImageType>(_maskFile, mask);
    addNoiseFilter->SetMaskImage(mask);
    }

  utlException(_NoiseSigmaArg.isSet() && _SNRArg.isSet(), "Only one of these is allowed: --noisesigma or --snr");
  utlException(!_NoiseSigmaArg.isSet() && !_SNRArg.isSet(), "Set --noisesigma or --snr");
  if (_NoiseSigmaArg.isSet())
    addNoiseFilter->SetSigma(_NoiseSigma);
  if (_SNRArg.isSet())
    addNoiseFilter->SetSigma(1.0/_SNR);

  if (_Noisetype=="GAUSSIAN")
    addNoiseFilter->SetNoisetype(AddNoiseFilterType::GAUSSIAN);
  if (_Noisetype=="RICIAN")
    addNoiseFilter->SetNoisetype(AddNoiseFilterType::RICIAN);

  // Generate NoiseGenerate
  // addNoiseFilter->Print(std::cout<<"addNoiseFilter\n");
  addNoiseFilter->Update();
  
  dwiNoise = addNoiseFilter->GetOutput();
  itk::SaveImage<ImageType>(dwiNoise,_OutputFile);
      
  
  return EXIT_SUCCESS;
}
