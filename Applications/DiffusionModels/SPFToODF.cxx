/**
 *       @file  SPFToODFConverter.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-08-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "utl.h"
#include "SPFToODFCLP.h"
#include "itkODFFromSPFImageFilter.h"
#include "itkImage.h"
#include "itkSPFScaleFromMeanDiffusivityImageFilter.h"



/**
 * \brief  DWI/EAP profile (represented by SH basis) converted from SPF coefficients.
 */
int 
main (int argc, char const* argv[])
{

  // GenerateCLP
  PARSE_ARGS;
  
  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3>  VectorImageType;
  typedef itk::Image<PrecisionType, 3>  ScalarImageType;
    
  VectorImageType::Pointer spf=0;
  ScalarImageType::Pointer maskImage=0;
  itk::ReadImage<VectorImageType>(_InputSPFFile, spf);
  if (_MaskFileArg.isSet())
    itk::ReadImage<ScalarImageType>(_MaskFile, maskImage);

  ScalarImageType::Pointer mdImage=0, scaleImage=0;
  if (_MDImageFileArg.isSet())
    {
    itk::ReadImage<ScalarImageType>(_MDImageFile, mdImage);

    typedef itk::SPFScaleFromMeanDiffusivityImageFilter<ScalarImageType, ScalarImageType>  ScaleFromMDfilterType;
    ScaleFromMDfilterType::Pointer scaleFromMDfilter = ScaleFromMDfilterType::New();
    scaleFromMDfilter->SetMD0(_MD0);
    scaleFromMDfilter->SetTau(_Tau);
      scaleFromMDfilter->SetIsOriginalBasis(true);
    scaleFromMDfilter->SetInput(mdImage);
    scaleFromMDfilter->Update();
    scaleImage = scaleFromMDfilter->GetOutput();
    }

  typedef itk::FeaturesFromSPFImageFilter<VectorImageType, VectorImageType> FeaturesFromSPFFilterType;
  // FeaturesFromSPFFilterType::Pointer featureFromSPFFilter=NULL;

  typedef itk::ODFFromSPFImageFilter<VectorImageType, VectorImageType> ODFFromSPFFilterType;
  // ProfileFromSPFFilterType::Pointer featureFromSPFFilter = ProfileFromSPFFilterType::New();
  ODFFromSPFFilterType::Pointer featureFromSPFFilter = ODFFromSPFFilterType::New();
  if (_MaskFileArg.isSet())
    featureFromSPFFilter->SetMaskImage(maskImage);
  featureFromSPFFilter->SetSHRank(_SHRank);
  featureFromSPFFilter->SetRadialRank(_RadialRank);
  featureFromSPFFilter->SetMD0(_MD0);
  featureFromSPFFilter->SetTau(_Tau);
  featureFromSPFFilter->SetBasisScale(_Scale);
    featureFromSPFFilter->SetBasisType(FeaturesFromSPFFilterType::SPF);
  if (_MDImageFileArg.isSet())
    {
    featureFromSPFFilter->SetScaleImage(scaleImage);
    }
  if (_Debug)
    featureFromSPFFilter->DebugOn();
  featureFromSPFFilter->SetInput(spf);
  featureFromSPFFilter->SetBMax(_BMax);
  featureFromSPFFilter->SetODFOrder(_ODFOrder);
  if (_NumberOfThreads>0)
    featureFromSPFFilter->SetNumberOfThreads(_NumberOfThreads);
  std::cout << "ODF estimation starts" << std::endl << std::flush;
  featureFromSPFFilter->Update();
  std::cout << "ODF estimation ends" << std::endl << std::flush;
  VectorImageType::Pointer odf = featureFromSPFFilter->GetOutput();
  itk::SaveImage<VectorImageType>(odf, _OutputFile);

  return 0;
}

