/**
 *       @file  SPFToScalarMap.cxx
 *      @brief  calculate RTO map from SPF coefficients
 *     Created  "03-18-2014
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "SPFToScalarMapCLP.h"
#include "itkScalarMapFromSPFImageFilter.h"

#include "itkCommandProgressUpdate.h"

/**
 * \brief calculate RTO map from SPF coefficients  
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;
  
  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3>  VectorImageType;
  typedef itk::Image<PrecisionType, 3>  ScalarImageType;
    
  VectorImageType::Pointer spf=NULL;
  ScalarImageType::Pointer maskImage=NULL;
  itk::ReadImage<VectorImageType>(_InputSPFFile, spf);
  if (_MaskFileArg.isSet())
    itk::ReadImage<ScalarImageType>(_MaskFile, maskImage);

  ScalarImageType::Pointer mdImage= ScalarImageType::New();
  ScalarImageType::Pointer scaleImage= ScalarImageType::New();
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

  typedef itk::FeaturesFromSPFImageFilter<VectorImageType, ScalarImageType> FeaturesFromSPFFilterType;

  typedef itk::ScalarMapFromSPFImageFilter<VectorImageType, ScalarImageType> FilterType;
  FilterType::Pointer featureFromSPFFilter = FilterType::New();
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
  if (_MapType=="RTO")
    featureFromSPFFilter->SetMapType(FilterType::RTO);
  else if(_MapType=="MSD")
    featureFromSPFFilter->SetMapType(FilterType::MSD);
  else if(_MapType=="PFA")
    featureFromSPFFilter->SetMapType(FilterType::PFA);
  if (_Debug)
    featureFromSPFFilter->DebugOn();
  featureFromSPFFilter->SetInput(spf);
  if (_NumberOfThreads>0)
    featureFromSPFFilter->SetNumberOfThreads(_NumberOfThreads);

  itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  if (_ShowProgressArg.isSet())
    featureFromSPFFilter->AddObserver( itk::ProgressEvent(), observer );

  std::cout << "Scalar map estimation starts" << std::endl << std::flush;
  featureFromSPFFilter->Update();
  std::cout << "Scalar map estimation ends" << std::endl << std::flush;

  ScalarImageType::Pointer p0 = featureFromSPFFilter->GetOutput();
  itk::SaveImage<ScalarImageType>(p0, _OutputFile);
  
  return 0;
}
