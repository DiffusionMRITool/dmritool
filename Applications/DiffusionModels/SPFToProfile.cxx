/**
 *       @file  SPFToProfileConverter.cxx
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
#include "SPFToProfileCLP.h"
#include "itkProfileFromSPFImageFilter.h"
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

  typedef itk::ProfileFromSPFImageFilter<VectorImageType, VectorImageType> ProfileFromSPFFilterType;
  ProfileFromSPFFilterType::Pointer featureFromSPFFilter = ProfileFromSPFFilterType::New();
  // featureFromSPFFilter = ProfileFromSPFFilterType::New();
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
  featureFromSPFFilter->SetIsInQSpace(!_rSpaceArg.isSet());
  if (_OrientationsFileArg.isSet())
    {
    ProfileFromSPFFilterType::MatrixPointer grad = utl::ReadGrad<double>(_OrientationsFile, DIRECTION_NODUPLICATE, CARTESIAN_TO_SPHERICAL);
    featureFromSPFFilter->SetOrientations(grad); 
    }
  if (_RadiusVectorFileArg.isSet())
    {
    ProfileFromSPFFilterType::STDVectorPointer radiusVec(new ProfileFromSPFFilterType::STDVectorType());
    utl::ReadVector(_RadiusVectorFile, *radiusVec);
    featureFromSPFFilter->SetRadiusVector(radiusVec); 
    }
  if (_NumberOfThreads>0)
    featureFromSPFFilter->SetNumberOfThreads(_NumberOfThreads);
  if (_DebugArg.isSet())
    featureFromSPFFilter->DebugOn();
  featureFromSPFFilter->SetInput(spf);
  featureFromSPFFilter->SetRadius(_Radius);
  featureFromSPFFilter->SetIsFourier(_IsFourier);
  std::cout << "EAP profile estimation starts" << std::endl << std::flush;
  featureFromSPFFilter->Update();
  std::cout << "EAP profile estimation ends" << std::endl << std::flush;
  VectorImageType::Pointer eap = featureFromSPFFilter->GetOutput();
  itk::SaveImage<VectorImageType>(eap, _OutputFile);

  return 0;
}
