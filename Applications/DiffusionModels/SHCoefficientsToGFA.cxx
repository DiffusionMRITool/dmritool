/**
 *       @file  SHCoefficientsTOGFA.cxx
 *      @brief  
 *     Created  "02-11-2013
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "SHCoefficientsToGFACLP.h"

#include "utl.h"
#include "itkSHCoefficientsToGFAImageFilter.h"
#include "itkSHCoefficientsPowerImageFilter.h"

/**
 * \brief  Calculate GFA from SH coefficients.
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  typedef double  TScalarType;
  typedef itk::VectorImage<TScalarType,3> VectorImageType;
  typedef itk::Image<TScalarType,3> ImageType;

  VectorImageType::Pointer shImage = VectorImageType::New();
  itk::ReadVectorImage(_InputSHFile, shImage);

  typedef itk::SHCoefficientsToGFAImageFilter<VectorImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  if (_Power==1)
    {
    filter->SetInput(shImage);
    }
  else
    {
    typedef itk::SHCoefficientsPowerImageFilter<VectorImageType, VectorImageType>  FilterType;
    FilterType::Pointer shFilter = FilterType::New();


    shFilter->SetInput(shImage);
    shFilter->SetPower(_Power);
    shFilter->SetSHRank(utl::DimToRankSH(shImage->GetNumberOfComponentsPerPixel()));

    shFilter->Update();

    VectorImageType::Pointer shPowerImage = shFilter->GetOutput();
    filter->SetInput(shPowerImage);
    }

  filter->Update();
  ImageType::Pointer gfaImage = filter->GetOutput();

  itk::SaveImage(gfaImage, _OutputFile);
  
  return 0;
}
