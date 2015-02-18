/**
 *       @file  SHCoefficientsTOGFA.cxx
 *      @brief  
 *     Created  "02-11-2013
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "utl.h"
#include "itkSHCoefficientsToGFAImaeFilter.h"
#include "SHCoefficientsToGFACLP.h"

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

  typedef itk::SHCoefficientsToGFAImaeFilter<VectorImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput(shImage);
  filter->Update();
  ImageType::Pointer gfaImage;
  gfaImage = filter->GetOutput();

  itk::SaveImage(gfaImage, _OutputFile);
  
  return 0;
}
