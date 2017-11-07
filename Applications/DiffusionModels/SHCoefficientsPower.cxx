/**
 *       @file  SHCoefficientsPower.cxx
 *      @brief  
 *     Created  "02-06-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "SHCoefficientsPowerCLP.h"

#include "itkCommandProgressUpdate.h"
#include "utl.h"
#include "itkSHCoefficientsPowerImageFilter.h"

/**
 * \brief In each voxel,  calculate a power function of a given spherical function. 
 * Both input function and output function are represented using SH coefficients. 
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;
  
  utlGlobalException(_Power==1, "need to set power which is not 1");

  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3> ImageType;

  typedef itk::SHCoefficientsPowerImageFilter<ImageType, ImageType>  FilterType;
  FilterType::Pointer filter = FilterType::New();

  ImageType::Pointer shImage = ImageType::New();
  itk::ReadImage<ImageType>(_InputFile, shImage);

  filter->SetInput(shImage);
  filter->SetPower(_Power);
  filter->SetSHRank(_SHRank>0 ? _SHRank : utl::DimToRankSH(shImage->GetNumberOfComponentsPerPixel()));
  
  itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  if (_ShowProgressArg.isSet())
    filter->AddObserver( itk::ProgressEvent(), observer );
  if (_Debug)
    filter->DebugOn();
  if (_NumberOfThreads>0)
    filter->SetNumberOfThreads(_NumberOfThreads);
  
  filter->Update();

  ImageType::Pointer shPowerImage = filter->GetOutput();
  itk::SaveImage<ImageType>(shPowerImage, _OutputFile);
  
  return 0;
}
