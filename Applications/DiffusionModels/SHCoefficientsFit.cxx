/**
 *       @file  SHCoefficientsFit.cxx
 *      @brief  
 *     Created  "10-06-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "itkCommandProgressUpdate.h"
#include "itkSHCoefficientsFit.h"
#include "utl.h"

#include "SHCoefficientsFitCLP.h"
#include "itkUnaryFunctorVectorImageFilter.h"

/**
 * \brief  fit spherical function samples using SH basis
 */
int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;
  utl::LogLevel = _Verbose;

  utlGlobalException(_SHRank<0, "need to set _SHRank");

  typedef double PrecisionType; 
  typedef itk::VectorImage<PrecisionType, 3> ImageType;
  typedef itk::Functor::SHCoefficientsFit<double>  FunctorType;
  typedef itk::UnaryFunctorVectorImageFilter<ImageType, ImageType, FunctorType> UnaryFunctorFilterType;

  ImageType::Pointer sfImage = ImageType::New();
  itk::ReadImage<ImageType>(_InputFile, sfImage);

  FunctorType::MatrixPointer grad = utl::ReadGrad<double>(_OrientationFile,DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  
  FunctorType shFit;
  shFit.SetSHRank(_SHRank);
  if (_Lambda>0)
    shFit.SetLambda(_Lambda);
  shFit.SetOrientations(grad);
  shFit.SetPower(_Power);
  shFit.SetLogLevel(utl::LogLevel);

  typename UnaryFunctorFilterType::Pointer filter = UnaryFunctorFilterType::New();

  if (_MaskImageFile!="")
    filter->SetMaskImage(_MaskImageFile);
  filter->SetInput(sfImage);
  filter->SetFunctor(shFit);

  filter->SetDebug(_Verbose>=LOG_DEBUG);
  filter->SetLogLevel(utl::LogLevel);
  if (_NumberOfThreads>0)
    filter->SetNumberOfThreads(_NumberOfThreads);
  itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  if (_ShowProgressArg.isSet())
    filter->AddObserver( itk::ProgressEvent(), observer );

  filter->Update();

  ImageType::Pointer shImage = filter->GetOutput();
  itk::SaveImage<ImageType>(shImage, _OutputFile);
  
  return 0;
}
