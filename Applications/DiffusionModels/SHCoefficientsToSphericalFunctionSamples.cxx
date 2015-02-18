/**
 *       @file  SHCoefficientsToSphericalFunctionSamples.cxx
 *      @brief  
 *     Created  "05-11-2013
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "utl.h"
#include "itkMultiplyByConstantMatrixVectorImageFilter.h"
#include "SHCoefficientsToSphericalFunctionSamplesCLP.h"

/**
 * \brief  calculate samples of a spherical function from its SH coefficients
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;
  
  typedef double  TScalarType;
  typedef utl::NDArray<TScalarType,2> MatrixType;
  typedef itk::VectorImage<TScalarType,3> ImageType;

  ImageType::Pointer inputImage = ImageType::New();
  itk::ReadVectorImage(_InputSHFile, inputImage);

  unsigned shDim = inputImage->GetNumberOfComponentsPerPixel();
  unsigned shRank = utl::DimToRankSH(shDim);

  utl_shared_ptr<MatrixType> grad = utl::ReadGrad<TScalarType>(_dataOrientationsFile, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  utl_shared_ptr<MatrixType> basisMatrix = utl::ComputeSHMatrix(shRank, *grad, CARTESIAN_TO_SPHERICAL);

  if (_Debug)
    std::cout << "basisMatrix: \n" << *basisMatrix;

  typedef itk::MultiplyByConstantMatrixVectorImageFilter<ImageType, MatrixType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput(inputImage);
  filter->SetConstantMatrix(*basisMatrix);
  filter->Update();

  ImageType::Pointer outputImage = filter->GetOutput();
  itk::SaveImage(outputImage, _OutputFile);

  return 0;
}
