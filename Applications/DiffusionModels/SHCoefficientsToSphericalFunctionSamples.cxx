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

#include "itkFunctors.h"
#include "itkUnaryFunctorVectorImageFilter.h"

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

  ImageType::Pointer outputImage, sfImage;
  sfImage = filter->GetOutput();
  if (std::fabs(_Power-1.0)>1e-10)
    {
    itk::Functor::RPOWER<double, double> p0;
    p0.SetArgument(_Power);
    utl::Functor::ScalarFunctorWrapper<itk::Functor::RPOWER<double,double>> pp(p0);
    itk::UnaryVectorOPImage<ImageType, ImageType>(sfImage, outputImage, pp);
    }
  else
    outputImage = sfImage;

  itk::SaveImage(outputImage, _OutputFile);

  return 0;
}
