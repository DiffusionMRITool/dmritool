/**
 *       @file  TensorFormatConverter.cxx
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "TensorFormatConverterCLP.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "utl.h"

/**
 * \brief  Convert different tensor format
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  utlGlobalException(_InputFormat==_OutputFormat, "input format == output format");

  typedef itk::VectorImage<double,3> ImageType;
  ImageType::Pointer inputImage = ImageType::New();
  itk::ReadVectorImage(_InputFile, inputImage);

  ImageType::Pointer outputImage = ImageType::New();
  itk::CopyImageInformation(inputImage, outputImage);
  int outDim = 6;
  if (_OutputFormat=="9D")
    outDim=9;
  else
    outDim=6;
  outputImage->SetNumberOfComponentsPerPixel(outDim);
  outputImage->Allocate();

  itk::ImageRegionConstIteratorWithIndex<ImageType> itInput(inputImage, inputImage->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> itOutput(outputImage, outputImage->GetLargestPossibleRegion());

  utl::Vector<double> inputVec;
  utl::Matrix<double> tensor;

  ImageType::PixelType inputPixel, outputPixel(outDim);
  for (itInput.GoToBegin(), itOutput.GoToBegin(); 
      !itInput.IsAtEnd(); 
      ++itInput, ++itOutput) 
    {
    inputPixel = itInput.Get();
    
    if (_InputFormat=="9D" && _OutputFormat=="6D_UPPER")
      utl::ConvertTensor9DTo6D(inputPixel, outputPixel, TENSOR_UPPER_TRIANGULAR);
    else if (_InputFormat=="6D_UPPER" && _OutputFormat=="9D")
      utl::ConvertTensor6DTo9D(inputPixel, outputPixel, TENSOR_UPPER_TRIANGULAR);

    else if (_InputFormat=="6D_EMBED" && _OutputFormat=="6D_UPPER")
      utl::ConvertTensor6DTo6D(tensor, outputPixel, TENSOR_EMBED6D, TENSOR_UPPER_TRIANGULAR);
    else if (_InputFormat=="6D_UPPER" && _OutputFormat=="6D_EMBED")
      utl::ConvertTensor6DTo6D(tensor, outputPixel, TENSOR_UPPER_TRIANGULAR, TENSOR_EMBED6D);

    else if (_InputFormat=="6D_UPPER" && _OutputFormat=="6D_LOWER")
      utl::ConvertTensor6DTo6D(inputPixel, outputPixel, TENSOR_UPPER_TRIANGULAR, TENSOR_LOWER_TRIANGULAR);
    else if (_InputFormat=="6D_LOWER" && _OutputFormat=="6D_UPPER")
      utl::ConvertTensor6DTo6D(inputPixel, outputPixel, TENSOR_LOWER_TRIANGULAR, TENSOR_UPPER_TRIANGULAR);

    else if (_InputFormat=="6D_DIAGONAL_FIRST" && _OutputFormat=="6D_LOWER")
      utl::ConvertTensor6DTo6D(inputPixel, outputPixel, TENSOR_DIAGONAL_FIRST, TENSOR_LOWER_TRIANGULAR);
    else if (_InputFormat=="6D_LOWER" && _OutputFormat=="6D_DIAGONAL_FIRST")
      utl::ConvertTensor6DTo6D(inputPixel, outputPixel, TENSOR_LOWER_TRIANGULAR, TENSOR_DIAGONAL_FIRST);

    else if (_InputFormat=="6D_DIAGONAL_FIRST" && _OutputFormat=="6D_UPPER")
      utl::ConvertTensor6DTo6D(inputPixel, outputPixel, TENSOR_DIAGONAL_FIRST, TENSOR_UPPER_TRIANGULAR);
    else if (_InputFormat=="6D_UPPER" && _OutputFormat=="6D_DIAGONAL_FIRST")
      utl::ConvertTensor6DTo6D(inputPixel, outputPixel, TENSOR_UPPER_TRIANGULAR, TENSOR_DIAGONAL_FIRST);

    else
      utlGlobalException(true, "TODO");

    itOutput.Set(outputPixel);
    }
  
  itk::SaveImage(outputImage, _OutputFile);
  
  return 0;
}
