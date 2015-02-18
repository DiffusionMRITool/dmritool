/**
 *       @file  ImageMultiplication.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */
#include "utlITK.h"

#include "ImageMultiplicationCLP.h"
#include "itkMultiplyByConstantVectorImageFilter.h"


/**
 * \brief  mutiple a const to a image 
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  typedef float FloatType;

  // NOTE: use 4D to work for VectorImage<FloatType,3> and Image<FloatType,4>
  typedef itk::VectorImage<FloatType,4> InputImageType;
  typedef itk::VectorImage<FloatType,4> OutputImageType;
  typedef InputImageType::PixelType  ContainerType;

  InputImageType::Pointer input = InputImageType::New();
  itk::ReadImage(_InputFile, input);
  int pixelDimension = input->GetNumberOfComponentsPerPixel();

  typedef itk::MultiplyByConstantVectorImageFilter<InputImageType, ContainerType, OutputImageType> MultiplyFilterType;
  MultiplyFilterType::Pointer filter = MultiplyFilterType::New();


  ContainerType vec;
  vec.SetSize(pixelDimension);
  for ( int i = 0; i < pixelDimension; i += 1 ) 
    vec[i] = _Scale;

  filter->SetInput(input);
  filter->SetConstantVector(vec);

  filter->Update();
  
  OutputImageType::Pointer output = filter->GetOutput();

  itk::SaveImage(output, _OutputFile);
  
  return 0;
}
