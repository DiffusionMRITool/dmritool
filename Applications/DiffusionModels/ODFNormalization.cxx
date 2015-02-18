/**
 *       @file  ODFNormalization.cxx
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

#include "itkMultiplyByConstantVectorImageFilter.h"
#include "ODFNormalizationCLP.h"

#include "itkNormalizeODFImageFilter.h"

/**
 * \brief  Normalize ODF represented by SH basis or uniform spherical samples
 *
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  typedef float FloatType;

  typedef itk::VectorImage<FloatType,3> InputImageType;
  typedef itk::VectorImage<FloatType,3> OutputImageType;
  
  InputImageType::Pointer input = InputImageType::New();
  itk::ReadImage(_InputFile, input);


  typedef itk::NormalizeODFImageFilter<InputImageType, OutputImageType>  FilterType;
  FilterType::Pointer filter = FilterType::New();
  
  filter->SetInput(input);

  if (_Type=="SH")
    filter->SetODFType(FilterType::FunctorType::SH);
  else if (_Type=="SAMPLE")
    filter->SetODFType(FilterType::FunctorType::SAMPLE);

  filter->Update();
  
  OutputImageType::Pointer output = filter->GetOutput();
  itk::SaveImage(output, _OutputFile);

  return 0;
}
