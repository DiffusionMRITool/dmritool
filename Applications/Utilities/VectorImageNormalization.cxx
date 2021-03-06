/**
 *       @file  VectorImageNormalization.cxx
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
#include "itkNormalizeVectorImageFilter.h"
#include "utlITK.h"

#include "VectorImageNormalizationCLP.h"

#include "itkSpatiallyDenseSparseVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImageFileReader.h"
#include "itkSpatiallyDenseSparseVectorImageFileWriter.h"

/**
 * \brief  Normalize each voxel in a VectorImage or SparseVectorImage
 *
 * \author Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  utlGlobalException(_NormalizationType=="NONE", "need to set --type");

  typedef float PixelType;
  bool isSparseImage = itk::IsSparseImage(_InputFile);

  if (isSparseImage)
    {
    typedef itk::SpatiallyDenseSparseVectorImage<PixelType, 3> InputImageType;
    typedef itk::SpatiallyDenseSparseVectorImage<PixelType, 3> OutputImageType;

    InputImageType::Pointer input = InputImageType::New();
    OutputImageType::Pointer output = OutputImageType::New();
    typedef itk::SpatiallyDenseSparseVectorImageFileReader<InputImageType> ReaderType;
    typedef itk::SpatiallyDenseSparseVectorImageFileWriter<OutputImageType> WriterType;

    itk::ReadImage<InputImageType, ReaderType>(_InputFile, input);

    typedef itk::NormalizeVectorImageFilter<InputImageType, OutputImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput(input);
    if (_NormalizationType=="SUM")
      filter->SetNormalizeType(FilterType::FunctorType::SUM);
    else if (_NormalizationType=="L1NORM")
      filter->SetNormalizeType(FilterType::FunctorType::L1NORM);
    else if (_NormalizationType=="L2NORM")
      filter->SetNormalizeType(FilterType::FunctorType::L2NORM);

    // NOTE: SparseVectorImage only works for single thread
    filter->SetNumberOfThreads(1);
    filter->Update();

    output = filter->GetOutput();

    itk::SaveImage<OutputImageType, WriterType>(output, _OutputFile);
    }
  else
    {
    typedef itk::VectorImage<PixelType, 3> InputImageType;
    typedef itk::VectorImage<PixelType, 3> OutputImageType;

    InputImageType::Pointer input = InputImageType::New();
    OutputImageType::Pointer output = OutputImageType::New();

    itk::ReadImage(_InputFile, input);

    typedef itk::NormalizeVectorImageFilter<InputImageType, OutputImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput(input);
    if (_NormalizationType=="SUM")
      filter->SetNormalizeType(FilterType::FunctorType::SUM);
    else if (_NormalizationType=="L1NORM")
      filter->SetNormalizeType(FilterType::FunctorType::L1NORM);
    else if (_NormalizationType=="L2NORM")
      filter->SetNormalizeType(FilterType::FunctorType::L2NORM);

    filter->Update();

    output = filter->GetOutput();

    itk::SaveImage(output, _OutputFile);
    }

  return 0;
}
