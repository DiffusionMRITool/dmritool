/*=========================================================================

 Program:   Vector Image to 4D Image Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include <iostream>
#include <string.h>
#include <stdio.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFileWriter.h"
#include "itkImageIOBase.h"
#include "VectorTo4DImageConverterCLP.h"
#include "itkVectorImageToMultiVolumeImageFilter.h"


int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  
  // Define Variables
//  typedef float PixelType;
  typedef double PixelType;
  typedef itk::VectorImage<PixelType, 3> VectorImageType;
  typedef itk::Image<PixelType, 4> MultiVolumeImageType;
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typedef itk::CastImageFileWriter<MultiVolumeImageType> WriterType;

  VectorImageType::Pointer inputImage;
  MultiVolumeImageType::Pointer outputImage = MultiVolumeImageType::New();
  WriterType::Pointer writer = WriterType::New();
  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(_InputFile);
  try
    {
    std::cout << "Reading file: " << _InputFile << std::endl;
    reader->Update();
    }
  catch (itk::ExceptionObject & err)
    {
    std::cerr << "ExceptionObject caught!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  inputImage = reader->GetOutput();
  
  typedef itk::VectorImageToMultiVolumeImageFilter<PixelType, PixelType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(inputImage);
  outputImage = convertor->GetOutput();
  convertor->Update();

  // Write Output
  reader->UpdateOutputInformation();
  itk::ImageIOBase::Pointer inputImageIOBase = reader->GetImageIO();
  inputImageIOBase->ReadImageInformation();
  itk::ImageIOBase::IOComponentType inputComponentType =
      inputImageIOBase->GetComponentType();

  try
    {
    std::cout << "Writing file: " << _OutputFile << std::endl;
    writer->SetFileName( _OutputFile );
    writer->SetInput( outputImage );
    writer->SetComponentType( inputComponentType );
    writer->Update();
    }
  catch ( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

