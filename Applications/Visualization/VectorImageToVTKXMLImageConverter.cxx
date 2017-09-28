/*=========================================================================

 Program:   Image to VTK XML Image Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include <iostream>
#include <string.h>

#include "itkVectorImage.h"
#include "itkImageFileReader.h"

#include "vtkSmartPointer.h"
#include "vtkXMLImageDataWriter.h"
#include "itkImageToVTKImageDataFilter.h"

#include "utl.h"
#include "VectorImageToVTKXMLImageConverterCLP.h"

int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  
  // Define Variables
  const static unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::VectorImage<PixelType, Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef vtkSmartPointer<vtkXMLImageDataWriter> WriterType;
  WriterType writer = WriterType::New();

  // Read data
  ImageType::Pointer image = ImageType::New();
  itk::ReadVectorImage(_InputFile, image);

  // Converter
  typedef itk::ImageToVTKImageDataFilter<ImageType> ConverterType;
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInput( image );
  converter->Update();

  // Write file
  try
    {
    std::cout << "Writing file: " << _OutputFile << std::endl;
    writer->SetInputData( converter->GetOutput() );
    writer->SetFileName( _OutputFile.c_str() );
    writer->Write();
    }
  catch ( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

