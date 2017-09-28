/*=========================================================================

 Program:   VTK XML Image to Image Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include <iostream>
#include <string.h>

#include "itkVectorImage.h"
#include "itkImageFileWriter.h"

#include "vtkSmartPointer.h"
#include "vtkXMLImageDataReader.h"
#include "vtkImageData.h"
#include "itkVTKImageDataToImageFilter.h"

#include "VTKXMLImageToVectorImageConverterCLP.h"

int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  
  // Define Variables
  const static unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::VectorImage<PixelType, Dimension> OutputImageType;
  typedef vtkSmartPointer<vtkXMLImageDataReader> ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  ReaderType reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  // Read data
  reader->SetFileName( _InputFile.c_str() );
  try
    {
    std::cout << "Reading file: " << _InputFile << std::endl;
    reader->Update();
    reader->GetOutput()->Register(reader);
//    reader->GetOutput()->GetPointData()->SetActiveScalars( "scalars" );
//    reader->GetOutput()->GetPointData()->SetActiveAttribute(0, vtkDataSetAttributes::SCALARS);
    }
  catch (itk::ExceptionObject & err)
    {
    std::cerr << "ExceptionObject caught!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Converter
  typedef itk::VTKImageDataToImageFilter<OutputImageType> ConverterType;
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInputData( reader->GetOutput() );
  converter->Update();

  // Write file
  try
    {
    std::cout << "Writing file: " << _OutputFile << std::endl;
    writer->SetInput( converter->GetOutput() );
    writer->SetFileName( _OutputFile );
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

