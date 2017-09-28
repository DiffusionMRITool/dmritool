/*=========================================================================

 Program:   Sparse to Vector Image Converter

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"
#include "itkVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImageFileReader.h"
#include "itkCastVectorImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "SparseToVectorImageConverterCLP.h"

int
main(int argc, char *argv[])
{
  PARSE_ARGS;

  // Define Variables
  typedef double PixelType;
  typedef itk::SpatiallyDenseSparseVectorImage<PixelType, 3> InputImageType;
  typedef itk::VectorImage<PixelType, 3> OutputImageType;
  typedef itk::SpatiallyDenseSparseVectorImageFileReader<InputImageType> ReaderType;

  InputImageType::Pointer inputImage;
  OutputImageType::Pointer outputImage = OutputImageType::New();
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::CastVectorImageFileWriter<OutputImageType> OutputImageWriterType;
  OutputImageWriterType::Pointer writer = OutputImageWriterType::New(); 
  
  // Read input image
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
  
  // Allocate output image
  outputImage->CopyInformation(inputImage);
  outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
  outputImage->Allocate();

  // Iterator for the input image
  itk::ImageRegionConstIterator<InputImageType> inputIt(inputImage, inputImage->GetLargestPossibleRegion() );
  
  // Iterator for the output image
  itk::ImageRegionIterator<OutputImageType> outputIt(outputImage, inputImage->GetLargestPossibleRegion() );

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  
  // Transfer data
  while( !inputIt.IsAtEnd() )
    {
    outputIt.Set(inputIt.Get());
    
    ++inputIt;
    ++outputIt;
    }
  
  // Write Output
  if (_OutputFileArg.isSet())
    {
    reader->UpdateOutputInformation();
    itk::ImageIOBase::Pointer inputImageIOBase = reader->GetImageIO();
    inputImageIOBase->ReadImageInformation();
    itk::ImageIOBase::IOComponentType inputImageComponentType =
        inputImageIOBase->GetComponentType();

    try
      {
      std::cout << "Writing file: " << _OutputFile << std::endl;
      writer->SetFileName( _OutputFile );
      writer->SetInput( outputImage );
      writer->SetComponentType( inputImageComponentType );
      writer->Update();
      }
    catch ( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught!" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;  
}

