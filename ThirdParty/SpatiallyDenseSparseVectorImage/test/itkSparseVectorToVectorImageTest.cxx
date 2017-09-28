#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"
#include "itkVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

int
itkSparseVectorToVectorImageTest(int argc, char *argv[])
{
  if (argc!=3)
    {
    std::cout << str << "Convert a sparse vector image to a vector image" << std::endl << std::flush;
    std::cout << str << argv[0] << " <inputImage> <outputImage>" << std::endl << std::flush;
    return EXIT_FAILURE;
    }

  std::string _InputFile(argv[1]);
  std::string _OutputFile(argv[2]);

  // Define Variables
  typedef float PixelType;
  typedef itk::VectorImage<PixelType, 3> VectorImageType;
  typedef itk::SpatiallyDenseSparseVectorImage<PixelType, 3> SparseVectorImageType;
  typedef itk::SpatiallyDenseSparseVectorImageFileReader<SparseVectorImageType> ReaderType;

  SparseVectorImageType::Pointer inputImage;
  VectorImageType::Pointer vectorImage = VectorImageType::New();
  ReaderType::Pointer reader = ReaderType::New();
  
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
  vectorImage->CopyInformation(inputImage);
  vectorImage->SetRegions(inputImage->GetLargestPossibleRegion());
  vectorImage->Allocate();

  // Iterator for the input image
  itk::ImageRegionConstIterator<SparseVectorImageType> inputIt(inputImage, inputImage->GetLargestPossibleRegion() );
  
  // Iterator for the output image
  itk::ImageRegionIterator<VectorImageType> outputIt(vectorImage, vectorImage->GetLargestPossibleRegion() );

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
  try
    {
    typedef itk::ImageFileWriter<VectorImageType> OutputImageWriterType;
    OutputImageWriterType::Pointer writer = OutputImageWriterType::New(); 
    std::cout << "Writing file: " << _OutputFile << std::endl;
    writer->SetFileName( _OutputFile );
    writer->SetInput( vectorImage );
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
