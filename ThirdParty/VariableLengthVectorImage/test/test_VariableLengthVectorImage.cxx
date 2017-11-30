/**
 *       @file  test_VariableLengthVectorImage.cxx
 *      @brief  
 *     Created  "11-09-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "itkImage.h"
#include "itkVariableLengthVector.h"
#include "utlITK.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVariableLengthVectorImageFileReader.h"
#include "itkVariableLengthVectorImageFileWriter.h"

int 
main (int argc, char const* argv[])
{
  const unsigned int Dimension = 3;
  const unsigned int VectorLength = 6;
  typedef double PixelType;
  typedef itk::Image< itk::VariableLengthVector< PixelType >, Dimension > VariableLengthVectorImageType;
  typedef itk::VariableLengthVector< PixelType > InternalPixelType;
  typedef itk::VariableLengthVectorImageFileReader< VariableLengthVectorImageType > VariableLengthVectorImageReader;
  typedef itk::VariableLengthVectorImageFileWriter< VariableLengthVectorImageType > VariableLengthVectorImageWriter;

  VariableLengthVectorImageType::Pointer image = VariableLengthVectorImageType::New();
  VariableLengthVectorImageType::IndexType start, indTmp;
  InternalPixelType f( VectorLength ), f1, f2(3);
  f2[0]=1, f2[1]=0, f2[2]=-1;
  VariableLengthVectorImageType::SizeType size;
  for( unsigned int i=0; i<VectorLength; i++ ) { f[i] = i+1; }
  start.Fill(0);
  size.Fill(2);
  VariableLengthVectorImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer( f );
  indTmp[0]=0, indTmp[1]=1, indTmp[2]=0;
  image->SetPixel( indTmp, f2);


  image->Print(std::cout<<"image = \n");

  itk::ImageRegionIteratorWithIndex<VariableLengthVectorImageType> it(image, image->GetLargestPossibleRegion() );
  for ( it.GoToBegin();  !it.IsAtEnd(); ++it ) 
    {
    f1 = it.Get();
    std::cout << "index = " << it.GetIndex() << std::endl << std::flush;
    std::cout << "pixel = " << f1 << std::endl << std::flush;
    }

  itk::SaveImage<VariableLengthVectorImageType, VariableLengthVectorImageWriter>(image, "test.vlv");

    {
    VariableLengthVectorImageType::Pointer image2 = VariableLengthVectorImageType::New();
    itk::ReadImage<VariableLengthVectorImageType, VariableLengthVectorImageReader>("test.vlv", image2);
    image2->Print(std::cout<<"image2 = \n");
    itk::ImageRegionIteratorWithIndex<VariableLengthVectorImageType> it(image2, image2->GetLargestPossibleRegion() );
    for ( it.GoToBegin();  !it.IsAtEnd(); ++it ) 
      {
      f1 = it.Get();
      std::cout << "index = " << it.GetIndex() << std::endl << std::flush;
      std::cout << "pixel = " << f1 << std::endl << std::flush;
      }
    }
  
  return 0;
}
