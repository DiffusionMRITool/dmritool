/**
 *       @file  itkVectorImageRegionIteratorWithIndexTest.cxx
 *      @brief  
 *     Created  "07-29-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */
// #include "utlCom"
#include "utlITK.h"
#include "utlNDArray.h"
#include "utl.h"
#include "itkVectorImageRegionIteratorWithIndex.h"
#include "itkVectorImageRegionIterator.h"

int 
main (int argc, char const* argv[])
{
  typedef itk::Image<double,4> NDImageType;
  typedef itk::Image<double,3> ScalarImageType;
  typedef itk::VectorImage<double,3> VectorImageType;

  NDImageType::SizeType size4d;
  VectorImageType::SizeType size3d;
  int vecSize=5;
  size4d[0]=size3d[0]=2, size4d[1]=size3d[1]=3, size4d[2]=size3d[2]=4, size4d[3]=vecSize;
  int totalNum = size4d[0]*size4d[1]*size4d[2]*size4d[3];
  int total3d = size3d[0]*size3d[1]*size3d[2];

  utl::NDArray<double, 4> array(2,3,4,5);
  unsigned ss3d[3];
  ss3d[0]=2, ss3d[1]=3, ss3d[2]=4;
  utl::NDArray<double, 3> array3d(ss3d);
  for ( int i = 0; i < totalNum; ++i ) 
    array[i] = i;
  for ( int i = 0; i < total3d; ++i ) 
    array3d[i] = i;

  // row-major, the last index changes fast. 
  utl::PrintUtlNDArray(array, "array");
  array.PrintWithIndex(std::cout<< "array");

    {

    // colume-major, the first index changes fast. 
    NDImageType::Pointer image4d = itk::GenerateImage<NDImageType>(size4d);
    NDImageType::PixelContainer::Pointer pixelContainer = NDImageType::PixelContainer::New();
    pixelContainer->SetImportPointer(array.GetData(), array.Size());
    image4d->SetPixelContainer(pixelContainer);
    image4d->Print(std::cout<< "image4d=\n");
    itk::PrintImage4D(image4d, "image4d");
    
    NDImageType::IndexType index;
    itk::VariableLengthVector<double> vec, vec1;

      {
      // t-axis
      std::cout << "\nVector in t-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<NDImageType> it(image4d, image4d->GetLargestPossibleRegion());
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        std::cout << "index = " << index << std::endl << std::flush;
        it.GetVector(vec);
        itk::PrintVariableLengthVector(vec, "vec");
        it.SetVector(2*vec);
        it.GetVector(vec1);
        itk::PrintVariableLengthVector(vec1, "vec1");
        it.SetVector(vec);
        }
      }

      {
      // x-axis
      std::cout << "\nVector in x-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<NDImageType> it(image4d, image4d->GetLargestPossibleRegion(),0);
      std::cout << "it.GetRegion() = " << it.GetRegion() << std::endl << std::flush;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        for ( int i = 0; i < size4d[3]; ++i ) 
          {
          index = it.GetIndex();
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          it.SetVector(2*vec, i);
          it.GetVector(vec1, i);
          itk::PrintVariableLengthVector(vec1, "vec1");
          it.SetVector(vec, i);
          }
        }
      }
      
      {
      // x-axis
      std::cout << "\nVector in x-axis (half region in y-axis)" << std::endl << std::flush;
      auto region = image4d->GetLargestPossibleRegion();
      auto regionIndex = region.GetIndex();
      auto regionSize = region.GetSize();
      regionIndex[1]=1;
      regionSize[1]=1;
      region.SetIndex(regionIndex);
      region.SetSize(regionSize);
      std::cout << "region = " << region << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<NDImageType> it(image4d, region, 0);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        for ( int i = 0; i < size4d[3]; ++i ) 
          {
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          it.SetVector(2*vec, i);
          it.GetVector(vec1, i);
          itk::PrintVariableLengthVector(vec1, "vec1");
          it.SetVector(vec, i);
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          }
        }
      }
      
      {
      // y-axis
      std::cout << "\nVector in y-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<NDImageType> it(image4d, image4d->GetLargestPossibleRegion(),1);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        for ( int i = 0; i < size4d[3]; ++i ) 
          {
          index = it.GetIndex();
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          }
        }
      }

    }
    

    {

    // colume-major, the first index changes fast. 
    VectorImageType::Pointer image4d = itk::GenerateImage<VectorImageType>(size3d, vecSize);
    VectorImageType::PixelContainer::Pointer pixelContainer = VectorImageType::PixelContainer::New();
    pixelContainer->SetImportPointer(array.GetData(), array.Size());
    image4d->SetPixelContainer(pixelContainer);
    image4d->Print(std::cout<< "image4d=\n");
    itk::PrintVectorImage(image4d, "image4d");
    
    VectorImageType::IndexType index;
    itk::VariableLengthVector<double> vec, vec1;

      {
      // t-axis
      std::cout << "\nVector in t-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<VectorImageType> it(image4d, image4d->GetLargestPossibleRegion());
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        std::cout << "index = " << index << std::endl << std::flush;
        it.GetVector(vec);
        itk::PrintVariableLengthVector(vec, "vec");
        it.SetVector(2*vec);
        it.GetVector(vec1);
        itk::PrintVariableLengthVector(vec1, "vec1");
        it.SetVector(vec);
        }
      }
      {
      // x-axis
      std::cout << "\nVector in t-axis (half region in y-axis)" << std::endl << std::flush;
      auto region = image4d->GetLargestPossibleRegion();
      auto regionIndex = region.GetIndex();
      auto regionSize = region.GetSize();
      regionIndex[1]=1;
      regionSize[1]=1;
      region.SetIndex(regionIndex);
      region.SetSize(regionSize);
      std::cout << "region = " << region << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<VectorImageType> it(image4d, region);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        for ( int i = 0; i < 1; ++i ) 
          {
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          it.SetVector(2*vec, i);
          it.GetVector(vec1, i);
          itk::PrintVariableLengthVector(vec1, "vec1");
          it.SetVector(vec, i);
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          }
        }
      }
      
      {
      // x-axis
      std::cout << "\nVector in x-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<VectorImageType> it(image4d, image4d->GetLargestPossibleRegion(), 0);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        for ( int i = 0; i < size4d[3]; ++i ) 
          {
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          it.SetVector(2*vec, i);
          it.GetVector(vec1, i);
          itk::PrintVariableLengthVector(vec1, "vec1");
          it.SetVector(vec, i);
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          }
        }
      }
      
      {
      // x-axis
      std::cout << "\nVector in x-axis (half region in y-axis)" << std::endl << std::flush;
      auto region = image4d->GetLargestPossibleRegion();
      auto regionIndex = region.GetIndex();
      auto regionSize = region.GetSize();
      regionIndex[1]=1;
      regionSize[1]=1;
      region.SetIndex(regionIndex);
      region.SetSize(regionSize);
      std::cout << "region = " << region << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<VectorImageType> it(image4d, region, 0);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        for ( int i = 0; i < size4d[3]; ++i ) 
          {
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          it.SetVector(2*vec, i);
          it.GetVector(vec1, i);
          itk::PrintVariableLengthVector(vec1, "vec1");
          it.SetVector(vec, i);
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          }
        }
      }

      {
      // y-axis
      std::cout << "\nVector in y-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<VectorImageType> it(image4d, image4d->GetLargestPossibleRegion(), 1);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        for ( int i = 0; i < size4d[3]; ++i ) 
          {
          std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
          it.GetVector(vec, i);
          itk::PrintVariableLengthVector(vec, "vec");
          it.SetVector(2*vec,i);
          it.GetVector(vec1,i);
          itk::PrintVariableLengthVector(vec1, "vec1");
          it.SetVector(vec,i);
          itk::PrintVariableLengthVector(vec, "vec");
          }
        }
      }
    }


    {

    std::cout << "\n\n\n***********************\n" << std::endl << std::flush;
    std::cout << "\n\n\ntest 3d image \n\n" << std::endl << std::flush;

    // colume-major, the first index changes fast. 
    vecSize=1;
    size4d[3]=1;

    NDImageType::Pointer image4d = itk::GenerateImage<NDImageType>(size4d);
    NDImageType::PixelContainer::Pointer pixelContainer = NDImageType::PixelContainer::New();
    pixelContainer->SetImportPointer(array.GetData(), total3d);
    image4d->SetPixelContainer(pixelContainer);
    image4d->Print(std::cout<< "image4d=\n");
    itk::PrintImage4D(image4d, "image4d");
    
    NDImageType::IndexType index;
    itk::VariableLengthVector<double> vec, vec1;

      {
      // t-axis
      std::cout << "\nVector in t-axis" << std::endl << std::flush;
      itk::VectorImageRegionIteratorWithIndex<NDImageType> it(image4d, image4d->GetLargestPossibleRegion(),-1);
      for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
        {
        index = it.GetIndex();
        std::cout << "index = " << index << std::endl << std::flush;
        std::cout << "val = " << it.Get() << std::endl << std::flush;
        it.GetVector(vec);
        itk::PrintVariableLengthVector(vec, "vec");
        it.SetVector(2*vec);
        it.GetVector(vec1);
        itk::PrintVariableLengthVector(vec1, "vec1");
        it.SetVector(vec);
        }
      }


    }
  
  return 0;
}
