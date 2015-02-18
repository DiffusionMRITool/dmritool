/**
 *       @file  mexITK.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-06-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __mexITK_h
#define __mexITK_h

#include <mex.h>
#include "utlITK.h"
#include "mexutils.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{

template <class T, unsigned int VImageDimension>
inline void 
GetMXArrayFromITKImage ( const SmartPointer<Image<T,VImageDimension> >& image, mxArray*& pr  )
{
  typedef Image<T,VImageDimension> ImageType;
  utlException(IsImageEmpty(image), "the image is empty");
  
  typename ImageType::RegionType region = image->GetLargestPossibleRegion();
  typename ImageType::SizeType size = region.GetSize();
  
  utlException(VImageDimension>4, "the dimension of image is too large! VImageDimension=" << VImageDimension);
  pr = utl::Create4DImage<T>(VImageDimension>=1?size[0]:1, VImageDimension>=2?size[1]:1, VImageDimension>=3?size[2]:1, VImageDimension==4?size[3]:1);
  T * data = (T*)mxGetData(pr);

  ImageRegionConstIterator<ImageType> imageIt(image, region);
  unsigned int count = 0;
  for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, count++)
    {
    data[count] = imageIt.Get();
    }
}

template <class T, unsigned int VImageDimension>
inline void 
GetITKImageFromMXArray ( const mxArray* pr, SmartPointer<Image<T,VImageDimension> > & image  )
{
  typedef Image<T, VImageDimension> ImageType;
  if (!image)
    image = ImageType::New();
  
  mwSize dimArray = mxGetNumberOfDimensions(pr);
  const mwSize* dims = mxGetDimensions(pr);
  
  utlException(dimArray>VImageDimension, "pr has more dimensions than " << VImageDimension);
  
  typename ImageType::RegionType imageRegion;
  typename ImageType::SizeType imageSize;
  typename ImageType::IndexType imageIndex;
  
  for ( int i = 0; i < VImageDimension; i += 1 ) 
    {
    if (i<dimArray)
      imageSize[i] = static_cast<unsigned int>(dims[i]);
    else
      imageSize[i] = 1;
    imageIndex[i] = 0;
    }

  imageRegion.SetSize(imageSize);
  imageRegion.SetIndex(imageIndex);
  image->SetRegions(imageRegion);

  image->Allocate();

  ImageRegionIterator<ImageType> imageIt(image, imageRegion);
  unsigned int count = 0;
  T * data = (T*)mxGetData(pr);
  for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, count++)
    {
    imageIt.Set(data[count]);
    }
}

template <class T, unsigned int VImageDimension>
inline void 
GetMXArrayFromITKVectorImage ( const SmartPointer<VectorImage<T,VImageDimension> >& image, mxArray*& pr  )
{
  typedef VectorImage<T,VImageDimension> ImageType;
  utlException(IsImageEmpty(image), "the image is empty");
  
  typename ImageType::RegionType region = image->GetLargestPossibleRegion();
  typename ImageType::SizeType size = region.GetSize();
  int numberOfComponent = image->GetNumberOfComponentsPerPixel();
  int sizeEachVolume = 1;
  for ( int i = 0; i < VImageDimension; i += 1 ) 
    sizeEachVolume *= size[i];
  
  utlException(VImageDimension>3, "the dimension of image is too large! VImageDimension=" << VImageDimension);
  pr = utl::Create4DImage<T>(VImageDimension>=1?size[0]:1, VImageDimension>=2?size[1]:1, VImageDimension>=3?size[2]:1, numberOfComponent);
  T * data = (T*)mxGetData(pr);

  typename ImageType::PixelType pixel;
  ImageRegionConstIterator<ImageType> imageIt(image, region);
  unsigned int count = 0;
  for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, count++)
    {
    pixel = imageIt.Get();
    for ( int i = 0; i < numberOfComponent; i += 1 ) 
      {
      data[count+i*sizeEachVolume] = pixel[i];
      }
    }
}

template <class T>
inline void 
GetITKVectorImageFromMXArray ( const mxArray* pr, SmartPointer<VectorImage<T,3> >& image  )
{
  typedef VectorImage<T, 3> ImageType;
  if (!image)
    image = ImageType::New();

  mwSize VImageDimension = mxGetNumberOfDimensions(pr);
  const mwSize* dims = mxGetDimensions(pr);
  int numberOfComponent = dims[VImageDimension-1];
  
  utlException(VImageDimension>4, "pr has more dimensions than 4");
  
  typename ImageType::RegionType imageRegion;
  typename ImageType::SizeType imageSize;
  typename ImageType::IndexType imageIndex;
  
  unsigned int realDim = utl::min(3, (int)VImageDimension-1);
  int sizeEachVolume = 1;
  for ( int i = 0; i < 3; i += 1 ) 
    {
    if (i<realDim)
      {
      sizeEachVolume *= dims[i];
      imageSize[i] = static_cast<unsigned int>(dims[i]);
      }
    else
      imageSize[i] = 1;
    imageIndex[i] = 0;
    }

  imageRegion.SetSize(imageSize);
  imageRegion.SetIndex(imageIndex);
  image->SetRegions(imageRegion);
  image->SetNumberOfComponentsPerPixel(numberOfComponent);

  image->Allocate();

  typename ImageType::PixelType pixel;
  pixel.SetSize(numberOfComponent);
  ImageRegionIterator<ImageType> imageIt(image, imageRegion);
  unsigned int count = 0;
  T * data = (T*)mxGetData(pr);
  for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, count++)
    {
    for ( int i = 0; i < numberOfComponent; i += 1 ) 
      pixel[i] = data[count+i*sizeEachVolume];
    imageIt.Set(pixel);
    }
}

}

#endif 

