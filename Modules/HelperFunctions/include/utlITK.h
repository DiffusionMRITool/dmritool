/**
 *       @file  utlITK.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-01-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */



#ifndef __utlITK_h
#define __utlITK_h

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkVariableLengthVector.h>
#include <itkNumericTraits.h>

#include "utlCore.h"
#include "utlITKConceptChecking.h"

#include "utlITKMacro.h"

namespace itk 
{

/** @addtogroup utlHelperFunctions
@{ */

template <class T>
void
PrintVariableLengthVector(const VariableLengthVector<T>vec, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  char tmp[1024];
  int NSize=vec.GetSize();
  std::vector<double> st = utl::GetContainerStats(vec.GetDataPointer(), vec.GetDataPointer()+vec.Size());
  sprintf(tmp, "%-8s(%p):  size = %lu,  stat = { %g, %g [%g], %g } : ", str==""?"vector":str.c_str(), &vec, (long unsigned int)NSize, st[0], st[2], st[3], st[1] );
  std::string strr(tmp);
  if (NSize>0)
    {
    os << strr << " = [ ";
    for ( int i = 0; i < NSize-1; i += 1 ) 
      os << vec[i] << separate;
    os << vec[NSize-1] << " ];\n";
    }
  else
    os << strr << " is empty vector" << std::endl;

}

template <class T>
VariableLengthVector<T>
VnlVectorToVariableLengthVector ( const vnl_vector<T>& vec )
{
  VariableLengthVector<T> v(vec.size());
  for ( int i = 0; i < vec.size(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
vnl_vector<T> 
VariableLengthVectorToVnlVector ( const VariableLengthVector<T>& vec )
{
  vnl_vector<T> v(vec.GetSize());
  for ( int i = 0; i < vec.GetSize(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

/** Read information of the image, not the data in image  */
template <class ImageType>
bool 
ReadImageInformation (const std::string filename, SmartPointer<ImageType>& image) 
{
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(filename);
  try 
    {
    reader->UpdateOutputInformation(); 
    } 
  catch (itk::ExceptionObject & err) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return false;
    }
  image = reader->GetOutput();
  return true;
}

/** Read Image  */
template <class ImageType>
bool 
ReadImage (const std::string filename, SmartPointer<ImageType>& image, const std::string printInfo="Reading Image:") 
{
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(filename);
  try 
    {
    if (utl::IsLogNormal())
      std::cout << printInfo << " " << filename << std::endl;
    reader->Update(); 
    } 
  catch (itk::ExceptionObject & err) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return false;
    }
  image = reader->GetOutput();
  return true;
}

template <class ImageType, class ReaderType>
bool 
ReadImage (const std::string filename, SmartPointer<ImageType>& image, const std::string printInfo="Reading Image:") 
{
  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(filename);
  try 
    {
    if (utl::IsLogNormal())
      std::cout << printInfo << " " << filename << std::endl;
    reader->Update(); 
    } 
  catch (itk::ExceptionObject & err) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return false;
    }
  image = reader->GetOutput();
  return true;
}

template <class ImageType>
bool 
SaveImage (const SmartPointer<ImageType>& image, const std::string filename, const std::string printInfo="Writing Image:")
{
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  
  writer->SetFileName(filename);
  writer->SetInput(image);
  try 
    {
    if (utl::IsLogNormal())
      std::cout << printInfo << " " << filename << std::endl;
    writer->Update(); 
    } 
  catch (itk::ExceptionObject & err) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return false;
    }
  return true;
}

template <class ImageType, class WriterType>
bool 
SaveImage (const SmartPointer<ImageType>& image, const std::string filename, const std::string printInfo="Writing Image:")
{
  typename WriterType::Pointer writer = WriterType::New();
  
  writer->SetFileName(filename);
  writer->SetInput(image);
  try 
    {
    if (utl::IsLogNormal())
      std::cout << printInfo << " " << filename << std::endl;
    writer->Update(); 
    } 
  catch (itk::ExceptionObject & err) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return false;
    }
  return true;
}

inline bool
IsVectorImage(const std::string filename)
{
  typedef itk::VectorImage<float, 4> MultiVolumeVectorImageType;
  typedef itk::ImageFileReader<MultiVolumeVectorImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  reader->UpdateOutputInformation();
  MultiVolumeVectorImageType::Pointer image = reader->GetOutput();

  unsigned int numberOfComponentsPerPixel = image->GetNumberOfComponentsPerPixel();
  return numberOfComponentsPerPixel > 1;
}

inline bool
Is3DImage(const std::string filename)
{
  if (IsVectorImage(filename))
    return false;

  typedef itk::Image<float, 4> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  reader->UpdateOutputInformation();
  ImageType::Pointer image = reader->GetOutput();

  typename ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  return size[3]==1;
}

template <class ImageType>
int 
GetVectorImageVectorSize(const SmartPointer<ImageType>& image, const int axis=-1)
{
  std::string name = image->GetNameOfClass();
  if (name=="VectorImage" || name=="SpatiallyDenseSparseVectorImage")
    {
    if (axis<0 || axis==ImageType::ImageDimension)
      return image->GetNumberOfComponentsPerPixel(); 
    else
      {
      typename ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
      return size[axis]; 
      }
    }
  else if (name=="Image")
    {
    typename ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    if (axis==ImageType::ImageDimension || (axis<0 && ImageType::ImageDimension<=3))
      return 1;
    else if (axis<0 || axis==ImageType::ImageDimension-1)
      return size[ImageType::ImageDimension-1];
    else
      return size[axis];
    }
  else
    {
    utlGlobalException(true, "not supported type: " + name);
    return -1;
    }
}

template <class ImageType>
void
SetVectorImageVectorSize(const SmartPointer<ImageType>& image, const int vecsize)
{
  std::string name = image->GetNameOfClass();
  if (name=="VectorImage" || name=="SpatiallyDenseSparseVectorImage")
    {
    image->SetNumberOfComponentsPerPixel(vecsize); 
    }
  else if (name=="Image")
    {
    typename ImageType::RegionType region = image->GetLargestPossibleRegion();
    typename ImageType::SizeType size = region.GetSize();
    size[ImageType::ImageDimension-1]=vecsize;
    region.SetSize(size);
    image->SetRegions(region);
    }
  else
    utlGlobalException(true, "not supported type: " + name);
}

template <class ImageType>
std::vector<int>
GetVectorImage3DVolumeSize(const SmartPointer<ImageType>& image)
{
  std::string name = image->GetNameOfClass();
  std::vector<int> size(3, 1);
  typename ImageType::SizeType imagesize = image->GetLargestPossibleRegion().GetSize();
  if (name=="VectorImage" || name=="SpatiallyDenseSparseVectorImage")
    {
    for ( int i = 0; i < ImageType::ImageDimension; ++i ) 
      {
      utlException(i>2 && imagesize[i]!=1, "the image has more than 3D. Image size = " << imagesize);
      size[i] = imagesize[i];
      }
    }
  else if (name=="Image")
    {
    for ( int i = 0; i < utl::min<int>(ImageType::ImageDimension, 3); ++i ) 
      size[i] = imagesize[i];
    for ( int i = 4; i < ImageType::ImageDimension; ++i ) 
      utlException(imagesize[i]!=1, "the image has more than 3D. Image size = " << imagesize);
    }
  else
    utlGlobalException(true, "not supported type: " + name);
  return size;
}

template <class ImageType>
void
SetVectorImage3DVolumeSize(SmartPointer<ImageType>& image, const std::vector<int>& size)
{
  std::string name = image->GetNameOfClass();
  utlException(size.size()!=3, "need to have size 3");
  typename ImageType::RegionType region = image->GetLargestPossibleRegion();
  typename ImageType::SizeType imagesize = region.GetSize();
  for ( int i = 0; i < utl::min<int>(ImageType::ImageDimension, size.size()); ++i ) 
    imagesize[i] = size[i];
  region.SetSize(imagesize);
  image->SetRegions(region);
}

/** Get 4d size from 3D vector image or 4D scalar image */
template <class ImageType>
std::vector<int>
GetVectorImageFullSize(const SmartPointer<ImageType>& image)
{
  std::string name = image->GetNameOfClass();
  std::vector<int> size;
  typename ImageType::SizeType imagesize = image->GetLargestPossibleRegion().GetSize();
  const int dim = ImageType::ImageDimension;
  if (name=="VectorImage" || name=="SpatiallyDenseSparseVectorImage")
    {
    size.resize(dim+1);
    for ( int i = 0; i < dim; ++i ) 
      size[i] = imagesize[i];
    size[dim] = image->GetNumberOfComponentsPerPixel();
    }
  else if (name=="Image")
    {
    size.resize(dim);
    for ( int i = 0; i < dim; ++i ) 
      size[i] = imagesize[i];
    }
  else
    utlGlobalException(true, "not supported type: " + name);
  return size;
}

/** Set 4d size to 3D vector image or 4D scalar image */
template <class ImageType>
void
SetVectorImageFullSize(SmartPointer<ImageType>& image, const std::vector<int>& size)
{
  std::string name = image->GetNameOfClass();
  typename ImageType::RegionType region = image->GetLargestPossibleRegion();
  typename ImageType::SizeType imagesize = region.GetSize();
  const int dim = ImageType::ImageDimension;
  if (name=="VectorImage" || name=="SpatiallyDenseSparseVectorImage")
    {
    utlSAException(size.size()!=dim+1)(size.size())(dim).msg("wrong size");
    for ( int i = 0; i < dim; ++i ) 
      imagesize[i] = size[i];
    image->SetNumberOfComponentsPerPixel(size[dim]);
    }
  else if (name=="Image")
    {
    utlSAException(size.size()!=dim)(size.size())(dim).msg("wrong size");
    for ( int i = 0; i < dim; ++i ) 
      imagesize[i] = size[i];
    }
  else
    utlGlobalException(true, "not supported type: " + name);
  region.SetSize(imagesize);
  image->SetRegions(region);
}

inline bool
IsSparseImage(const std::string filename)
{
  std::string fileNoExt, ext;
  utl::GetFileExtension(filename, ext, fileNoExt);
  return ext=="spr";
}

template <class ImageType >
bool
IsImageEmpty(const SmartPointer<ImageType>& image)
{
  if (!image)
    return true;

  typename ImageType::RegionType region = image->GetLargestPossibleRegion();
  typename ImageType::SizeType size = region.GetSize();

  for ( int i = 0; i < ImageType::ImageDimension; i += 1 ) 
    {
    if (size[i]>0)
      return false;
    }
  return true;
}

/** Generate an image with given size and VectorLength.  */
template <class ImageType>
typename ImageType::Pointer
GenerateImage(const typename ImageType::SizeType& size, const int vectorLength=1)
{
  typename ImageType::Pointer image = ImageType::New();
  typename ImageType::RegionType region;
  region.SetSize(size);
  image->SetRegions(region);
  image->SetNumberOfComponentsPerPixel(vectorLength);
  image->Allocate();
  return image;
}

template <class ImageType>
typename ImageType::Pointer
GenerateImageFromSingleVoxel(const typename ImageType::PixelType& pixel)
{
  typename ImageType::SizeType sizeVoxel;
  typename ImageType::IndexType indexVoxel;
  typedef typename ImageType::PixelType PixelType;
  for ( int i = 0; i < sizeVoxel.Size(); ++i ) 
    {
    sizeVoxel[i]=1.0;
    indexVoxel[i]=0.0;
    }
  typename ImageType::Pointer image = itk::GenerateImage<ImageType>(sizeVoxel, itk::NumericTraits<PixelType>::GetLength(pixel));
  image->SetPixel(indexVoxel, pixel);
  return image;
}

template <class Image1Type, class Image2Type>
void
ImageToImage ( const SmartPointer<Image1Type>& image1, SmartPointer<Image2Type>& image2)
{
  utlGlobalException(!image2, "image2 cannot be null");
  image2->CopyInformation(image1);
  image2->SetRegions(image1->GetLargestPossibleRegion());
  image2->Allocate();
  itk::ImageRegionConstIterator<Image1Type> it1(image1, image1->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<Image2Type> it2(image2, image2->GetLargestPossibleRegion() );

  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    it2.Set(it1.Get());
    ++it2;
    ++it1;
    }
}

template <class Image1Type, class Image2Type>
SmartPointer<Image2Type>
ImageToImage ( const SmartPointer<Image1Type>& image1 )
{
  typename Image2Type::Pointer image2 = Image2Type::New();
  ImageToImage<Image1Type, Image2Type>(image1, image2);
  return image2;
}

template <unsigned dimIn, unsigned dimOut>
void CopyImageRegion(const ImageRegion<dimIn>& regionIn, ImageRegion<dimOut>& regionOut, const int numberOfComponens=-1)
{
  typename ImageRegion<dimIn>::SizeType sizeIn = regionIn.GetSize();
  typename ImageRegion<dimIn>::IndexType indexIn = regionIn.GetIndex();

  typename ImageRegion<dimOut>::SizeType sizeOut = regionOut.GetSize();
  typename ImageRegion<dimOut>::IndexType indexOut = regionOut.GetIndex();

  int dimension = utl::min((int)dimIn, (int)dimOut);

  for ( int i = 0; i < dimension; ++i ) 
    {
    sizeOut[i] = sizeIn[i];
    indexOut[i] = indexIn[i];
    }

  if ((int)dimIn < (int)dimOut)
    {
    for ( int i = dimension; i < dimOut; i += 1 ) 
      {
      if (dimension+1==dimOut )
        {
        utlSAException(numberOfComponens==-1)(numberOfComponens).msg("need to set numberOfComponens");
        sizeOut[i] = numberOfComponens;
        }
      else
        sizeOut[i] = 1;
      indexOut[i] = 0;
      }
    }

  regionOut.SetSize(sizeOut);
  regionOut.SetIndex(indexOut);
}


/** Copy Image Information from image with different image type */
template <class ImageWithInfoType, class ImageType>
void 
CopyImageInformation ( const SmartPointer<ImageWithInfoType>& imageFrom, SmartPointer<ImageType>& imageTo )
{
  // utlException(ImageWithInfoType::ImageDimension < ImageType::ImageDimension, "the dimensions are not compatible");
  int dimension = utl::min((int)ImageWithInfoType::ImageDimension, (int)ImageType::ImageDimension);
  // typename ImageWithInfoType::PixelType    inputImagePixelValue;
  // typename ImageWithInfoType::IndexType    inputImagePixelIndex;
  typename ImageWithInfoType::SpacingType  inputImageSpacing = imageFrom->GetSpacing();
  typename ImageWithInfoType::RegionType inRegion = imageFrom->GetLargestPossibleRegion();
  typename ImageWithInfoType::SizeType inputImageSize = inRegion.GetSize();
  typename ImageWithInfoType::IndexType inputImageIndex = inRegion.GetIndex();
  typename ImageWithInfoType::PointType inOrigin = imageFrom->GetOrigin();
  typename ImageWithInfoType::DirectionType inDirection = imageFrom->GetDirection();

  typename ImageType::RegionType imageRegion;
  typename ImageType::SizeType imageSize;
  typename ImageType::PixelType imagePixelValue;
  typename ImageType::IndexType imagePixelIndex;
  typename ImageType::SpacingType imageSpacing;
  typename ImageType::PointType imageOrigin;
  typename ImageType::DirectionType imageDirection;

  for ( int i = 0; i < dimension; i += 1 ) 
    {
    imageSpacing[i] = inputImageSpacing[i];
    imageSize[i] = inputImageSize[i];
    imageOrigin[i] = inOrigin[i];
    imagePixelIndex[i] = inputImageIndex[i];
    for ( int j = 0; j < dimension; j += 1 ) 
      imageDirection(i,j) = inDirection(i,j);
    }

  if ((int)ImageWithInfoType::ImageDimension < (int)ImageType::ImageDimension)
    {
    for ( int i = dimension; i < ImageType::ImageDimension; i += 1 ) 
      {
      if (dimension+1==ImageType::ImageDimension && std::string(imageFrom->GetNameOfClass())=="VectorImage" && std::string(imageTo->GetNameOfClass())=="Image" )
        imageSize[i] = imageFrom->GetNumberOfComponentsPerPixel();
      else
        imageSize[i] = 1;
      imageSpacing[i] = 1;
      imageOrigin[i] = 0;
      imagePixelIndex[i] = 0;
      imageDirection(i,i) = 1.0;
      }
    }

  imageRegion.SetSize(imageSize);
  imageRegion.SetIndex(imagePixelIndex);

  imageTo->SetSpacing(imageSpacing);
  imageTo->SetOrigin(imageOrigin);
  imageTo->SetRegions(imageRegion);
  imageTo->SetDirection(imageDirection);

  imageTo->SetNumberOfComponentsPerPixel(GetVectorImageVectorSize(imageFrom));
}

/** print itk::VectorImage  */
template <class TPixelType, unsigned int VImageDimension >
void
PrintVectorImage(const SmartPointer<VectorImage<TPixelType, VImageDimension> >& image, const std::string mse="", std::ostream& os=std::cout, bool isPrintHeader=false)
{
  typedef VectorImage<TPixelType, VImageDimension> VectorImageType;
  if (isPrintHeader)
    image->Print(os<<mse);
  ImageRegionIteratorWithIndex<VectorImageType> it(image, image->GetLargestPossibleRegion());
  typename VectorImageType::IndexType index;
  typename VectorImageType::PixelType pixel;
  int numberOfComponent = image->GetNumberOfComponentsPerPixel();
  it.GoToBegin();
  while(!it.IsAtEnd())
    {
    pixel = it.Get();
    std::vector<double> st = utl::GetContainerStats(pixel.GetDataPointer(), pixel.GetDataPointer()+pixel.GetSize());
    if (std::fabs(st[0])>1e-8 || std::fabs(st[1])>1e-8)
      {
      index = it.GetIndex();
      os << (mse==""?"VectorImage":mse) <<"(";
      for ( int i = 0; i < VectorImageType::ImageDimension; i += 1 ) 
        {
        if (i==VectorImageType::ImageDimension-1)
          os << index[i] << ")";
        else
          os << index[i] << ",";
        }
      char tmp[1024];
      sprintf(tmp, " : size = %lu,  stat = { %g, %g [%g], %g }  : ", (long unsigned)pixel.GetSize(), st[0], st[2], st[3], st[1] );
      os << std::string(tmp);
      os << "[ ";
      for ( int i = 0; i < numberOfComponent; i += 1 ) 
        {
        if (i==numberOfComponent-1)
          os << pixel[i] << " ];" << std::endl;
        else
          os << pixel[i] << ", ";
        }
      }

    ++it;
    }
}

/** print itk::Image with dimension no more than 3 */
template <class TPixelType, unsigned int VImageDimension >
void
PrintImage3D(const SmartPointer<Image<TPixelType,VImageDimension> >& image, const std::string mse="", std::ostream& os=std::cout, bool isPrintHeader=false)
{
  utlGlobalException(VImageDimension>3, "dimension should be no more than 3!");
  typedef Image<TPixelType, VImageDimension> ImageType;
  if (isPrintHeader)
    image->Print(os<<mse);

  ImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion());
  typename ImageType::IndexType index;
  TPixelType pixel;
  it.GoToBegin();
  while(!it.IsAtEnd())
    {
    pixel = it.Get();
    if (std::fabs(pixel)>1e-8)
      {
      index = it.GetIndex();
      os << (mse==""?"Image":mse) <<"(";
      for ( int i = 0; i < ImageType::ImageDimension; i += 1 ) 
        {
        if (i==ImageType::ImageDimension-1)
          os << index[i] << ")";
        else
          os << index[i] << ",";
        }
      os << " : [ ";
      os << pixel << " ];" << std::endl;
      }
    ++it;
    }
}

/** print 4D itk::Image   */
template <class TPixelType, unsigned int VImageDimension>
void
PrintImage4D(const SmartPointer<Image<TPixelType,VImageDimension> >& image, const std::string mse="", std::ostream& os=std::cout, bool isPrintHeader=false)
{
  utlGlobalException(VImageDimension!=4, "wrong image dimension");
  typedef Image<TPixelType, VImageDimension> ImageType;
  typename ImageType::RegionType region = image->GetLargestPossibleRegion();
  typename ImageType::SizeType size = region.GetSize();
  
  if (isPrintHeader)
    image->Print(os<<mse);

  typedef Image<TPixelType, 3> Image3DType;
  typename Image3DType::Pointer image3D = Image3DType::New();
  CopyImageInformation<ImageType, Image3DType>(image, image3D);

  ImageRegionIteratorWithIndex<Image3DType> it(image3D, image3D->GetLargestPossibleRegion());
  typename Image3DType::IndexType index3D;
  typename ImageType::IndexType index;
  it.GoToBegin();
  while(!it.IsAtEnd())
    {
    index3D = it.GetIndex();
    std::vector<TPixelType> pixel;
    for ( int i = 0; i < 3; i += 1 ) 
      index[i] = index3D[i];
    for ( int i = 0; i < size[3]; i += 1 ) 
      {
      index[3] = i;
      pixel.push_back( image->GetPixel(index) );
      }

    std::vector<double> st = utl::GetContainerStats(pixel.begin(), pixel.end());
    if (std::fabs(st[0])>1e-8 || std::fabs(st[1])>1e-8)
      {
      os << (mse==""?"Image":mse) <<"(";
      for ( int i = 0; i < 3; i += 1 ) 
        {
        if (i==2)
          os << index3D[i] << ")";
        else
          os << index3D[i] << ",";
        }
      char tmp[1024];
      sprintf(tmp, " : size = %lu,  stat = { %g, %g [%g], %g }  : ", (long unsigned)pixel.size(), st[0], st[2], st[3], st[1] );
      os << std::string(tmp);
      os << "[ ";
      for ( int i = 0; i < size[3]; i += 1 ) 
        {
        if (i==size[3]-1)
          os << pixel[i] << " ];" << std::endl;
        else
          os << pixel[i] << ", ";
        }
      }
    
    ++it;
    }
}

/** print itk::Image  */
template <class TPixelType, unsigned int VImageDimension >
void
PrintImage(const SmartPointer<Image<TPixelType,VImageDimension> > image, const std::string mse="", std::ostream& os=std::cout, bool isPrintHeader=false)
{
  if (VImageDimension==4)
    PrintImage4D<TPixelType, VImageDimension>(image, mse, os, isPrintHeader);
  else if (VImageDimension<=3)
    PrintImage3D<TPixelType, VImageDimension>(image, mse, os, isPrintHeader);
  else
    utlGlobalException(true, "image dimension is larger than 4");
}

/** Check Image Information */
template <class Image1Type, class Image2Type>
bool
VerifyImageSize ( const SmartPointer<Image1Type>& image1, const SmartPointer<Image2Type>& image2, const bool isMinimalDimension=false )
{
  int dimension;
  if (isMinimalDimension)
    dimension = utl::min((int)Image1Type::ImageDimension, (int)Image2Type::ImageDimension, 3);
  else
    {
    if ((int)Image1Type::ImageDimension!=(int)Image2Type::ImageDimension)
      return false;
    dimension = Image1Type::ImageDimension;
    }
  typename Image1Type::RegionType region1 = image1->GetLargestPossibleRegion();
  typename Image1Type::SizeType size1 = region1.GetSize();
  
  typename Image2Type::RegionType region2 = image2->GetLargestPossibleRegion();
  typename Image2Type::SizeType size2 = region2.GetSize();
  
  for ( int i = 0; i < dimension; i += 1 ) 
    {
    if ( size1[i] != size2[i] )
      return false;
    }
  return true;
}

/** Check Image Information */
template <class Image1Type>
bool
VerifyImageSize ( const SmartPointer<Image1Type>& image1, const std::string file2, const bool isMinimalDimension=false )
{
  if (isMinimalDimension)
    {
    static const unsigned int ImageDimension = 4;
    typedef double PixelType;
    typedef VectorImage<PixelType, ImageDimension> Image2Type;

    typedef itk::ImageFileReader<Image2Type> Reader2Type;
    typename Reader2Type::Pointer reader2 = Reader2Type::New();
    reader2->SetFileName(file2);
    reader2->UpdateOutputInformation();
    typename Image2Type::Pointer image2 = reader2->GetOutput();

    return VerifyImageSize<Image1Type, Image2Type> ( image1, image2, isMinimalDimension);
    }
  else
    {
    typedef Image1Type  Image2Type;
    typedef itk::ImageFileReader<Image2Type> Reader2Type;
    typename Reader2Type::Pointer reader2 = Reader2Type::New();
    reader2->SetFileName(file2);
    reader2->UpdateOutputInformation();
    typename Image2Type::Pointer image2 = reader2->GetOutput();

    return VerifyImageSize<Image1Type, Image2Type> ( image1, image2, isMinimalDimension);
    }
}

/** Check Image Information */
inline bool
VerifyImageSize ( const std::string file1, const std::string file2, const bool isMinimalDimension=false )
{
  static const unsigned int ImageDimension = 4;
  typedef double PixelType;
  typedef VectorImage<PixelType, ImageDimension> Image1Type;
  typedef VectorImage<PixelType, ImageDimension> Image2Type;

  typedef itk::ImageFileReader<Image1Type> Reader1Type;
  Reader1Type::Pointer reader1 = Reader1Type::New();
  reader1->SetFileName(file1);
  reader1->UpdateOutputInformation();
  Image1Type::Pointer image1 = reader1->GetOutput();
  
  typedef itk::ImageFileReader<Image2Type> Reader2Type;
  Reader2Type::Pointer reader2 = Reader2Type::New();
  reader2->SetFileName(file2);
  reader2->UpdateOutputInformation();
  Image2Type::Pointer image2 = reader2->GetOutput();

  return VerifyImageSize<Image1Type, Image2Type> ( image1, image2, isMinimalDimension);
}

/** Check Image Information */
template <class Image1Type, class Image2Type>
bool
VerifyImageInformation ( const SmartPointer<Image1Type>& image1, const SmartPointer<Image2Type>& image2, const bool isMinimalDimension=false )
{
  int dimension;
  if (isMinimalDimension)
    {
    dimension = utl::min((int)Image1Type::ImageDimension, (int)Image2Type::ImageDimension, 3);
    }
  else
    {
    if ((int)Image1Type::ImageDimension!=(int)Image2Type::ImageDimension)
      return false;
    dimension = Image1Type::ImageDimension;
    }

  typename Image1Type::RegionType region1 = image1->GetLargestPossibleRegion();
  typename Image1Type::SizeType size1 = region1.GetSize();
  typename Image1Type::SpacingType spacing1 = image1->GetSpacing();
  typename Image1Type::PointType origin1 = image1->GetOrigin();
  typename Image1Type::DirectionType direction1 = image1->GetDirection();
  
  typename Image2Type::RegionType region2 = image2->GetLargestPossibleRegion();
  typename Image2Type::SizeType size2 = region2.GetSize();
  typename Image2Type::SpacingType spacing2 = image2->GetSpacing();
  typename Image2Type::PointType origin2 = image2->GetOrigin();
  typename Image2Type::DirectionType direction2 = image2->GetDirection();

  // image1->Print(std::cout<<"image1");
  // image2->Print(std::cout<<"image2");
  for ( int i = 0; i < dimension; i += 1 ) 
    {
    if ( (std::fabs(spacing1[i]-spacing2[i])>1e-6) || (std::fabs(size1[i]-size2[i])>1e-6) || (std::fabs(origin1[i]-origin2[i])>1e-6) )
      return false;
    for ( int j = 0; j < dimension; j += 1 ) 
      {
      if ( std::fabs(direction1(i,j) - direction2(i,j))>1e-6 )
        return false;
      }
    }
  return true;
}

/** Check Image Information */
template <class Image1Type>
bool
VerifyImageInformation ( const SmartPointer<Image1Type>& image1, const std::string file2, const bool isMinimalDimension=false )
{
  if (isMinimalDimension)
    {
    static const unsigned int ImageDimension = 4;
    typedef double PixelType;
    typedef VectorImage<PixelType, ImageDimension> Image2Type;

    typedef itk::ImageFileReader<Image2Type> Reader2Type;
    typename Reader2Type::Pointer reader2 = Reader2Type::New();
    reader2->SetFileName(file2);
    reader2->UpdateOutputInformation();
    typename Image2Type::Pointer image2 = reader2->GetOutput();

    return VerifyImageInformation<Image1Type, Image2Type> ( image1, image2, isMinimalDimension);
    }
  else
    {
    typedef Image1Type  Image2Type;
    typedef itk::ImageFileReader<Image2Type> Reader2Type;
    typename Reader2Type::Pointer reader2 = Reader2Type::New();
    reader2->SetFileName(file2);
    reader2->UpdateOutputInformation();
    typename Image2Type::Pointer image2 = reader2->GetOutput();

    return VerifyImageInformation<Image1Type, Image2Type> ( image1, image2, isMinimalDimension);
    }
}

/** Check Image Information */
inline bool
VerifyImageInformation ( const std::string file1, const std::string file2, const bool isMinimalDimension=false )
{
  static const unsigned int ImageDimension = 4;
  typedef double PixelType;
  typedef VectorImage<PixelType, ImageDimension> Image1Type;
  typedef VectorImage<PixelType, ImageDimension> Image2Type;

  typedef itk::ImageFileReader<Image1Type> Reader1Type;
  Reader1Type::Pointer reader1 = Reader1Type::New();
  reader1->SetFileName(file1);
  reader1->UpdateOutputInformation();
  Image1Type::Pointer image1 = reader1->GetOutput();
  
  typedef itk::ImageFileReader<Image2Type> Reader2Type;
  Reader2Type::Pointer reader2 = Reader2Type::New();
  reader2->SetFileName(file2);
  reader2->UpdateOutputInformation();
  Image2Type::Pointer image2 = reader2->GetOutput();

  return VerifyImageInformation<Image1Type, Image2Type> ( image1, image2, isMinimalDimension);
}


template <class PointsContainer, class VnlValueType>
void 
PointsContainerToVnlMatrix ( const PointsContainer& points, vnl_matrix<VnlValueType>& matrix )
{
  matrix.clear();
  typedef typename PointsContainer::ConstIterator PointsIterator;
  typedef typename PointsContainer::Element PointType;

  const unsigned int pointDimension = PointType::PointDimension;
  unsigned int numberOfPoints = points.Size();
  if (numberOfPoints==0)
    return;

  matrix.set_size(numberOfPoints, pointDimension);
  PointsIterator iterator = points.Begin();
  PointsIterator end = points.End();
  
  unsigned int count=0;
  while( iterator != end )
    {
    PointType orientation = iterator.Value();
    for (unsigned int k=0; k<pointDimension; k++)
      matrix(count,k) = orientation[k];
    iterator++;
    count++;
    }
}

template <class VnlValueType, class PointsContainer>
void 
VnlMatrixToPointsContainer ( const vnl_matrix<VnlValueType>& matrix, PointsContainer& points )
{
  points.Initialize();
  typedef typename PointsContainer::Element PointType;

  const unsigned int pointDimension = PointType::PointDimension;
  utlGlobalException(pointDimension!=matrix.columns(), "wrong size of matrix or point dimension. pointDimension="<< pointDimension << ", matrix.columns()="<<matrix.columns());

  for ( int i = 0; i < matrix.rows(); i += 1 ) 
    {
    PointType point;
    for ( int j = 0; j < matrix.columns(); j += 1 ) 
      {
      point[j] = matrix(i,j);
      }
    points.InsertElement(i, point);
    }
}

    /** @} */

}

#endif 

