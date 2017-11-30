/*=========================================================================

 Program:   Image Info

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

/**
 * \todo Report data filename (as in .mhd file)
 * 
 */

#include <iostream>
#include <stdio.h>
#include <vector>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include <ImageInfoCLP.h>

#include "utlITK.h"

enum
{
  OP_NULL,
  OP_DIMENSION,
  OP_SPACING,
  OP_COMPONENT_TYPE,
  OP_PIXEL_TYPE,
  OP_NUMBER_OF_COMPONENTS,
  OP_ENDIANNESS,
  OP_SIZE_IN_BYTES,
  OP_SIZE_IN_COMPONENTS,
  OP_SIZE_IN_PIXELS,
  OP_ORIGIN,
  OP_META_DATA,
  OP_GRADIENTS,
  OP_BVALUE,
  OP_BVALUES
};

void
SetOperationWithChecking(int &operation, int value)
{
  if (operation == OP_NULL)
    {
    operation = value;
    }
  else
    {
    std::cerr << "Only one operation allowed!" << std::endl;
    exit(EXIT_FAILURE);
    }
}

void
IntRangeCheck(int x, int LowerLimit, int UpperLimit)
{
  if (x < LowerLimit || x > UpperLimit)
    {
    std::cerr << "Parameter out of range!" << std::endl;
    exit(EXIT_FAILURE);
    }
}

int
show_file(const std::string& file)
{
  std::string line;
  std::ifstream myfile (file);
  if (myfile.is_open())
    {
    while ( getline (myfile,line) )
      {
      std::cout << line << '\n';
      }
    myfile.close();
    }
  return 0;
}

int
main(int argc, char *argv[])
{
  PARSE_ARGS;

  if (itk::GetImageType(_InputFile)==IMAGE_SPARSE)
    return show_file(_InputFile);
  else if (itk::GetImageType(_InputFile)==IMAGE_VARIABLELENGTH)
    return show_file(_InputFile);

  int operation = OP_NULL;

  if (_ImageDimensionArg.isSet())
    {
    SetOperationWithChecking(operation, OP_DIMENSION);
    IntRangeCheck(_ImageDimension, 0, 3);
    }

  if (_ImageSpacingArg.isSet())
    {
    SetOperationWithChecking(operation, OP_SPACING);
    IntRangeCheck(_ImageSpacing, 0, 3);
    }
  if (_ComponentTypeArg.isSet())
    SetOperationWithChecking(operation, OP_COMPONENT_TYPE);
  if (_PixelTypeArg.isSet())
    SetOperationWithChecking(operation, OP_PIXEL_TYPE);
  if (_NumberOfComponentsPerVoxelArg.isSet())
    SetOperationWithChecking(operation, OP_NUMBER_OF_COMPONENTS);
  if (_ImageSizeInBytesArg.isSet())
    SetOperationWithChecking(operation, OP_SIZE_IN_BYTES);
  if (_ImageSizeInComponentsArg.isSet())
    SetOperationWithChecking(operation, OP_SIZE_IN_COMPONENTS);
  if (_ImageSizeInPixelsArg.isSet())
    SetOperationWithChecking(operation, OP_SIZE_IN_PIXELS);
  if (_OriginArg.isSet())
    SetOperationWithChecking(operation, OP_ORIGIN);
  if (_MetaDataArg.isSet())
    SetOperationWithChecking(operation, OP_META_DATA);
  if (_GradientsArg.isSet())
    SetOperationWithChecking(operation, OP_GRADIENTS);
  if (_DiffusionWeightingArg.isSet())
    SetOperationWithChecking(operation, OP_BVALUE);
  if (_DiffusionWeightingsArg.isSet())
    SetOperationWithChecking(operation, OP_BVALUES);
  
  // Define Variables
  typedef float PixelType;
  static const unsigned int ImageDimension = 4;
  typedef itk::VectorImage<PixelType, ImageDimension> MultiVolumeVectorImageType;
  typedef itk::ImageFileReader<MultiVolumeVectorImageType> ReaderType;

  MultiVolumeVectorImageType::Pointer inputImage;
  ReaderType::Pointer reader = ReaderType::New();

  // Read Input File Info
  reader->SetFileName(_InputFile);
  reader->UpdateOutputInformation();
  inputImage = reader->GetOutput();
  itk::ImageIOBase::Pointer inputImageIOBase = reader->GetImageIO();
  inputImageIOBase->ReadImageInformation();

  itk::ImageIOBase::IOPixelType inputPixelType = inputImageIOBase->GetPixelType();
  itk::ImageIOBase::IOComponentType inputComponentType =
      inputImageIOBase->GetComponentType();
  itk::ImageIOBase::SizeType inputImageSizeInBytes =
      inputImageIOBase->GetImageSizeInBytes();
  itk::ImageIOBase::SizeType inputImageSizeInComponents =
      inputImageIOBase->GetImageSizeInComponents();
  itk::ImageIOBase::SizeType inputImageSizeInPixels =
    inputImageIOBase->GetImageSizeInPixels();

  itk::ImageIOBase::ByteOrder inputByteOrder = inputImageIOBase->GetByteOrder();
  //	bool useCompression = inputImageIOBase->GetUseCompression();

  std::string inputPixelTypeAsString(
    inputImageIOBase->GetPixelTypeAsString(inputPixelType));
  std::string inputComponentTypeAsString(
    inputImageIOBase->GetComponentTypeAsString(inputComponentType));
  std::string inputByteOrderAsString(inputImageIOBase->GetByteOrderAsString(inputByteOrder));

  unsigned int inputNumberOfComponentsPerPixel =
      inputImage->GetNumberOfComponentsPerPixel();

  MultiVolumeVectorImageType::RegionType inputRegion =
      inputImage->GetLargestPossibleRegion();
  MultiVolumeVectorImageType::SizeType inputSize = inputRegion.GetSize();
  MultiVolumeVectorImageType::SpacingType inputSpacing = inputImage->GetSpacing();
  MultiVolumeVectorImageType::PointType inputOrigin = inputImage->GetOrigin();
  MultiVolumeVectorImageType::DirectionType inputDirection = inputImage->GetDirection();

  // Meta Data
  itk::MetaDataDictionary MetaDictionary = inputImage->GetMetaDataDictionary();
  std::vector < std::string > MetaKeys = MetaDictionary.GetKeys();
  std::vector<std::string>::const_iterator itKey;
  std::string metaString;
  
  double x, y, z;
  double b0 = 0;

  std::cout.precision(10);
  switch (operation)
  {
    case OP_NULL:
      std::cout << "Dimension =" << std::flush;
      for (unsigned int k=0; k<ImageDimension; k++)
        {
        std::cout << " " << inputSize[k] << std::flush;
        }
      std::cout << std::endl;
      
      std::cout << "Spacing =" << std::flush;
      for (unsigned int k=0; k<ImageDimension; k++)
        {
        std::cout << " " << inputSpacing[k] << std::flush;
        }
      std::cout << std::endl;
      
      std::cout << "Origin =" << std::flush;
      for (unsigned int k=0; k<ImageDimension; k++)
        {
        std::cout << " " << inputOrigin[k] << std::flush;
        }
      std::cout << std::endl;
    
      std::cout << "ComponentType = " << inputComponentTypeAsString
        << std::endl;
      std::cout << "PixelType = " << inputPixelTypeAsString << std::endl;
      std::cout << "NumberOfComponents = " << inputNumberOfComponentsPerPixel
        << std::endl;
      std::cout << "ImageSizeInBytes = " << inputImageSizeInBytes << std::endl;
      std::cout << "ImageSizeInComponents = " << inputImageSizeInComponents
        << std::endl;
      std::cout << "ImageSizeInPixels = " << inputImageSizeInPixels << std::endl;
//      std::cout << "UseCompression = " << useCompression << std::endl;
//      std::cout << "ByteOrder = " << inputByteOrderAsString << std::endl;
      std::cout << "DirectionType = [" << std::flush;
      for (unsigned int r=0; r<ImageDimension; r++)
        {
        for (unsigned int c=0; c<ImageDimension; c++)
          {
          std::cout << inputDirection(r,c) << std::flush;
          if (c < ImageDimension - 1)
            {
            std::cout << " " << std::flush;
            }
          }
        if (r < ImageDimension - 1)
          {
          std::cout << "; " << std::flush;
          }
        }
      std::cout << "]" << std::endl;
    
      break;
    case OP_DIMENSION:
      std::cout << inputSize[_ImageDimension] << std::flush;
      break;
    
    case OP_SPACING:
//      std::cout << vcl_abs(inputSpacing[_ImageSpacing]);
      std::cout << inputSpacing[_ImageSpacing] << std::flush;
      break;
    
    case OP_COMPONENT_TYPE:
      std::cout << inputComponentTypeAsString;
      break;
    
    case OP_PIXEL_TYPE:
      std::cout << inputPixelTypeAsString;
      break;
    
    case OP_NUMBER_OF_COMPONENTS:
      std::cout << inputNumberOfComponentsPerPixel;
      break;
    
    case OP_SIZE_IN_BYTES:
      std::cout << inputImageSizeInBytes;
      break;
    
    case OP_SIZE_IN_COMPONENTS:
      std::cout << inputImageSizeInComponents;
      break;
    
    case OP_SIZE_IN_PIXELS:
      std::cout << inputImageSizeInPixels;
      break;
    
    case OP_ORIGIN:
      std::cout << inputOrigin;
      break;
    
    case OP_META_DATA:
      for (itKey = MetaKeys.begin(); itKey != MetaKeys.end(); ++itKey)
      {
        itk::ExposeMetaData < std::string
        > (MetaDictionary, *itKey, metaString);
        std::cout << *itKey << " ---> " << metaString << std::endl;
      }
      break;
    
    case OP_GRADIENTS:
      for (itKey = MetaKeys.begin(); itKey != MetaKeys.end(); ++itKey)
      {
        itk::ExposeMetaData < std::string> (MetaDictionary, *itKey, metaString);
        if (itKey->find("DWMRI_gradient") != std::string::npos)
        {
          sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
          //          std::cout << x << " " << y << " " << " " << z << std::endl;
          printf("%lf %lf %lf\n",x,y,z); // To ensure proper display formatting
        }
      }
      break;
    
    case OP_BVALUE:
      for (itKey = MetaKeys.begin(); itKey != MetaKeys.end(); ++itKey)
      {
        itk::ExposeMetaData < std::string
        > (MetaDictionary, *itKey, metaString);
        if (itKey->find("DWMRI_b-value") != std::string::npos)
        {
          b0 = atof(metaString.c_str());
          //          std::cout << b0 << std::endl;
          printf("%lf\n",b0); // To ensure proper display formatting
        }
      }
      break;
    
    case OP_BVALUES:
      for (itKey = MetaKeys.begin(); itKey != MetaKeys.end(); ++itKey)
      {
        itk::ExposeMetaData < std::string
        > (MetaDictionary, *itKey, metaString);
        if (itKey->find("DWMRI_b-value") != std::string::npos)
        {
          b0 = atof(metaString.c_str());
          //          printf("%lf\n",b0); // To ensure proper display formatting
        }
      }
      for (itKey = MetaKeys.begin(); itKey != MetaKeys.end(); ++itKey)
      {
        itk::ExposeMetaData < std::string
        > (MetaDictionary, *itKey, metaString);
        if (itKey->find("DWMRI_gradient") != std::string::npos)
        {
          sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
          //          printf("%lf %lf %lf\n",x,y,z); // To ensure proper display formatting
          printf("%lf\n", (x*x+y*y+z*z)*b0); // To ensure proper display formatting
        }
      }
      break;
    
    default:
      break;
  }
  
  return EXIT_SUCCESS;
}
