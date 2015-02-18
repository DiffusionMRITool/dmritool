/**
 *       @file  ITKImageWrite.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-14-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#include "mex.h" 
#include "utlMEX.h"

template <class ImageType>
void 
SetITKImageInformation ( itk::SmartPointer<ImageType>& image, const mxArray* originArray, const mxArray* spacingArray )
{
  typename ImageType::PointType origin;
  for ( int i = 0; i < ImageType::ImageDimension; i += 1 ) 
    origin[i]=0;
  if (originArray)
    {
    mwSize dimOrgin = mxGetNumberOfElements(originArray);
    double * data = mxGetPr(originArray);
    for ( int i = 0; i < utl::min((int)dimOrgin, (int)ImageType::ImageDimension); i += 1 ) 
      origin[i] = data[i];
    }
  typename ImageType::SpacingType spacing;
  for ( int i = 0; i < ImageType::ImageDimension; i += 1 ) 
    spacing[i]=1;
  if (spacingArray)
    {
    mwSize dimSpacing = mxGetNumberOfElements(spacingArray);
    double * data = mxGetPr(spacingArray);
    for ( int i = 0; i < utl::min((int)dimSpacing, (int)ImageType::ImageDimension); i += 1 ) 
      spacing[i] = data[i];
    }
  image->SetOrigin(origin);
  image->SetSpacing(spacing);
}

template <typename T>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs, const bool isVectorImage) 
{
  std::string filename = utl::GetString(prhs[0]);

  mxArray* originArray=NULL;
  mxArray* spacingArray=NULL;
  if (nrhs==3)
    {
    originArray = utl::GetArrayStruct(prhs[2], "origin" );
    spacingArray = utl::GetArrayStruct(prhs[2], "spacing" );
    }
  
  mwSize dimArray = mxGetNumberOfDimensions(prhs[1]);
  const mwSize* dims = mxGetDimensions(prhs[1]);

  if (isVectorImage)
    {
    utlGlobalException(dimArray!=4, "The image has to have 4 dimension");
    typedef itk::VectorImage<T, 3> ImageType;
    typename ImageType::Pointer image = ImageType::New();
    itk::GetITKVectorImageFromMXArray(prhs[1], image);
    SetITKImageInformation(image, originArray, spacingArray);
    itk::SaveImage(image, filename);
    }
  else
    {
    typedef itk::Image<T, 4> ImageType;
    typename ImageType::Pointer image = ImageType::New();
    itk::GetITKImageFromMXArray(prhs[1], image);
    SetITKImageInformation(image, originArray, spacingArray);
    itk::SaveImage(image, filename);
    }

}


void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs!=2 && nrhs!=3, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=0, "the program should no output");

  utlGlobalException(mxGetClassID(prhs[0])!=mxCHAR_CLASS, "the first input has to be a string");
  
  bool isVectorImage = utl::GetScalarStructDef<bool>(nrhs==3?prhs[2]:NULL,"vectorImage",false);
  
  if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS) 
    {
    std::cout << "mxSINGLE_CLASS, float" << std::endl << std::flush;
    callFunction<float>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxLOGICAL_CLASS)
    {
    std::cout << "mxLOGICAL_CLASS, mxLogical" << std::endl << std::flush;
    callFunction<mxLogical>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxUINT8_CLASS)
    {
    std::cout << "mxUINT8_CLASS, unsigned char" << std::endl << std::flush;
    callFunction<unsigned char>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxUINT16_CLASS)
    {
    std::cout << "mxUINT16_CLASS, unsigned short int" << std::endl << std::flush;
    callFunction<unsigned short int>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxUINT32_CLASS)
    {
    std::cout << "mxUINT32_CLASS, unsigned int" << std::endl << std::flush;
    callFunction<unsigned int>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxUINT64_CLASS)
    {
    std::cout << "mxUINT64_CLASS, uint64_T" << std::endl << std::flush;
    callFunction<uint64_T>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxINT8_CLASS)
    {
    std::cout << "mxINT8_CLASS, signed char" << std::endl << std::flush;
    callFunction<signed char>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxINT16_CLASS)
    {
    std::cout << "mxINT16_CLASS, short int" << std::endl << std::flush;
    callFunction<short int>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxINT32_CLASS)
    {
    std::cout << "mxINT32_CLASS, int" << std::endl << std::flush;
    callFunction<int>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxINT64_CLASS)
    {
    std::cout << "mxINT64_CLASS, int64_T" << std::endl << std::flush;
    callFunction<int64_T>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else if (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS)
    {
    std::cout << "mxDOUBLE_CLASS, double" << std::endl << std::flush;
    callFunction<double>(plhs,prhs,nlhs,nrhs, isVectorImage);
    }
  else
    utlGlobalException(true, "wrong type!!");
}

