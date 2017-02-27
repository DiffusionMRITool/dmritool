/**
 *       @file  ITKImageRead.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-11-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#include "mex.h" 
#include "utlMEX.h"

#include "itkImageFileReader.h"

template <typename T, unsigned int dim, bool isVectorImage>
   inline void callFunction(mxArray* plhs[], const mxArray* prhs[],
         const int nlhs,const int nrhs) 
{

  std::string filename = utl::GetString(prhs[0]);

  double * out_origin=0;
  double * out_spacing=0;
  if (nlhs>1)
    {
    plhs[1] = utl::CreateMatrix<T>(1, dim);
    plhs[2] = utl::CreateMatrix<T>(1, dim);
    out_origin = mxGetPr(plhs[1]);
    out_spacing = mxGetPr(plhs[2]);
    }
  
  if (isVectorImage)
    {
    typedef itk::VectorImage<T, dim> ImageType;
    typename ImageType::Pointer imageVec = ImageType::New();
    itk::ReadImage(filename, imageVec);
    itk::GetMXArrayFromITKVectorImage(imageVec, plhs[0]);

    if (nlhs>1)
      {
      typename ImageType::SpacingType  spacing = imageVec->GetSpacing();
      typename ImageType::PointType origin = imageVec->GetOrigin();
      for ( int i = 0; i < dim; i += 1 ) 
        {
        out_origin[i] = origin[i];
        out_spacing[i] = spacing[i];
        }
      }
    }
  else
    {
    typedef itk::Image<T, dim> ImageType;
    typename ImageType::Pointer image4D = ImageType::New();
    itk::ReadImage(filename, image4D);
    itk::GetMXArrayFromITKImage(image4D, plhs[0]);

    if (nlhs>1)
      {
      typename ImageType::SpacingType  spacing = image4D->GetSpacing();
      typename ImageType::PointType origin = image4D->GetOrigin();
      for ( int i = 0; i < dim; i += 1 ) 
        {
        out_origin[i] = origin[i];
        out_spacing[i] = spacing[i];
        }
      }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  utlGlobalException(nrhs!=1, "Bad number of inputs arguments");
  utlGlobalException(nlhs!=1 && nlhs!=3, "Bad number of outputs arguments");
  
  utlGlobalException(mxGetClassID(prhs[0])!=mxCHAR_CLASS, "the first input has to be a string");
  std::string filename = utl::GetString(prhs[0]);

  const unsigned int Dimension = 4;

  typedef itk::VectorImage<float, Dimension> MultiVolumeVectorImageType;
  typedef itk::ImageFileReader<MultiVolumeVectorImageType> ReaderType;
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  reader->UpdateOutputInformation();
  MultiVolumeVectorImageType::Pointer image = reader->GetOutput();

  unsigned int numberOfComponentsPerPixel = image->GetNumberOfComponentsPerPixel();
  bool isVectorImage = numberOfComponentsPerPixel>1;

  // image->Print(std::cout<<"image=");
  // itk::PrintVectorImage(image,"image");
  MultiVolumeVectorImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  int dim = Dimension;
  // if size=[1,1,1,1], dim=3. dim is 3 or 4.
  for ( int i = Dimension-1; i > 2 ; i-- ) 
    {
    if (size[i]==1)
      dim--;
    else
      break;
    }

  if (dim==4)
    isVectorImage? callFunction<double,4,true>(plhs,prhs,nlhs,nrhs) : callFunction<double,4,false>(plhs,prhs,nlhs,nrhs);
  else if (dim==3)
    isVectorImage? callFunction<double,3,true>(plhs,prhs,nlhs,nrhs) : callFunction<double,3,false>(plhs,prhs,nlhs,nrhs);
  // else if (dim==2)
  //   isVectorImage? callFunction<double,2,true>(plhs,prhs,nlhs,nrhs) : callFunction<double,2,false>(plhs,prhs,nlhs,nrhs);
  // else if (dim==1)
  //   isVectorImage? callFunction<double,1,true>(plhs,prhs,nlhs,nrhs) : callFunction<double,1,false>(plhs,prhs,nlhs,nrhs);
} 

