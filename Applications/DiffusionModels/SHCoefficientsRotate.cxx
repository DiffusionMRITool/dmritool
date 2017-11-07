/**
 *       @file  SHCoefficientsRotate.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "12-30-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#include "utl.h"
#include "itkRotateSHCoefficientsImageFilter.h"
#include "utlRotationMatrixFromVectors.h"
#include "SHCoefficientsRotateCLP.h"



int
main(int argc, char *argv[])
{
  // GenerateCLP
  PARSE_ARGS;
  
  utlException(_InputFile=="" || _OutputFile=="", "no input or no output!");

  // Time Probe
  itk::TimeProbe clock;
  
  // Input and Output Image Types
  typedef itk::VectorImage<double, 3> ImageType;
  ImageType::Pointer shImage, shRotatedImage;

  itk::ReadImage<ImageType>(_InputFile, shImage);
  // int dimension = shImage->GetNumberOfComponentsPerPixel();

  typedef itk::RotateSHCoefficientsImageFilter<ImageType, ImageType>  RotateFilterType;  
  RotateFilterType::Pointer rotateFilter = RotateFilterType::New();

  utlException(_VectorFromTo.size()!=4 && _VectorFromTo.size()!=6, "wrong size! _VectorFromTo.size()="<<_VectorFromTo.size());
  RotateFilterType::VectorType vFrom(3), vTo(3);
  if (_VectorFromTo.size()==6)
    {
    for ( int i = 0; i < 3; i += 1 ) 
      vFrom[i] = _VectorFromTo[i];
    for ( int i = 0; i < 3; i += 1 ) 
      vTo[i] = _VectorFromTo[i+3];
    vFrom /= vFrom.GetTwoNorm();
    vTo /= vFrom.GetTwoNorm();
    }
  else
    {
    vFrom[0]=vTo[0]=1.0;
    vFrom[1]=_VectorFromTo[0]*M_PI/180.0, vFrom[2]=_VectorFromTo[1]*M_PI/180.0;
    vTo[1]=_VectorFromTo[2]*M_PI/180.0, vTo[2]=_VectorFromTo[3]*M_PI/180.0;
    utl::spherical2Cartesian(vFrom[0], vFrom[1], vFrom[2]);
    utl::spherical2Cartesian(vTo[0], vTo[1], vTo[2]);
    }

  utl::PrintUtlVector(vFrom, "vFrom");
  utl::PrintUtlVector(vTo, "vTo");
  
  RotateFilterType::MatrixType rotationMatrix(3,3);
  utl::RotationMatrixFromUnitNormVectors<RotateFilterType::VectorType, RotateFilterType::MatrixType >(vFrom, vTo, rotationMatrix);
  utl::PrintUtlMatrix(rotationMatrix, "rotationMatrix");

  rotateFilter->SetInput(shImage);
  rotateFilter->SetRotationMatrix(rotationMatrix);
  // rotateFilter->InPlaceOn();
  if (_NumberOfThreads>0)
    rotateFilter->SetNumberOfThreads(_NumberOfThreads);
  clock.Start();
  rotateFilter->Update();
  clock.Stop();
  std::cout << clock.GetMean() << "s elapsed" << std::endl;

  shRotatedImage = rotateFilter->GetOutput();

  itk::SaveImage<ImageType>(shRotatedImage, _OutputFile);
  
}
