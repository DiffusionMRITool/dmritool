/**
 *       @file  itkDWIReaderTest.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-18-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkDWIReader.h"
#include "itkDWIWriter.h"

#include "utl.h"
#include "itkDWIReaderTestCLP.h"



/**
 * \brief  test DWIReader
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  typedef itk::DWIReader<double> DWIReaderType;
  DWIReaderType::Pointer reader = DWIReaderType::New();
  reader->SetConfigurationFile(_InputFile);
  if (_MaskFileArg.isSet())
    {
    typedef DWIReaderType::MaskImageType MaskImageType;
    MaskImageType::Pointer maskImage;
    itk::ReadImage<MaskImageType>(_MaskFile, maskImage);
    reader->SetMaskImage(maskImage);
    }
  if (_DebugArg.isSet())
    reader->DebugOn();
  if (_NoNormalizeDWIArg.isSet())
    reader->SetNormalizeDWI(false);
  if (_IsVectorImageArg.isSet())
    reader->SetIsInput4DImage(false);
  if (_B0FileArg.isSet())
    {
    typedef DWIReaderType::B0ImageType B0ImageType;
    B0ImageType::Pointer b0Image = B0ImageType::New();
    itk::ReadImage(_B0File, b0Image);
    reader->SetB0Image(b0Image);
    }
  reader->Update();
  reader->Print(std::cout<<"reader");

  DWIReaderType::DWIImageType::Pointer  dwiImage;
  DWIReaderType::B0ImageType::Pointer b0Image;
  dwiImage = reader->GetOutput();
  b0Image = reader->GetB0Image();

  // dwiImage->Print(std::cout<<"dwiImage = ");
  // b0Image->Print(std::cout<<"b0Image");

  if (_PrintArg.isSet())
    {
    if (b0Image)
      itk::PrintImage(b0Image, "b0Image");

    itk::PrintVectorImage(dwiImage, "dwiImage");
    }
  // itk::PrintVectorImage<double,3>(dwiImage, "dwiImage");
  // if (b0Image)
  //   itk::PrintImage<double,3>(b0Image, "b0Image");

  // DWIReaderType::DWIImageType::IndexType  dwiIndex;
  // dwiIndex[0]=26, dwiIndex[1]=12, dwiIndex[2]=1;
  // std::cout << "dwiImage(26,12,1) = " << dwiImage->GetPixel(dwiIndex) << std::endl << std::flush;
  // std::cout << "b0Image(26,12,1) = " << b0Image->GetPixel(dwiIndex) << std::endl << std::flush;

  if (_OutputFileArg.isSet())
    {
    typedef itk::DWIWriter<double> DWIWriterType;
    DWIWriterType::Pointer writer = DWIWriterType::New();
    writer->SetInput(dwiImage);
    if (_OutputEachShellArg.isSet())
      writer->SetOutputEachShell(true);
    writer->SetSamplingSchemeQSpace(reader->GetSamplingSchemeQSpace());
    writer->SetConfigurationFile(_OutputFile);
    writer->SetBFile("b_test.txt");
    writer->SetOrientationFile("grad_test.txt");
    writer->SetDWIFile("dwi_test.nii.gz");
    writer->Update();
    }

  return 0;
}
