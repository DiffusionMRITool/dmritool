/**
 *       @file  DWIPreprocess.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-03-2013
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
#include "DWIPreprocessCLP.h"

/**
 * \brief  Read DWI images, gradients, b values, and output the formal normalized formalt in a configuration file.
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  // use float to save memory 
  typedef float FloatType;

  typedef itk::DWIReader<FloatType> DWIReaderType;
  DWIReaderType::Pointer reader = DWIReaderType::New();
  reader->SetConfigurationFile(_InputFile);
  if (_MaskFileArg.isSet())
    {
    typedef DWIReaderType::MaskImageType MaskImageType;
    MaskImageType::Pointer maskImage;
    itk::ReadImage<MaskImageType>(_MaskFile, maskImage);
    reader->SetMaskImage(maskImage);
    }

  reader->SetShowWarnings(_WarnArg.isSet());
  if (_DebugArg.isSet())
    reader->DebugOn();
  if (_NoNormalizeDWIArg.isSet())
    reader->SetNormalizeDWI(false);
  if (_IsVectorImageArg.isSet())
    reader->SetIsInput4DImage(false);
  if (_BThresholdArg.isSet())
    reader->GetSamplingSchemeQSpace()->SetBThresholdSingleShell(_BThreshold);
  reader->Update();

  DWIReaderType::DWIImageType::Pointer  dwiImage;
  DWIReaderType::B0ImageType::Pointer b0Image;
  dwiImage = reader->GetOutput();
  b0Image = reader->GetB0Image();

  if (_B0ImageFileArg.isSet())
    itk::SaveImage(b0Image, _B0ImageFile);

  typedef itk::DWIWriter<FloatType> DWIWriterType;
  DWIWriterType::Pointer writer = DWIWriterType::New();
  writer->SetInput(dwiImage);
  if (_OutputEachShellArg.isSet())
    writer->SetOutputEachShell(true);
  writer->SetSamplingSchemeQSpace(reader->GetSamplingSchemeQSpace());
  writer->SetConfigurationFile(_OutputFile);

  std::string _OutputFile_ext, _OutputFile_noext;
  utl::GetFileExtension(_OutputFile, _OutputFile_ext, _OutputFile_noext);
  if (_BFileArg.isSet())
    writer->SetBFile(_BFile);
  else
    writer->SetBFile(_OutputFile_noext + "_b." + _OutputFile_ext);
  if (_GradFileArg.isSet())
    writer->SetOrientationFile(_GradFile);
  else
    writer->SetOrientationFile(_OutputFile_noext + "_grad." + _OutputFile_ext);
  if (_DWIFileArg.isSet())
    writer->SetDWIFile(_DWIFile);
  else
    writer->SetDWIFile(_OutputFile_noext + "_dwi.nii.gz");

  writer->Update();

  return 0;
}
