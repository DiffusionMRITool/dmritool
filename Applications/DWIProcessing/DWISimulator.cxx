/*=========================================================================

 Program:   DWI Simulator

 Copyright (c) Pew-Thian Yap, Jian Cheng. All rights reserved.
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

#include <iostream>
#include <fstream>
#include "DWISimulatorCLP.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkDWIGenerator.h"
#include "itkDWIWriter.h"
#include "utl.h"

int
main(int argc, char *argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  // Input and Output Image Types
  typedef itk::VectorImage<float, 3> OutputImageType;
  OutputImageType::Pointer outputImage;

  typedef itk::Image<float, 3> ScalarImageType;
  ScalarImageType::Pointer outputScalarImage;

  // DWI Generator 
  typedef itk::DWIGenerator<OutputImageType, ScalarImageType> GeneratorType;
  GeneratorType::Pointer dwiGenerator = GeneratorType::New();

  dwiGenerator->SetDebug(_DebugArg.isSet());
 
  // Generate Data
  dwiGenerator->SetFileName( _ParameterFile );
  utlGlobalException(_OutputDWIFile=="" && _OutputODFFile=="" && _OutputEAPFile=="" 
    && _OutputMSDFile=="" && _OutputRTOFile=="", "need to set output dwi/odf/eap/peak/rto/msd");
  dwiGenerator->SetIsOutputDWI(_OutputDWIFile!="");
  dwiGenerator->SetIsOutputEAP(_OutputEAPFile!="");
  dwiGenerator->SetIsOutputODF(_OutputODFFile!="");
  dwiGenerator->SetIsOutputRTO(_OutputRTOFile!="");
  dwiGenerator->SetIsOutputMSD(_OutputMSDFile!="");
  if (_MaxNumberOfPeaks>0)
    {
    dwiGenerator->SetMaxNumberOfPeaks(_MaxNumberOfPeaks);
    dwiGenerator->SetPeakType(itk::PeakContainerHelper::GetPeakType(_PeakType));
    }

  if ( _SNRArg.isSet() )
    dwiGenerator->SetSNR( _SNR );
  if (_NoiseSigmaArg.isSet())
    dwiGenerator->SetNoiseSigma( _NoiseSigma );
  
  utlGlobalException(_SNRArg.isSet() && _NoiseSigmaArg.isSet(), "Only one of these is allowed: --snr or --noisesigma");
  
  dwiGenerator->GetSamplingSchemeQSpace()->SetTau(_Tau);
  dwiGenerator->GetSamplingSchemeRSpace()->SetTau(_Tau);
  
  utlGlobalException(!_QSpaceOrientationsArg.isSet() && !_RSpaceOrientationsArg.isSet(), "no gradient files in Q space and R space");
  
  GeneratorType::STDVectorPointer bvec(new GeneratorType::STDVectorType());
  GeneratorType::STDVectorPointer rvec(new GeneratorType::STDVectorType());

  if (_QSpaceOrientationsArg.isSet())
    {
    utlGlobalException(!_BValuesArg.isSet() && !_BFileArg.isSet(), "need to set --bvalues or --bfile");
    utlGlobalException(_BValuesArg.isSet() && _BFileArg.isSet(), "Only one of these is allowed: --bvalues or --bfile");

    if (_BValuesArg.isSet())
      *bvec = _BValues;
    if (_BFileArg.isSet())
      utl::ReadVector(_BFile,*bvec);

    GeneratorType::MatrixPointer qMatrix =  utl::ReadGrad<double>(_QSpaceOrientations,DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
    
    utlGlobalException(_BFileArg.isSet() && _QSpaceOrientationsArg.isSet() && bvec->size()!=qMatrix->Rows(), "The size in --bfile and --qorientations are not the same");
    
    if (bvec->size()>0 && qMatrix->Rows()>0 && bvec->size()!=qMatrix->Rows())
      utl::MatchBVectorAndGradientMatrix(*bvec, *qMatrix);
    dwiGenerator->GetSamplingSchemeQSpace()->SetOrientationsCartesian( qMatrix );
    dwiGenerator->GetSamplingSchemeQSpace()->SetBVector(bvec);
    }

  if (_RSpaceOrientationsArg.isSet())
    {
    utlGlobalException(!_RValuesArg.isSet() && !_RFileArg.isSet(), "need to set --rvalues or --rfile");
    utlGlobalException(_RValuesArg.isSet() && _RFileArg.isSet(), "Only one of these is allowed: --rvalues or --rfile");

    if (_RValuesArg.isSet())
      *rvec = _RValues;
    if (_RFileArg.isSet())
      utl::ReadVector(_RFile,*rvec);
    
    GeneratorType::MatrixPointer rMatrix =  utl::ReadGrad<double>(_RSpaceOrientations,DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
    
    utlGlobalException(_RFileArg.isSet() && _RSpaceOrientationsArg.isSet() && rvec->size()!=rMatrix->Rows(), "The size in --rfile and --rorientations are not the same");

    if (rvec->size()>0 && rMatrix->Rows()>0 && rvec->size()!=rMatrix->Rows())
      utl::MatchBVectorAndGradientMatrix(*rvec, *rMatrix);
    dwiGenerator->GetSamplingSchemeRSpace()->SetOrientationsCartesian( rMatrix );
    dwiGenerator->GetSamplingSchemeRSpace()->SetRadiusVector(rvec);
    }


  if (_BackgroundDiffusionParametersArg.isSet())
    dwiGenerator->SetBackgroundDiffusionParameterValues(_BackgroundDiffusionParameters);

  utlGlobalException(_CylinderParameters.size()!=3, "wrong size in --cylinderparam");
  GeneratorType::CylinderModelPointer cylinder = GeneratorType::CylinderModelType::New();
  cylinder->SetLength(_CylinderParameters[0]);
  cylinder->SetRadius(_CylinderParameters[1]);
  cylinder->SetD0(_CylinderParameters[2]);
  dwiGenerator->SetCylinderModel(cylinder);

  dwiGenerator->Update();


  // Write Output 

  if (_OutputDWIFile!="")
    {
    std::cout << "Writing DWI data ... " << std::endl;
    outputImage = dwiGenerator->GetDWIImage();
    if (_OutputDWIType!="4D")
      {
      typedef itk::DWIWriter<float> DWIWriterType;
      DWIWriterType::Pointer writer = DWIWriterType::New();
      writer->SetInput(outputImage);
      if (_OutputDWIType=="EACHSHELL")
        writer->SetOutputEachShell(true);
      writer->SetSamplingSchemeQSpace(dwiGenerator->GetSamplingSchemeQSpace());

      std::string outputFile_noext, outputFile_ext;
      utl::GetFileExtension(_OutputDWIFile, outputFile_ext, outputFile_noext);
      writer->SetConfigurationFile(outputFile_noext + ".txt");
      writer->SetBFile(outputFile_noext + "_b.txt");
      writer->SetOrientationFile(outputFile_noext+"_grad.txt");
      writer->SetDWIFile(outputFile_noext + ".nii.gz");
      writer->Update();
      }
    else
      {
      itk::SaveImage<OutputImageType>(outputImage, _OutputDWIFile);
      }

    if ( _B0FileArg.isSet() )
      itk::SaveImage<ScalarImageType>(dwiGenerator->GetB0Image(), _B0File);
    }

  if (_OutputEAPFile!="")
    {
    std::cout << "Writing EAP data ... " << std::endl;
    outputImage = dwiGenerator->GetEAPImage();
    utlException(itk::IsImageEmpty(outputImage), "logical error! eap image does not exist. Need to set --rorientations --rvalues");
    itk::SaveImage<OutputImageType>(outputImage, _OutputEAPFile);
    }
  
  if (_OutputODFFile!="")
    {
    std::cout << "Writing ODF data ... " << std::endl;
    outputImage = dwiGenerator->GetODFImage();
    utlException(itk::IsImageEmpty(outputImage), "logical error! odf image does not exist. Need to set --rorientations --rvalues");
    itk::SaveImage<OutputImageType>(outputImage, _OutputODFFile);
    }

  if (_OutputPeakFile!="")
    {
    std::cout << "Writing Peak data ... " << std::endl;
    utlSAGlobalException(_MaxNumberOfPeaks<=0).msg("need to set _MaxNumberOfPeaks in --peaks");
    outputImage = dwiGenerator->GetPeakImage();
    utlException(itk::IsImageEmpty(outputImage), "logical error! peak image does not exist");
    itk::SaveImage<OutputImageType>(outputImage, _OutputPeakFile);
    }
  
  if (_OutputRTOFile!="")
    {
    std::cout << "Writing RTO data ... " << std::endl;
    outputScalarImage = dwiGenerator->GetRTOImage();
    utlException(itk::IsImageEmpty(outputScalarImage), "logical error! rto image does not exist. Need to set --rorientations --rvalues");
    itk::SaveImage<ScalarImageType>(outputScalarImage, _OutputRTOFile);
    }

  if (_OutputMSDFile!="")
    {
    std::cout << "Writing MSD data ... " << std::endl;
    outputScalarImage = dwiGenerator->GetMSDImage();
    utlException(itk::IsImageEmpty(outputScalarImage), "logical error! msd image does not exist. Need to set --rorientations --rvalues");
    itk::SaveImage<ScalarImageType>(outputScalarImage, _OutputMSDFile);
    }
    
  return EXIT_SUCCESS;
}

