/**
 *       @file  DWISingleVoxelSimulator.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-24-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkDWISingleVoxelGenerator.h"
#include "DWISingleVoxelSimulatorCLP.h"
#include "itkDWIWriter.h"

/**
 * \brief  Simulate DWI data with fixed single voxel configuration.
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  typedef itk::VectorImage<float, 3> OutputImageType;
  OutputImageType::Pointer outputImage;
  typedef itk::Image<float, 3> B0ImageType;
  
  typedef itk::DWISingleVoxelGenerator<OutputImageType, B0ImageType> GeneratorType;
  GeneratorType::Pointer dwiGenerator = GeneratorType::New();

  utlGlobalException(_DiffusionParameters.size()==0, "need to set --params");
  dwiGenerator->SetDiffusionParameterValues(_DiffusionParameters);

  dwiGenerator->SetDebug(_DebugArg.isSet());
  utlGlobalException(_OutputDWIFile=="" && _OutputODFFile=="" && _OutputEAPFile=="" 
    && _OutputMSDFile=="" && _OutputRTOFile=="", "need to set output dwi/odf/eap/peak/rto/msd");
  dwiGenerator->SetIsOutputDWI(_OutputDWIFile!="");
  dwiGenerator->SetIsOutputEAP(_OutputEAPFile!="");
  dwiGenerator->SetIsOutputODF(_OutputODFFile!="");
  dwiGenerator->SetIsOutputRTO(_OutputRTOFile!="");
  dwiGenerator->SetIsOutputMSD(_OutputMSDFile!="");
  
  if ( _SNRArg.isSet() )
    dwiGenerator->SetSNR( _SNR );
  if (_NoiseSigmaArg.isSet())
    dwiGenerator->SetNoiseSigma( _NoiseSigma );
  utlGlobalException(_SNRArg.isSet() && _NoiseSigmaArg.isSet(), "Only one of these is allowed: --snr or --noisesigma");
  
  if (_TauArg.isSet())
    {
    dwiGenerator->GetSamplingSchemeQSpace()->SetTau(_Tau);
    dwiGenerator->GetSamplingSchemeRSpace()->SetTau(_Tau);
    }
  
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

    if (rvec->size()>0 && rMatrix->Rows()>0 && rvec->size()!=rMatrix->Rows())
      utl::MatchBVectorAndGradientMatrix(*rvec, *rMatrix);
    dwiGenerator->GetSamplingSchemeRSpace()->SetOrientationsCartesian( rMatrix );
    dwiGenerator->GetSamplingSchemeRSpace()->SetRadiusVector(rvec);
    }
  
  if (_B0ScaleArg.isSet())
    dwiGenerator->SetB0Scale(_B0Scale);

  if (_ModelType=="TENSOR_SPHERICAL")
    dwiGenerator->SetModelType(GeneratorType::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS);
  else if (_ModelType=="TENSOR_CARTESIAN")
    dwiGenerator->SetModelType(GeneratorType::SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS);
  else if (_ModelType=="CYLINDER_SPHERICAL")
    dwiGenerator->SetModelType(GeneratorType::CYLINDER_SPHERICAL_MODEL);

  if (_RandomType=="FIXED")
    dwiGenerator->SetRandomType(GeneratorType::FIXED);
  else if (_RandomType=="UNIFORM")
    {
    dwiGenerator->SetRandomType(GeneratorType::UNIFORM);
    utlSAGlobalException(_TessellationOrder>0 && _ElecNumber>0).msg("only one of --elecnumber and --tessorder can be set");
    if (_TessellationOrder>0)
      dwiGenerator->SetUniformFromTessellation(_TessellationOrder);
    if (_ElecNumber>0)
      dwiGenerator->SetUniformFromElectrostaticEnergy(_ElecNumber);
    }

  utlGlobalException(_Size.size()!=3, "--size should have 3 dimension");
  OutputImageType::SizeType size;
  for ( int i = 0; i < 3; i += 1 ) 
    size[i] = _Size[i];
  dwiGenerator->SetOutputSize(size);

  if (_MaxNumberOfPeaks>0)
    {
    dwiGenerator->SetMaxNumberOfPeaks(_MaxNumberOfPeaks);
    dwiGenerator->SetPeakType(itk::PeakContainerHelper::GetPeakType(_PeakType));
    }
  
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
      itk::SaveImage<B0ImageType>(dwiGenerator->GetB0Image(), _B0File);
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
    
  return EXIT_SUCCESS;

}

