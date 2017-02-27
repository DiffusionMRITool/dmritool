/**
 *       @file  SamplingSchemeQSpaceRandomEstimation.cxx
 *      @brief  
 *     Created  "01-17-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#include "utl.h"
#include "itkSamplingSchemeQSpaceWriter.h"
#include "itkSamplingSchemeQSpace.h"

#include "SamplingSchemeQSpaceRandomEstimationCLP.h"

/**
 * \brief  randome generator for single and multiple shell sampling scheme
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  
  PARSE_ARGS;

  typedef itk::SamplingSchemeQSpace<double>  SamplingType;
  typedef SamplingType::Pointer  SamplingPointer;

  utlGlobalException(_NumberOfSamples.size()==0, "should give the numbers of samples in shells");

  SamplingPointer scheme = SamplingType::New();
  scheme->GenerateFromRandomPoints(_NumberOfSamples);

  typedef itk::SamplingSchemeQSpaceWriter<SamplingType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetSampling(scheme);
  writer->SetOrientationFile(_OutputOrientations);
  writer->SaveSingleShellOn();
  writer->SaveAllShellsInOneFileOff();
  writer->Update();
  
  return 0;
}

