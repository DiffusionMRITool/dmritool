/**
 *       @file  SamplingSchemeQSpaceIMOCEstimation.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "10-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#include "itkSamplingSchemeQSpaceIMOCEstimationFilter.h"
#include "utl.h"
#include "itkSamplingSchemeQSpaceWriter.h"
#include "itkSamplingSchemeQSpace.h"

#include "SamplingSchemeQSpaceIMOCEstimationCLP.h"

/**
 * \brief  incrementatl estimation for single and multiple shell sampling scheme
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  
  PARSE_ARGS;

  typedef itk::SamplingSchemeQSpace<double>  SamplingType;
  typedef SamplingType::Pointer  SamplingPointer;
  typedef itk::SamplingSchemeQSpaceIMOCEstimationFilter<SamplingType> EstimationType;
  EstimationType::Pointer estimator = EstimationType::New();

  utlGlobalException(_NumberOfSamples.size()==0, "should give the numbers of samples in shells");
  EstimationType::IndexVectorType numbers(_NumberOfSamples.size());
  for ( unsigned int i = 0; i < _NumberOfSamples.size(); i += 1 ) 
    numbers[i] = _NumberOfSamples[i];

  EstimationType::MatrixPointer fineOrientations(new EstimationType::MatrixType());

  estimator->SetDebug(_Debug);
  estimator->SetNumbersInShell(numbers);
  estimator->SetTessellationOrder(_Order);
  if(_FineOrientationsArg.isSet())
    {
    fineOrientations = utl::ReadGrad<double>(_FineOrientations, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
    estimator->SetFineOrientations(fineOrientations);
    }
  estimator->SetWeightForSingleShell(_Weight);
  estimator->SetAngleMinChange(_MinChange);
  estimator->SetChooseMinimalCoverageShell(true);

  itk::TimeProbe clock;
  clock.Start();
  estimator->UpdateOutputData();
  clock.Stop();
  std::cout << "Elapsed time : " << clock.GetTotal() << std::endl;

  SamplingPointer output = estimator->GetOutputOrientations();
  // output->Print(std::cout<<"output");

  typedef itk::SamplingSchemeQSpaceWriter<SamplingType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetSampling(output);
  writer->SetOrientationFile(_OutputOrientations);
  writer->SaveSingleShellOn();
  writer->SaveAllShellsInOneFileOff();
  writer->Update();
  
  return 0;
}

