/**
 *       @file  SamplingSchemeQSpaceIncrementalEstimation.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-15-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "itkSamplingSchemeQSpaceIncrementalEstimationFilter.h"
#include "utl.h"
#include "itkSamplingSchemeQSpaceWriter.h"
#include "itkSamplingSchemeQSpace.h"

#include "SamplingSchemeQSpaceIncrementalEstimationCLP.h"

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
  typedef itk::SamplingSchemeQSpaceIncrementalEstimationFilter<SamplingType> EstimationType;
  EstimationType::Pointer estimator = EstimationType::New();

  utlGlobalException(_NumberOfSamples.size()==0, "should give the numbers of samples in shells");
  EstimationType::IndexVectorType numbers(_NumberOfSamples.size());
  for ( unsigned int i = 0; i < _NumberOfSamples.size(); i += 1 ) 
    numbers[i] = _NumberOfSamples[i];

  EstimationType::MatrixPointer fineOrientations(new EstimationType::MatrixType());

  estimator->SetDebug(_Debug);
  estimator->SetNumbersInShell(numbers);
  estimator->SetTessellationOrder(_Order);
  if (_Criteria=="DISTANCE")
    estimator->SetCriteriaType(EstimationType::DISTANCE);
  else if (_Criteria=="ELECTROSTATIC")
    estimator->SetCriteriaType(EstimationType::ELECTROSTATIC);
  if(_FineOrientationsArg.isSet())
    {
    fineOrientations = utl::ReadGrad<double>(_FineOrientations, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
    estimator->SetFineOrientations(fineOrientations);
    }
  if(_InitialOrientationsArg.isSet())
    {
    SamplingPointer initialSamples = SamplingType::New();
    initialSamples->ReadOrientationFile(_InitialOrientations, DIRECTION_NODUPLICATE);
    estimator->SetInitialOrientations(initialSamples);
    }
  estimator->SetElectrostaticOrder(_EnergyOrder);
  estimator->SetWeightForSingleShell(_Weight);

  itk::TimeProbe clock;
  clock.Start();
  estimator->UpdateOutputData();
  clock.Stop();
  std::cout << "Elapsed time : " << clock.GetTotal() << std::endl;

  SamplingPointer output = estimator->GetOutputOrientations();
  // output->Print(std::cout<<"output");
  // output->RemoveSamplesNotIndexed();

  typedef itk::SamplingSchemeQSpaceWriter<SamplingType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetSampling(output);
  writer->SetOrientationFile(_OutputOrientations);
  writer->SaveSingleShellOn();
  writer->SaveAllShellsInOneFileOff();
  writer->Update();
  
  return 0;
}

