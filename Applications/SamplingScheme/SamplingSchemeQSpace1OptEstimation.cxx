/**
 *       @file  SamplingSchemeQSpace1OptEstimation.cxx
 *      @brief  
 *     Created  "01-17-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "itkSamplingSchemeQSpace1OptEstimationFilter.h"
#include "utl.h"
#include "itkSamplingSchemeQSpaceWriter.h"
#include "itkSamplingSchemeQSpace.h"

#include "SamplingSchemeQSpace1OptEstimationCLP.h"

/**
 * \brief  incrementatl estimation for single and multiple shell sampling scheme
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  
  PARSE_ARGS;
  
  utl::LogLevel = _Verbose;

  typedef itk::SamplingSchemeQSpace<double>  SamplingType;
  typedef SamplingType::Pointer  SamplingPointer;
  typedef itk::SamplingSchemeQSpace1OptEstimationFilter<SamplingType> EstimationType;
  EstimationType::Pointer estimator = EstimationType::New();

  utlGlobalException(_InitialOrientations.size()==0, "need to set --initial");
  SamplingPointer initialSamples = SamplingType::New();
  initialSamples->ReadOrientationFileList(_InitialOrientations, DIRECTION_NODUPLICATE);

  EstimationType::MatrixPointer fineOrientations(new EstimationType::MatrixType());

  estimator->SetDebug(_Verbose>=LOG_DEBUG);
  // estimator->SetNumbersInShell(numbers);
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
  estimator->SetInitialOrientations(initialSamples);
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
