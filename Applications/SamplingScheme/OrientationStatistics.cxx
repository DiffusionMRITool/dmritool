/**
 *       @file  OrientationStatistics.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "07-25-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "OrientationStatisticsCLP.h"
#include "itkSamplingScheme3D.h"

#include "utl.h"

inline double 
processOrientation ( const itk::SamplingScheme3D<double>::Pointer scheme, const bool _Asymmetric )
{
  typedef itk::SamplingScheme3D<double> SamplingType;
  // scheme->Print(std::cout<<"scheme");
  
  int size = scheme->GetNumberOfSamples(); 
  std::cout << "size = " << size << std::endl << std::flush;

  int numberUniqueSamples=0, numberAntipodalSamples=0, numberRepeatedSamples=0;
  scheme->GetNumbers(numberUniqueSamples, numberAntipodalSamples, numberRepeatedSamples);
  std::cout << numberUniqueSamples << " unique samples, " << numberAntipodalSamples << " antipodal samples, " << numberRepeatedSamples << " repeated samples."  << std::endl << std::flush;

  std::vector<double> angleVec = scheme->CalculateMinDistance(!_Asymmetric);
  std::vector<double> stats = utl::GetContainerStats(angleVec.begin(), angleVec.end());
  std::cout << "minimal angle = " << stats[0]*180.0/M_PI << ", radian=" << stats[0] << " (covering radius)" << std::endl << std::flush;
  std::cout << "maximal angle = " << stats[1]*180.0/M_PI << ", radian=" << stats[1] << std::endl << std::flush;
  std::cout << "mean angle = " << stats[2]*180.0/M_PI << std::endl << std::flush;
  std::cout << "std angle = " << stats[3]*180.0/M_PI << std::endl << std::flush;

  double ud2 = SamplingType::CalculateMinDistanceUpperBound(2*size);
  double ud = SamplingType::CalculateMinDistanceUpperBound(size);
  std::cout << "upper bound ("<< 2*size <<" points) = " << ud2*180/M_PI << ", radian=" << ud2 << std::endl << std::flush;
  std::cout << "upper bound ("<< size <<" points) = " << ud*180/M_PI << ", radian=" << ud << std::endl << std::flush;

  if (numberRepeatedSamples==0)
    {
    // double isSymmetric = numberAntipodalSamples==0 && !_Asymmetric;

    double electrostaticEnergy = scheme->CalculateElectrostaticEnergy(2.0, !_Asymmetric);
    std::cout << "electrostaticEnergy (order=2) = " << electrostaticEnergy << std::endl << std::flush;
    electrostaticEnergy = scheme->CalculateElectrostaticEnergy(1.0, !_Asymmetric);
    std::cout << "electrostaticEnergy (order=1) = " << electrostaticEnergy << std::endl << std::flush;

    double packing_density = scheme->CalculatePackingDensity(!_Asymmetric);
    std::cout << "Spherical packing density = " << packing_density << std::endl << std::flush;
    }

  return stats[0];
}

/**
 * \brief  obtain energies, minimal angles from a given set of orientations
 *
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;
  
  typedef itk::SamplingScheme3D<double> SamplingType;
  SamplingType::Pointer schemeAll = SamplingType::New();

  double cost=0;
  for ( int i = 0; i < _InputOrientationFile.size(); i += 1 ) 
    {
    std::cout << "file: " << _InputOrientationFile[i] << std::endl << std::flush;
    SamplingType::Pointer scheme = SamplingType::New();
    scheme->ReadOrientationFile(_InputOrientationFile[i]);
    // scheme->Print(std::cout<<"scheme=");
    cost += _Weight/_InputOrientationFile.size() * processOrientation(scheme,_Asymmetric);
    std::cout << std::endl << std::flush;
    if (_Combine && _InputOrientationFile.size()>1)
      {
      for ( int j = 0; j < scheme->size(); j += 1 ) 
        schemeAll->AppendOrientation((*scheme)[j]);
      }
    }

  if (_Combine && _InputOrientationFile.size()>1)
    {
    std::cout << "Combine " << _InputOrientationFile.size() << " orientations" << std::endl << std::flush;
    cost += (1-_Weight) * processOrientation(schemeAll,_Asymmetric);
    }

  std::cout << "Spherical code cost function = " << cost << std::endl << std::flush;

  return 0;
}
