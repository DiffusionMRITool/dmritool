/**
 *       @file  CorrectBValues.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-01-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include "utlCore.h"
#include "CorrectBValuesCLP.h"


/**
 * \brief  b values whose distance is smallter the threshold will be considered in 
 *  the same shell, and they will be replaced as their mean b value. 
 *
 *  \author  Jian Cheng
 */
int 
main (int argc, char const* argv[])
{

  // GenerateCLP
  PARSE_ARGS;

  typedef std::vector<double>  STDVectorType;
  STDVectorType bVector;
  std::cout << "Read b values from " << _InputFile << std::endl << std::flush;
  utl::ReadVector(_InputFile,bVector);

  std::vector<STDVectorType> bVectors; 
  std::vector<std::vector<int> > bIndices;
  STDVectorType bMax, bMin;
  for ( int i = 0; i < bVector.size(); i += 1 ) 
    {
    double b = bVector[i];
    int j=0;
    for ( j = 0; j < bVectors.size(); j += 1 ) 
      {
      if (b>=bMin[j]-_BThreshold && b<=bMax[j]+_BThreshold)
        {
        bVectors[j].push_back(b);
        bIndices[j].push_back(i);
        if (b<bMin[j])
          bMin[j] = b;
        else if (b>bMax[j])
          bMax[j] = b;
        break;
        }
      }
    if (j==bVectors.size())
      {
      STDVectorType bVecTemp;
      bVecTemp.push_back(b);
      bVectors.push_back(bVecTemp);
      std::vector<int> bIndexTemp;
      bIndexTemp.push_back(i);
      bIndices.push_back(bIndexTemp);
      bMin.push_back(b);
      bMax.push_back(b);
      }
    }

  for ( int j = 0; j < bVectors.size(); j += 1 ) 
    {
    utlAssert(bMax[j]-bMin[j]<100.0, "the range of b values is larger than 100, which can not be in the same shell");
    std::cout << "shell " << j << ": bMin=" << bMin[j] << " ,bMax=" << bMax[j] << ", " << bVectors[j].size() << " b values." << std::endl << std::flush;
    double bMean=0;
    for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
      bMean += bVectors[j][k];
    bMean /= bVectors[j].size();
    if (bMean<_BThreshold)
      bMean=0.0;
    for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
      {
      bVector[bIndices[j][k] ] = bMean;
      }
    }

  std::cout << "Save corrected b values to " << _OutputFile << std::endl << std::flush;
  utl::SaveVector(bVector, _OutputFile);

  return 0;
}
