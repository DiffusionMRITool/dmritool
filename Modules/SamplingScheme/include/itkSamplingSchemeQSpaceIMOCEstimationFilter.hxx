/**
 *       @file  itkSamplingSchemeQSpaceIMOCEstimationFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "10-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpaceIMOCEstimationFilter_hxx
#define __itkSamplingSchemeQSpaceIMOCEstimationFilter_hxx


#include "itkSamplingSchemeQSpaceIMOCEstimationFilter.h"
#include "utl.h"

#include <numeric>

namespace itk
{

template< class TSamplingType >
SamplingSchemeQSpaceIMOCEstimationFilter< TSamplingType >
::SamplingSchemeQSpaceIMOCEstimationFilter() 
  : Superclass(),
  m_FineOrientations(new MatrixType()) 
{
}

template< class TSamplingType >
void
SamplingSchemeQSpaceIMOCEstimationFilter< TSamplingType >
::Initialization()
{
  // SamplingType* input = const_cast< SamplingType * >( this->GetInput() );
  if (this->m_FineOrientations->Size()==0)
    {
    this->m_FineOrientations = utl::ReadGrad<double>(this->m_TessellationOrder, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); // catesian
    }
  m_FineScheme->SetOrientationsCartesian(m_FineOrientations);
  // utl::PrintUtlMatrix(*m_FineOrientations, "m_FineOrientations");

  // assume that m_FineScheme is uniformly distributed in sphere
  m_MinDistanceInFineScheme = m_FineScheme->CalculateMinDistance(0,true);
  m_MinDistanceInFineScheme = std::sqrt(2-2*std::cos(m_MinDistanceInFineScheme));

  // build kd tree
  m_Sample = SampleType::New();
  m_Sample->SetMeasurementVectorSize( 3 );
  unsigned int numberOfSamples = m_FineOrientations->Rows();
  MeasurementVectorType mv;
  for (unsigned int i = 0; i < numberOfSamples; ++i )
    {
    mv = (*m_FineScheme)[i];
    m_Sample->PushBack( mv );
    for ( int j = 0; j < 3; ++j ) 
      mv[j] = -mv[j];
    m_Sample->PushBack( mv );
    }

  m_TreeGenerator = TreeGeneratorType::New();
  m_TreeGenerator->SetSample( m_Sample );
  m_TreeGenerator->SetBucketSize( 16 );
  m_TreeGenerator->Update();

  // m_TreeGenerator->Print(std::cout<<"treeGenerator");
  m_KDTree = m_TreeGenerator->GetOutput();
}

template< class TSamplingType >
bool
SamplingSchemeQSpaceIMOCEstimationFilter< TSamplingType >
::IsSatisfiedSeparationAngles(const std::vector<double>& angles)
{
  std::vector<double> euclideanDistances(angles.size()), overlapBallsInnerProduct(angles.size());
  for ( int i = 0; i < angles.size(); ++i ) 
    {
    overlapBallsInnerProduct[i] = std::fabs(std::cos(2*angles[i]));
    euclideanDistances[i] = std::sqrt(2-2*std::cos(angles[i]));
    }

  unsigned int numberOfShells = this->m_NumbersInShell.size();
  unsigned int numberOfSamples = m_FineScheme->GetNumberOfSamples();
  int totalNumberOfSamples = std::accumulate(this->m_NumbersInShell.begin(), this->m_NumbersInShell.end(), 0);
  Index2DVectorType sumOverlapVec(angles.size(), IndexVectorType(numberOfSamples,-1));

  std::vector<int> coverageInShells(numberOfShells,0);
  bool isSameAngle = true;
  for ( int j = 1; j < numberOfShells; ++j ) 
    {
    if (std::fabs(angles[0]-angles[j])>1e-4)
      {
      isSameAngle = false;
      break;
      }
    }
  bool needCoverageInshells = m_ChooseMinimalCoverageShell && !isSameAngle;

  MeasurementVectorType queryPoint, queryCandidate;
  typedef typename TreeType::InstanceIdentifierVectorType InstanceIdentifierVectorType;
  InstanceIdentifierVectorType neighbors, neighbors_1, candidates, candidates_1, neighbors_jj;
  InstanceIdentifierVectorType candidateNumbersInOverlap;
    
  Index2DVectorPointer indices = this->m_FineScheme->GetIndicesInShells();
  indices->clear();
    
  Index2DVectorType chosenIndices(numberOfShells);

  typedef utl_unordered_map<int, InstanceIdentifierVectorType> NeighborMapType;
  typedef typename NeighborMapType::const_iterator NeighborMapIterator;
  NeighborMapType neighborMap;
  int indexDistPair;

  if (numberOfShells==1)
    {
    IndexVectorType hasChosen(numberOfSamples,-1);
    std::set<int> nextCandidateIndices;

    int currentIndex = -1;
    if (chosenIndices[0].size()==0)
      currentIndex = 0;

    while (chosenIndices[0].size()< this->m_NumbersInShell[0] )
      {
      // std::cout << "\n chosenIndices[0.size()=" << chosenIndices[0.size() << std::endl << std::flush;

      hasChosen[currentIndex]=0;
      chosenIndices[0].push_back(currentIndex);
      queryPoint = (*m_FineScheme)[currentIndex];
      // std::cout << "queryPoint = " << queryPoint << std::endl << std::flush;
      
      indexDistPair=currentIndex;

      NeighborMapIterator mapIter = neighborMap.find(indexDistPair);
      if (mapIter==neighborMap.end())
        {
        m_KDTree->Search( queryPoint, euclideanDistances[0], neighbors );
        neighborMap[indexDistPair] = neighbors;
        }
      else
        neighbors = mapIter->second;
      for ( int i = 0; i < neighbors.size(); ++i ) 
        hasChosen[neighbors[i]/2] = 0;

      m_KDTree->Search( queryPoint, euclideanDistances[0]+1*m_MinDistanceInFineScheme, candidates );
      utlException(candidates.size()<=neighbors.size(), "the size of candidates should be larger than the size of neighbors.");
      // utlPrintVar2(true, neighbors.size(), candidates.size());
      for ( int i = 0; i < candidates.size(); ++i ) 
        {
        int jj = candidates[i]/2;
        if (hasChosen[jj] == -1)
          nextCandidateIndices.insert(jj);
        }

      int indexMaxOverlap=-1;
      int sumMaxOverlap=-1;
      int numberOfCandidates=0, numberSearch=0;
      // utl::Tic(std::cout<<"time");
// #pragma omp parallel shared(indexMaxOverlap, sumMaxOverlap) private(iter, neighbors_jj, queryCandidate) firstprivate(indexMaxOverlap_ii, sumMaxOverlap_ii)
      for ( typename std::set<int>::iterator iter = nextCandidateIndices.begin(); iter!=nextCandidateIndices.end(); ++iter ) 
        {
        int jj = *iter;
        if (hasChosen[jj]==-1)
          {
          queryCandidate = (*m_FineScheme)[jj];
          bool overlapBall = true;
          double innerProduct = 0;
          for ( int kk = 0; kk < 3; ++kk ) 
            innerProduct+= queryPoint[kk]*queryCandidate[kk];
          if (std::fabs(innerProduct) < overlapBallsInnerProduct[0])
            overlapBall = false;
            
          int sumOverlap = 0;
          if (sumOverlapVec[0][jj]>=0 && !overlapBall)
            sumOverlap = sumOverlapVec[0][jj];
          else
            {
          //  sumOverlapVec[0][jj]<0 || (sumOverlapVec[0][jj]>=0 && overlapBall) 
            // m_KDTree->Search( queryCandidate, euclideanDistances[0], neighbors_jj );
            indexDistPair=jj;
            NeighborMapIterator mapIter = neighborMap.find(indexDistPair);
            if (mapIter==neighborMap.end())
              {
              m_KDTree->Search( queryCandidate, euclideanDistances[0], neighbors_jj );
              neighborMap[indexDistPair] = neighbors_jj;
              }
            else
              neighbors_jj = mapIter->second;
            numberSearch++;
            for ( int k = 0; k < neighbors_jj.size(); ++k ) 
              {
              if (hasChosen[neighbors_jj[k]/2]==0)
                sumOverlap++;
              }
            sumOverlapVec[0][jj]=sumOverlap;
            }

          if (sumMaxOverlap<sumOverlap)
            {
            indexMaxOverlap = jj;
            sumMaxOverlap = sumOverlap;
            }
          numberOfCandidates++;
          }
        }
      // utl::Toc();
      // utlPrintVar2(true, numberOfCandidates, numberSearch);

      if (numberOfCandidates==0)
        break;

      currentIndex = indexMaxOverlap;
      }

    if (chosenIndices[0].size()==this->m_NumbersInShell[0])
      {
      *indices = chosenIndices;
      return true;
      }
    else
      return false;
    }
  else
    {

    // utl::PrintVector(this->m_NumbersInShell, "this->m_NumbersInShell");
    // utlPrintVar1(true, totalNumberOfSamples);
    // utl::PrintVector(euclideanDistances, "euclideanDistances 11");
    // utl::PrintVector(angles, "angles 11");
    for ( int s = 0; s < numberOfShells; ++s ) 
      {
      // It is impossible that the seperation angle in the combined shell is larger than the seperation angle in individual shell.
      if (euclideanDistances[s]<= euclideanDistances.back())
        return false;
      }

    Index2DVectorType hasChosen(numberOfSamples, IndexVectorType(angles.size(),-1));
    std::vector<std::set<int> > nextCandidateIndices(numberOfShells+1);

    // the first sample in each shell
    int currentIndex = 0;
    for ( int s = 0; s < numberOfShells; ++s ) 
      {
      hasChosen[currentIndex][s]=0;
      hasChosen[currentIndex][numberOfShells]=s;

      chosenIndices[s].push_back(currentIndex);
      queryPoint = (*m_FineScheme)[currentIndex];

      // utlPrintVar1(true, currentIndex);
      // std::cout << "queryPoint=" << queryPoint << std::endl << std::flush;
      
      indexDistPair=currentIndex + s*numberOfSamples;
      NeighborMapIterator mapIter = neighborMap.find(indexDistPair);
      if (mapIter==neighborMap.end())
        {
        m_KDTree->Search( queryPoint, euclideanDistances[s], neighbors );
        neighborMap[indexDistPair] = neighbors;
        }
      else
        neighbors = mapIter->second;
      for ( int i = 0; i < neighbors.size(); ++i ) 
        hasChosen[neighbors[i]/2][s] = 0;
      // update coverages
      if (needCoverageInshells)
        coverageInShells[s] += neighbors.size();

      indexDistPair=currentIndex + numberOfShells*numberOfSamples;
      mapIter = neighborMap.find(indexDistPair);
      if (mapIter==neighborMap.end())
        {
        m_KDTree->Search( queryPoint, euclideanDistances.back(), neighbors_1 );
        neighborMap[indexDistPair] = neighbors_1;
        }
      else
        neighbors_1 = mapIter->second;
      for ( int i = 0; i < neighbors_1.size(); ++i ) 
        hasChosen[neighbors_1[i]/2][numberOfShells] = s;
            
      m_KDTree->Search( queryPoint, euclideanDistances[s]+1*m_MinDistanceInFineScheme, candidates );
      utlException(candidates.size()<=neighbors.size(), "the size of candidates should be larger than the size of neighbors.");

      m_KDTree->Search( queryPoint, euclideanDistances.back()+1*m_MinDistanceInFineScheme, candidates_1 );
      utlException(candidates_1.size()<=neighbors_1.size(), "the size of candidates should be larger than the size of neighbors.");

      for ( int i = 0; i < candidates.size(); ++i ) 
        {
        int jj = candidates[i]/2;
        if ( hasChosen[jj][s]==-1 && hasChosen[jj][numberOfShells]==-1)
          nextCandidateIndices[s].insert(jj);
        }

      for ( int i = 0; i < candidates_1.size(); ++i ) 
        {
        int jj = candidates_1[i]/2;
        if ( hasChosen[jj][numberOfShells]==-1)
          {
          for ( int sss = 0; sss < numberOfShells; ++sss ) 
            {
            if (sss!=s && hasChosen[jj][sss]==-1)
              nextCandidateIndices[sss].insert(jj);
            }
          nextCandidateIndices[numberOfShells].insert(jj);
          }
        }

      int indexMaxOverlap=-1;
      int sumMaxOverlap=-1;
      int numberOfCandidates=0, numberSearch=0;

      for ( typename std::set<int>::iterator iter = nextCandidateIndices[numberOfShells].begin(); iter!=nextCandidateIndices[numberOfShells].end(); ++iter ) 
        {
        int jj = *iter;
        if (hasChosen[jj][numberOfShells]==-1)
          {
          queryCandidate = (*m_FineScheme)[jj];
          bool overlapBall = true;
          double innerProduct = 0;
          for ( int kk = 0; kk < 3; ++kk ) 
            innerProduct+= queryPoint[kk]*queryCandidate[kk];
          if (std::fabs(innerProduct) < overlapBallsInnerProduct[numberOfShells])
            overlapBall = false;

          int sumOverlap = 0;
          if (sumOverlapVec[numberOfShells][jj]>=0 && !overlapBall)
            sumOverlap = sumOverlapVec[numberOfShells][jj];
          else
            {
          //  sumOverlapVec[numberOfShells][jj]<0 || (sumOverlapVec[numberOfShells][jj]>=0 && overlapBall) 
            // m_KDTree->Search( queryCandidate, euclideanDistances.back(), neighbors_jj );
            indexDistPair = jj + numberOfShells*numberOfSamples;
            mapIter = neighborMap.find(indexDistPair);
            if (mapIter==neighborMap.end())
              {
              m_KDTree->Search( queryCandidate, euclideanDistances.back(), neighbors_jj );
              neighborMap[indexDistPair] = neighbors_jj;
              }
            else
              neighbors_jj = mapIter->second;
            numberSearch++;
            for ( int k = 0; k < neighbors_jj.size(); ++k ) 
              {
              if (hasChosen[neighbors_jj[k]/2][numberOfShells]>=0)
                sumOverlap++;
              }
            sumOverlapVec[numberOfShells][jj]=sumOverlap;
            }

          if (sumMaxOverlap<sumOverlap)
            {
            indexMaxOverlap = jj;
            sumMaxOverlap = sumOverlap;
            }
          numberOfCandidates++;
          }
        }

      currentIndex = indexMaxOverlap;
      }

    // for ( int i = 0; i < chosenIndices.size(); ++i ) 
    //   utl::PrintVector(chosenIndices[i], "chosenIndices[i]");

    // all other samples 
    MeasurementVectorType queryPointLast(queryPoint);
    int indexShellMaxOverlap=-1, indexShellMaxOverlapLast=-1, shellIndexMinCoverage=-1;

    for ( unsigned int n = numberOfShells; n < totalNumberOfSamples; ++n ) 
      {
      // std::cout << "n = " << n << std::endl << std::flush;

      int indexMaxOverlap=-1;
      int sumMaxOverlap=-1;
      int numberOfCandidates=0, numberSearch=0;

      if (needCoverageInshells)
        {
        // set the coverage of the shell with requested number of samples as the largest value.
        for ( int s = 0; s < numberOfShells; ++s ) 
          {
          if (chosenIndices[s].size() >= this->m_NumbersInShell[s])
            coverageInShells[s] = numberOfSamples;
          }
        // find the shell with the minimal coverage and it has not been set requested number of samples.
        shellIndexMinCoverage = utl::argmin(coverageInShells.begin(), coverageInShells.end());
        }

      for ( int s = 0; s < numberOfShells; ++s ) 
        {

        if (needCoverageInshells && shellIndexMinCoverage!=s)
          {
          // only set the current sample to the shell with minimal coverage.
          continue;
          }

        if (chosenIndices[s].size() < this->m_NumbersInShell[s])
          {

          for ( typename std::set<int>::iterator iter = nextCandidateIndices[s].begin(); iter!=nextCandidateIndices[s].end(); ++iter ) 
            {
            int jj = *iter;
            if (hasChosen[jj][s]==-1 && hasChosen[jj][numberOfShells]==-1)
              {
              queryCandidate = (*m_FineScheme)[jj];
              bool overlapBall = true;
              if (n>numberOfShells)
                {
                double innerProduct = 0;
                for ( int kk = 0; kk < 3; ++kk ) 
                  innerProduct+= queryPointLast[kk]*queryCandidate[kk];
                // if ( ( std::fabs(innerProduct) < overlapBallsInnerProduct[numberOfShells]) 
                //   || (indexShellMaxOverlapLast==s && std::fabs(innerProduct) < overlapBallsInnerProduct[s]) )
                if ( std::fabs(innerProduct) < overlapBallsInnerProduct[s] )
                  overlapBall = false;
                }

              int sumOverlap = 0;
              if (sumOverlapVec[s][jj]>=0 && !overlapBall)
                sumOverlap = sumOverlapVec[s][jj];
              else
                {
                 // sumOverlapVec[s][jj]<0 || (sumOverlapVec[s][jj]>=0 && overlapBall) 
                indexDistPair=jj + s*numberOfSamples;
                NeighborMapIterator mapIter = neighborMap.find(indexDistPair);
                if (mapIter==neighborMap.end())
                  {
                  m_KDTree->Search( queryCandidate, euclideanDistances[s], neighbors_jj );
                  neighborMap[indexDistPair] = neighbors_jj;
                  }
                else
                  neighbors_jj = mapIter->second;
                numberSearch++;
                for ( int k = 0; k < neighbors_jj.size(); ++k ) 
                  {
                  if (hasChosen[neighbors_jj[k]/2][s]>=0 || hasChosen[neighbors_jj[k]/2][numberOfShells]>=0)
                    sumOverlap++;
                  }
                sumOverlapVec[s][jj]=sumOverlap;
                }

              // utlPrintVar3(true, s, jj, sumOverlap);
              if (sumMaxOverlap<sumOverlap)
                {
                indexShellMaxOverlap = s;
                indexMaxOverlap = jj;
                sumMaxOverlap = sumOverlap;
                }
              numberOfCandidates++;
              }
            }

          }
        }

      if (numberOfCandidates==0)
        break;

      currentIndex = indexMaxOverlap;
      indexShellMaxOverlapLast = indexMaxOverlap;

      hasChosen[currentIndex][indexShellMaxOverlap]=0;
      hasChosen[currentIndex][numberOfShells]=indexShellMaxOverlap;
      chosenIndices[indexShellMaxOverlap].push_back(currentIndex);
      queryPoint = (*m_FineScheme)[currentIndex];
      queryPointLast = queryPoint;

      indexDistPair=currentIndex+indexShellMaxOverlap*numberOfSamples;
      NeighborMapIterator mapIter = neighborMap.find(indexDistPair);
      if (mapIter==neighborMap.end())
        {
        m_KDTree->Search( queryPoint, euclideanDistances[indexShellMaxOverlap], neighbors );
        neighborMap[indexDistPair] = neighbors;
        }
      else
        neighbors = mapIter->second;
      for ( int i = 0; i < neighbors.size(); ++i ) 
        {
        // update coverages
        if (needCoverageInshells)
          {
          if (hasChosen[neighbors[i]/2][indexShellMaxOverlap]==-1)
            coverageInShells[indexShellMaxOverlap]++;
          }
        hasChosen[neighbors[i]/2][indexShellMaxOverlap] = 0;
        }

      indexDistPair=currentIndex+numberOfShells*numberOfSamples;
      mapIter = neighborMap.find(indexDistPair);
      if (mapIter==neighborMap.end())
        {
        m_KDTree->Search( queryPoint, euclideanDistances.back(), neighbors_1 );
        neighborMap[indexDistPair] = neighbors_1;
        }
      else
        neighbors_1 = mapIter->second;
      for ( int i = 0; i < neighbors_1.size(); ++i ) 
        hasChosen[neighbors_1[i]/2][numberOfShells] = indexShellMaxOverlap;

      m_KDTree->Search( queryPoint, euclideanDistances[indexShellMaxOverlap]+1*m_MinDistanceInFineScheme, candidates );
      utlException(candidates.size()<=neighbors.size(), "the size of candidates should be larger than the size of neighbors.");

      m_KDTree->Search( queryPoint, euclideanDistances.back()+1*m_MinDistanceInFineScheme, candidates_1 );
      utlException(candidates_1.size()<=neighbors_1.size(), "the size of candidates should be larger than the size of neighbors.");

      for ( int i = 0; i < candidates.size(); ++i ) 
        {
        int jj = candidates[i]/2;
        if ( hasChosen[jj][indexShellMaxOverlap]==-1 && hasChosen[jj][numberOfShells]==-1)
          nextCandidateIndices[indexShellMaxOverlap].insert(jj);
        }

      for ( int i = 0; i < candidates_1.size(); ++i ) 
        {
        int jj = candidates_1[i]/2;
        if ( hasChosen[jj][numberOfShells]==-1)
          {
          for ( int sss = 0; sss < numberOfShells; ++sss ) 
            {
            if (sss!=indexShellMaxOverlap && hasChosen[jj][sss]==-1)
              nextCandidateIndices[sss].insert(jj);
            }
          nextCandidateIndices[numberOfShells].insert(jj);
          }
        }

      }

    // for ( int i = 0; i < chosenIndices.size(); ++i ) 
    //   utl::PrintVector(chosenIndices[i], "chosenIndices[i]");
    for ( int i = 0; i < chosenIndices.size(); ++i ) 
      {
      if (chosenIndices[i].size()!=this->m_NumbersInShell[i])
        return false;
      }
    *indices = chosenIndices;
    return true;

    }

  return true;
}

template< class TSamplingType >
void
SamplingSchemeQSpaceIMOCEstimationFilter< TSamplingType >
::GenerateData()
{
  utlGlobalException(this->m_CriteriaType!=Self::DISTANCE, "m_CriteriaType should be DISTANCE");

  Initialization();

  // m_KDTree->Print(std::cout<<"m_KDTree");

  int numberOfShells = this->m_NumbersInShell.size();
  int totalNumberOfSamples = std::accumulate(this->m_NumbersInShell.begin(), this->m_NumbersInShell.end(), 0);
  int angleSize = numberOfShells==1?1:(numberOfShells+1);
  std::vector<double> angles(angleSize,-1.0), anglesUpperBound(angleSize,-1.0), anglesLowerBound(angleSize,0.0);
  std::vector<std::vector<double> > anglesTest(angleSize);

  for ( int i = 0; i < numberOfShells; ++i ) 
    anglesUpperBound[i] = SamplingType::CalculateMinDistanceUpperBound(2*this->m_NumbersInShell[i]);
  if (numberOfShells>1)
    anglesUpperBound[numberOfShells] = SamplingType::CalculateMinDistanceUpperBound(2*totalNumberOfSamples); 

  if (numberOfShells==1)
    {
    double angle_0=0, angle_1=anglesUpperBound[0];
    while (angle_0<angle_1)
      {
      angles[0]= (angle_0+angle_1)*0.5;
      bool successAngles = IsSatisfiedSeparationAngles(angles);
      // utlPrintVar2(true, angles[0], successAngles);

      if (successAngles)
        {
        angle_0 = angles[0];
        anglesTest[0].push_back(angles[0]);
        }
      else
        angle_1 = angles[0];
      if (anglesTest[0].size()>=2 && std::fabs(anglesTest[0].back()-anglesTest[0][anglesTest[0].size()-2])/anglesTest[0][anglesTest[0].size()-2]<m_AngleMinChange )
        {
        this->m_OutputOrientations = m_FineScheme;
        // this->m_OutputOrientations->Print(std::cout<<"output");
        break;
        }
      }
    }
  else
    {

    std::vector<double> angle_0(anglesLowerBound), angle_1(anglesUpperBound);
    while (true)
      {
      for ( int i = 0; i < angleSize; ++i ) 
        angles[i] = (angle_0[i]+angle_1[i])*0.5;
      bool successAngles = IsSatisfiedSeparationAngles(angles);
      // utl::PrintVector(angles, "angles");
      // utlPrintVar1(true, successAngles);
      
      if (successAngles)
        {
        for ( int i = 0; i < angleSize; ++i ) 
          {
          angle_0[i] = angles[i];
          anglesTest[i].push_back(angles[i]);
          }
        }
      else
        {
        for ( int i = 0; i < angleSize; ++i ) 
          angle_1[i] = angles[i];
        }

      bool minAngleStop = true;
      for ( int i = 0; i < angleSize; ++i ) 
        {
        if (anglesTest[i].size()<2 || std::fabs(anglesTest[i].back()-anglesTest[i][anglesTest[i].size()-2])/anglesTest[i][anglesTest[i].size()-2]>=m_AngleMinChange )
          {
          minAngleStop = false;
          break;
          }
        }
      if (minAngleStop)
        {
        this->m_OutputOrientations = m_FineScheme;
        // this->m_OutputOrientations->Print(std::cout<<"output");
        break;
        }
      }

    }
}

}


#endif 


