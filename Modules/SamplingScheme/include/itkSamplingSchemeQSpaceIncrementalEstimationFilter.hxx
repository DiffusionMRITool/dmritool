/**
 *       @file  itkSamplingSchemeQSpaceIncrementalEstimationFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-14-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpaceIncrementalEstimationFilter_hxx
#define __itkSamplingSchemeQSpaceIncrementalEstimationFilter_hxx

#include "itkSamplingSchemeQSpaceIncrementalEstimationFilter.h"
#include "utl.h"

namespace itk
{

template< class TSamplingType >
SamplingSchemeQSpaceIncrementalEstimationFilter< TSamplingType >
::SamplingSchemeQSpaceIncrementalEstimationFilter() 
  : Superclass(),
  m_FineOrientations(new MatrixType()) 
{
  m_TessellationOrder = 7; 
}

template< class TSamplingType >
void
SamplingSchemeQSpaceIncrementalEstimationFilter< TSamplingType >
::Initialization()
{
  // SamplingType* input = const_cast< SamplingType * >( this->GetInput() );
  if (this->m_FineOrientations->Size()==0)
    {
    this->m_FineOrientations = utl::ReadGrad<double>(this->m_TessellationOrder, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); // catesian
    }

  if (this->m_InitialOrientations)
    {
    MatrixPointer initialMatrix = this->m_InitialOrientations->GetOrientationsCartesian();
    *this->m_FineOrientations = utl::ConnectUtlMatrix(*initialMatrix, *this->m_FineOrientations, true);
    
    SamplingPointer output = this->GetOutputOrientations();
    output->SetOrientationsCartesian(this->m_FineOrientations);
    output->SetIndicesInShells(this->m_InitialOrientations->GetIndicesInShells());
    Index2DVectorPointer indices = output->GetIndicesInShells();
    
    unsigned int numberInitialSamples = this->m_InitialOrientations->GetNumberOfSamples();
    unsigned int num=0;
    for ( unsigned int i = 0; i < indices->size(); i += 1 ) 
      {
      num += indices->operator[](i).size();
      }
    utlGlobalException(num!=numberInitialSamples, "wrong intialization in m_InitialOrientations");
    }
}

// template< class TSamplingType >
// void
// SamplingSchemeQSpaceIncrementalEstimationFilter< TSamplingType >
// ::IndicesOfMaximalDistance(const MatrixType& mat, unsigned int& r1, unsigned int& r2)
// {
//   unsigned int row = mat.Rows();
//   double x,y,z,x1,y1,z1, dotMin=100, dot;
//   for ( unsigned int i = 1; i < row; i += 1 ) 
//     {
//     x=mat(i,0),  y=mat(i,1),  z=mat(i,2);
//     for ( unsigned int j = 0; j < i; j += 1 ) 
//       {
//       x1=mat(i,0),  y1=mat(i,1), z1=mat(i,2);
//       dot = x*x1 + y*y1 + z*z1;
//       if ( dot<0 )
//         dot = -dot;

//       if (dot<dotMin)
//         {
//         dotMin = dot;
//         r1=i, r2=j;
//         }
//       }
//     }
// }

template< class TSamplingType >
void
SamplingSchemeQSpaceIncrementalEstimationFilter< TSamplingType >
::GenerateData()
{
  Initialization();

  unsigned int numberOfShells = this->m_NumbersInShell.size();
  unsigned int numberInitialSamples = 0;
  if (this->m_InitialOrientations) 
    numberInitialSamples = this->m_InitialOrientations->GetNumberOfSamples();

  SamplingPointer output = this->GetOutputOrientations();
  output->SetOrientationsCartesian(this->m_FineOrientations);
  Index2DVectorPointer indices = output->GetIndicesInShells();

  unsigned int numberOfSamples = output->GetNumberOfSamples();
  MatrixPointer outputMatrix = output->GetOrientationsCartesian();

  // itk::TimeProbe clock;
  // clock.Start();
  // MatrixPointer dotMatrix (new MatrixType()); 
  // MatrixPointer elecMatrix (new MatrixType()); 
  // if (numberOfShells>1)
  //   {
  //   if (this->m_CriteriaType==Self::DISTANCE)
  //     dotMatrix = output->CalculateInnerProductMatrix(true);
  //   else if (this->m_CriteriaType==Self::ELECTROSTATIC)
  //     elecMatrix = output->CalculateElectrostaticEnergyMatrix(this->m_ElectrostaticOrder);
  //   }
  // clock.Stop();
  // std::cout << "after matrix calculation: " << clock.GetTotal() << std::endl;

  utlGlobalException(this->m_NumbersInShell.size()==0, "need to set this->m_NumbersInShell");
  unsigned int totalNumbers=0;
  for ( unsigned int i = 0; i < this->m_NumbersInShell.size(); i += 1 ) 
    {
    totalNumbers += this->m_NumbersInShell[i];
    }

  std::vector<int> chosenIndex;
  std::vector<char> hasChosen(numberOfSamples, 0);
  std::cout << "numberOfShells = " << numberOfShells << std::endl << std::flush;
  std::cout << "numberOfSamples = " << numberOfSamples << std::endl << std::flush;

  if (!this->m_InitialOrientations)
    {
    // first two points
    int d0, d1;
    // clock.Start();
    // if (numberOfShells>1)
    //   {
    //   if (this->m_CriteriaType==Self::DISTANCE)
    //     utl::ArgminSymmetricMatrix<MatrixType>(*dotMatrix, numberOfSamples, d0, d1, false);
    //   else if (this->m_CriteriaType==Self::ELECTROSTATIC)
    //     utl::ArgminSymmetricMatrix<MatrixType>(*elecMatrix, numberOfSamples, d0, d1, false);
    //   }

    IndexVectorType indexTmp;
    // utlPrintVar2(true, d0, d1);
    // clock.Stop();
    // std::cout << "after matrix minimal: " << clock.GetTotal() << std::endl;
    if (numberOfShells==1)
      {
      d0=0;
      indexTmp.push_back(d0);
      // indexTmp.push_back(d1);
      indices->push_back(indexTmp);
      }
    else
      {
      d0=0;
      indexTmp.push_back(d0);
      indices->push_back(indexTmp);
      // indexTmp[0] = d1;
      // indices->push_back(indexTmp);
      }
    chosenIndex.push_back(d0);
    // chosenIndex.push_back(d1);
    hasChosen[d0] = 1; 
    // hasChosen[d1] = 1;
    }

  for ( unsigned int i = 0; i < indices->size(); i += 1 ) 
    for ( unsigned int j = 0; j < (*indices)[i].size(); j += 1 ) 
      {
      int ii = (*indices)[i][j];
      chosenIndex.push_back(ii);
      hasChosen[ii] = 1;
      }

  if (numberOfShells==1)
    {

    if (this->m_CriteriaType==Self::DISTANCE)
      {
      double x,y,z,x1,y1,z1;
      std::vector<double> maxDots(numberOfSamples, -1); 
      // clock.Start();
      for ( unsigned int i = (this->m_InitialOrientations?numberInitialSamples:1); i < totalNumbers; i += 1 ) 
        {
        int index=-1;
        double minValue=std::numeric_limits<double>::max();
        for ( unsigned int j = 0; j < numberOfSamples; j += 1 ) 
          {
          if (hasChosen[j]==0)
            {
            double maxValue_j = -std::numeric_limits<double>::max();
            x = output->operator[](j)[0];
            y = output->operator[](j)[1];
            z = output->operator[](j)[2];
            if (maxDots[j]< 0 )
              {
              for ( unsigned int k = 0; k < chosenIndex.size(); k += 1 ) 
                {
                unsigned int kk = chosenIndex[k];
                // double value = dotMatrix->operator()(j,kk);
                x1 = output->operator[](kk)[0];
                y1 = output->operator[](kk)[1];
                z1 = output->operator[](kk)[2];
                double value = x*x1 + y*y1 + z*z1;
                if (value<0)
                  value = -value;
                if (value>maxValue_j)
                  maxValue_j = value;
                }
              maxDots[j] = maxValue_j;
              }
            else
              {
              unsigned int kk = chosenIndex.back();
              // double value = dotMatrix->operator()(j,kk);
              x1 = output->operator[](kk)[0];
              y1 = output->operator[](kk)[1];
              z1 = output->operator[](kk)[2];
              double value = x*x1 + y*y1 + z*z1;
              if (value<0)
                value = -value;
              if (value > maxDots[j] )
                maxDots[j] = value;
              maxValue_j = maxDots[j];
              }

            if (maxValue_j < minValue )
              {
              minValue = maxValue_j;
              index = j;
              }
            }
          }

        chosenIndex.push_back(index);
        hasChosen[index] = 1;
        indices->operator[](0).push_back(index);
        }
    // clock.Stop();
    // std::cout << "after loop: " << clock.GetTotal() << std::endl;
      }
    else if (this->m_CriteriaType==Self::ELECTROSTATIC)
      {
      double x,y,z,x1,y1,z1;
      bool orderEqual2 = std::abs(this->m_ElectrostaticOrder-2)<1e-8 ? true : false;
      std::vector<double> energy(numberOfSamples,0.0);
      // clock.Start();
      for ( unsigned int i = (this->m_InitialOrientations?numberInitialSamples:1); i < totalNumbers; i += 1 ) 
        {
        int index=-1;
        double minValue=std::numeric_limits<double>::max();
        for ( unsigned int j = 0; j < numberOfSamples; j += 1 ) 
          {
          if (hasChosen[j]==0)
            {
            x = output->operator[](j)[0];
            y = output->operator[](j)[1];
            z = output->operator[](j)[2];
            if (energy[j]<1e-10)
              {
              for ( unsigned int k = 0; k < chosenIndex.size(); k += 1 ) 
                {
                unsigned int kk = chosenIndex[k];
                // double value = elecMatrix->operator()(j,kk);
                x1 = output->operator[](kk)[0];
                y1 = output->operator[](kk)[1];
                z1 = output->operator[](kk)[2];
                double dot = x*x1 + y*y1 + z*z1;
                double result1 = 2+2*dot;
                double result2 = 2-2*dot;
                double value = (!orderEqual2) ? (1.0/std::pow(result1, this->m_ElectrostaticOrder*0.5) + 1.0/std::pow(result2, this->m_ElectrostaticOrder*0.5)) : (1.0/result1 + 1.0/result2);  
                energy[j] += value;
                }
              }
            else
              {
              unsigned int kk = chosenIndex.back();
              // double value = dotMatrix->operator()(j,kk);
              x1 = output->operator[](kk)[0];
              y1 = output->operator[](kk)[1];
              z1 = output->operator[](kk)[2];
              double dot = x*x1 + y*y1 + z*z1;
              double result1 = 2+2*dot;
              double result2 = 2-2*dot;
              double value = (!orderEqual2) ? (1.0/std::pow(result1, this->m_ElectrostaticOrder*0.5) + 1.0/std::pow(result2, this->m_ElectrostaticOrder*0.5)) : (1.0/result1 + 1.0/result2);  
              energy[j] += value;
              }
            if (energy[j] < minValue )
              {
              minValue = energy[j];
              index = j;
              }
            }
          }
        chosenIndex.push_back(index);
        hasChosen[index] = 1;
        indices->operator[](0).push_back(index);
        }
      // clock.Stop();
      // std::cout << "after loop: " << clock.GetTotal() << std::endl;
      }
    else
      utlGlobalException(true, "wrong m_CriteriaType");

    }
  else
    {

    if (this->m_CriteriaType==Self::DISTANCE)
      {

      double x,y,z,x1,y1,z1;
      std::vector< std::vector<double> > maxDotsVec;
      for ( int i = 0; i < numberOfShells+1; i += 1 ) 
        {
        std::vector<double> dotTemp(numberOfSamples,-1);
        maxDotsVec.push_back(dotTemp);
        }

      /**************************************************************************************/
      // first several points for the shell without initialization. 
      unsigned startI = (this->m_InitialOrientations?numberInitialSamples:1);
      for ( ; startI < numberOfShells; startI += 1 ) 
        {
        int index=-1;
        double minValue=std::numeric_limits<double>::max();
        for ( unsigned int j = 0; j < numberOfSamples; j += 1 ) 
          {
          if (hasChosen[j]==0)
            {
            double maxValue_j = -std::numeric_limits<double>::max();
            x = output->operator[](j)[0];
            y = output->operator[](j)[1];
            z = output->operator[](j)[2];
            if (maxDotsVec[numberOfShells][j]< 0 )
              {
              for ( unsigned int k = 0; k < chosenIndex.size(); k += 1 ) 
                {
                unsigned int kk = chosenIndex[k];
                // double value = dotMatrix->operator()(j,kk);
                x1 = output->operator[](kk)[0];
                y1 = output->operator[](kk)[1];
                z1 = output->operator[](kk)[2];
                double value = x*x1 + y*y1 + z*z1;
                if (value<0)
                  value = -value;
                if (value>maxValue_j)
                  maxValue_j = value;
                }
              maxDotsVec[numberOfShells][j] = maxValue_j;
              }
            else
              {
              unsigned int kk = chosenIndex.back();
              // double value = dotMatrix->operator()(j,kk);
              x1 = output->operator[](kk)[0];
              y1 = output->operator[](kk)[1];
              z1 = output->operator[](kk)[2];
              double value = x*x1 + y*y1 + z*z1;
              if (value<0)
                value = -value;
              if (value > maxDotsVec[numberOfShells][j] )
                maxDotsVec[numberOfShells][j] = value;
              maxValue_j = maxDotsVec[numberOfShells][j];
              }

            if (maxValue_j < minValue )
              {
              minValue = maxValue_j;
              index = j;
              }
            }
          }
        chosenIndex.push_back(index);
        hasChosen[index] = 1;
        IndexVectorType tmp;
        tmp.push_back(index);
        indices->push_back(tmp);
        }

      utlException(indices->size()!=numberOfShells, "Logic ERROR!");
      // utlPrintVar1(true, startI);
      // utl::Save2DVector(*indices, std::cout<<"indices = ");

      /**************************************************************************************/
      // other points
      int chosedShellIndex_final = -1;
      int chosedShellIndex_final_last = -1;
      for ( unsigned int i = startI; i < totalNumbers; i += 1 ) 
        {
        int index=-1;
        double maxDist = -std::numeric_limits<double>::max();
        int chosedShellIndex = -1;

        for ( unsigned int j = 0; j < numberOfSamples; j += 1 ) 
          {
          if (hasChosen[j]==0)
            {

            // maximal dots in single shells
            std::vector<double> maxValue_j(numberOfShells+1, -std::numeric_limits<double>::max()); 
            x = output->operator[](j)[0];
            y = output->operator[](j)[1];
            z = output->operator[](j)[2];

            // utlPrintVar1(true, chosedShellIndex_final);
            if (chosedShellIndex_final_last>=0)
              {
              // update stored maxDotsVec based on last added point
              unsigned int kk = chosenIndex.back();
              // double value = dotMatrix->operator()(j,kk);
              x1 = output->operator[](kk)[0];
              y1 = output->operator[](kk)[1];
              z1 = output->operator[](kk)[2];
              double value = x*x1 + y*y1 + z*z1;
              if (value<0)
                value = -value;
              if (value > maxDotsVec[chosedShellIndex_final_last][j] )
                maxDotsVec[chosedShellIndex_final_last][j] = value;
              if (value > maxDotsVec[numberOfShells][j] )
                maxDotsVec[numberOfShells][j] = value;
              }

            for ( unsigned int s = 0; s < numberOfShells; s += 1 ) 
              {
              int numbers = (*indices)[s].size();
              if (numbers < this->m_NumbersInShell[s])
                {
                // only add point if there is no enough number 
                if (maxDotsVec[s][j]<0)
                  {
                  // only calculate all dots in the first begining
                  for ( unsigned int k = 0; k < numbers; k += 1 ) 
                    {
                    unsigned int kk = (*indices)[s][k];
                    // double value = dotMatrix->operator()(j,kk);
                    x1 = output->operator[](kk)[0];
                    y1 = output->operator[](kk)[1];
                    z1 = output->operator[](kk)[2];
                    double value = x*x1 + y*y1 + z*z1;
                    if (value<0)
                      value = -value;
                    if (value>maxValue_j[s])
                      maxValue_j[s] = value;
                    }
                  maxDotsVec[s][j] = maxValue_j[s];
                  if (maxValue_j[s]>maxDotsVec[numberOfShells][j])
                    maxDotsVec[numberOfShells][j] = maxValue_j[s];
                  }
                else
                  {
                  maxValue_j[s] = maxDotsVec[s][j];
                  }
                }
              }

            // maximal dot in whole shell, which is the maxinmal value amoung maximal dots in each shell
            // for ( unsigned int s = 0; s < numberOfShells; s += 1 ) 
            //   {
            //   if (maxDotsVec[s][j] > maxValue_j[numberOfShells])
            //     maxValue_j[numberOfShells] = maxDotsVec[s][j];
            //   }
            
            // maximal dot in whole shell
            maxValue_j[numberOfShells] = maxDotsVec[numberOfShells][j];


            // choose the shell with minimal maximal dot
            chosedShellIndex = -1;
            double chosedShellDot = std::numeric_limits<double>::max();
            for ( unsigned int s = 0; s < numberOfShells; s += 1 ) 
              {
              if ( (*indices)[s].size()==this->m_NumbersInShell[s])
                continue;
              if (maxValue_j[s] < chosedShellDot)
                {
                chosedShellIndex = s;
                chosedShellDot = maxValue_j[s];
                }
              }
            // utlPrintVar2(true, j, chosedShellIndex);
            // utl::PrintVector(maxValue_j, "maxValue_j");

            // distance in cost function which combines whole shell and the chosed shell
            double costDist = (1-this->m_WeightForSingleShell)*std::acos(maxValue_j[numberOfShells]) 
                  + this->m_WeightForSingleShell*std::acos(maxValue_j[chosedShellIndex]);
            if (costDist > maxDist )
              {
              maxDist = costDist;
              index = j;
              chosedShellIndex_final = chosedShellIndex;
              }

            }
          }

        chosenIndex.push_back(index);
        hasChosen[index] = 1;
        indices->operator[](chosedShellIndex_final).push_back(index);
        chosedShellIndex_final_last = chosedShellIndex_final;
        }
      
      }
    else if (this->m_CriteriaType==Self::ELECTROSTATIC)
      {

      utlGlobalException(true, "TODO");
      }
    else
      utlGlobalException(true, "wrong m_CriteriaType");

    }
}

}


#endif 



