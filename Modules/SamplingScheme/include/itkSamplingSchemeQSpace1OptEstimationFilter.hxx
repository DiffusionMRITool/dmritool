/**
 *       @file  itkSamplingSchemeQSpace1OptEstimationFilter.hxx
 *      @brief  
 *     Created  "01-16-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#ifndef __itkSamplingSchemeQSpace1OptEstimationFilter_hxx
#define __itkSamplingSchemeQSpace1OptEstimationFilter_hxx

#include "itkSamplingSchemeQSpace1OptEstimationFilter.h"
#include "utl.h"

namespace itk
{

template< class TSamplingType >
SamplingSchemeQSpace1OptEstimationFilter< TSamplingType >
::SamplingSchemeQSpace1OptEstimationFilter() 
  : Superclass(),
  m_FineOrientations(new MatrixType()) 
{
  m_TessellationOrder = 7; 
  m_FineScheme = SamplingType::New();
}

template< class TSamplingType >
void
SamplingSchemeQSpace1OptEstimationFilter< TSamplingType >
::Initialization()
{
  if (this->m_FineScheme->GetNumberOfSamples()==0)
    {
    if (this->m_FineOrientations->Size()==0)
      this->m_FineOrientations = utl::ReadGrad<double>(this->m_TessellationOrder, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN); // catesian
    m_FineScheme->SetOrientationsCartesian(m_FineOrientations);
    }

  utlGlobalException(!this->IsSetInitialization(), "Need to set m_InitialOrientations for initialization!");

  this->m_OutputOrientations = this->m_InitialOrientations->Clone();
  
  Index2DVectorPointer indices = this->m_OutputOrientations->GetIndicesInShells();
  this->m_NumbersInShell.resize(indices->size());
  for ( int i = 0; i < indices->size(); ++i ) 
    {
    this->m_NumbersInShell[i] = ((*indices)[i].size());
    }
  
}


template< class TSamplingType >
void
SamplingSchemeQSpace1OptEstimationFilter< TSamplingType >
::GenerateData()
{
  Initialization();

  SamplingPointer output = this->GetOutputOrientations();
  Index2DVectorPointer indices = output->GetIndicesInShells();
  unsigned int numberOfSamples = output->GetNumberOfSamples();
  unsigned int N = m_FineScheme->GetNumberOfSamples();
  MatrixPointer outputMatrix = output->GetOrientationsCartesian();
  unsigned int numberOfShells = this->m_NumbersInShell.size();

  std::vector<char> hasChosen(N, 0);
  std::vector<int> chosenIndex;


  if (numberOfShells==1)
    {

    MatrixType dotMat;
    Index2DVectorType indexNear2(N);
    STDVectorType dotVec(numberOfSamples);

    for ( int k = 0; k < numberOfSamples; ++k ) 
      {
      dotVec[k] = output->CalculateMaxDot(k, true);
      }

    dotMat = (*output->GetOrientationsCartesian()) * m_FineScheme->GetOrientationsCartesian()->GetTranspose();
    dotMat = utl::Abs(dotMat);

    bool has1Opt=true;
    double dotMax_k=0, dotMax_k2=0;

    do 
      {
      double minDot = 100;
      int index_j=-1, index_k=-1;
        
      IndexVectorType index2(2,-1);
      for ( int i = 0; i < N; ++i ) 
        {
        if (hasChosen[i]>0)
          continue;
        utl::Vector<double> vec = dotMat.GetColumn(i);
        utl::argmax2(vec, vec.Size(), index2[0], index2[1]);
        indexNear2[i] = index2;
        }
      // for ( int i= 0; i < N; ++i ) 
      //   {
      //   utl::PrintVector(indexNear2[i], "", " ", std::cout<<"indexNear2[" <<i<<"]");
      //   utlPrintVar(true, std::acos(dotMat(indexNear2[i][0],i)), std::acos(dotMat(indexNear2[i][1],i)));
      //   }

      for ( unsigned int k = 0; k < numberOfSamples; k += 1 ) 
        {
        double dot_k=dotVec[k];

        for ( unsigned int j = 0; j < N; j += 1 ) 
          {
          if (hasChosen[j]>0)
            continue;
          // utlPrintVar(true, j, k);
          // PointType p = (*this->m_FineScheme)[j]; 
          double dot_k2=-1;

          int index2_0 = indexNear2[j][0];
          int index2_1 = indexNear2[j][1];
          if (index2_0!=k)
            dot_k2 = dotMat(index2_0,j);
          else
            dot_k2 = dotMat(index2_1,j);

          if (dot_k2<dot_k)
            {
            if (minDot>dot_k2)
              {
              minDot=dot_k2;
              index_j = j;
              index_k = k;
              dotMax_k = dot_k;
              dotMax_k2 = dot_k2;
              }
            }
          // (*output)[k] = p1;
          }
        }

      if (index_j>=0 && index_k>=0)
        {
        utlPrintVar(utl::IsLogDebug(), index_j, index_k, std::acos(dotMax_k), std::acos(dotMax_k2));
        hasChosen[index_j] = 1;
        (*output)[index_k] = (*this->m_FineScheme)[index_j]; 
        // dotVec[index_k] = dotMax_k2;
        for ( int k = 0; k < numberOfSamples; ++k ) 
          {
          dotVec[k] = output->CalculateMaxDot(k, true);
          }


        PointType p1 = (*output)[index_k]; 
        for ( int i = 0; i < N; ++i ) 
          {
          if (hasChosen[i]>0)
            continue;
          PointType p2 = (*this->m_FineScheme)[i]; 
          double dotTmp = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
          if (dotTmp<0)
            dotTmp = -dotTmp;
          dotMat(index_k, i) = dotTmp;
          }

        // dotMat = (*output->GetOrientationsCartesian(true)) * m_FineScheme->GetOrientationsCartesian()->GetTranspose();
        // dotMat = utl::Abs(dotMat);
        // dotMat = utl::Acos(utl::Abs(dotMat));

        }
      else
        has1Opt=false;

      } while ( has1Opt );

    }
  else
    {
    std::vector<MatrixType> dotMat(numberOfShells+1);
    std::vector<Index2DVectorType> indexNear2(numberOfShells+1, Index2DVectorType(N));

    for ( int s = 0; s < numberOfShells; ++s ) 
      {
      utl::ProductUtlMMt(*output->GetOrientationsCartesianInShell(s), *m_FineScheme->GetOrientationsCartesian(), dotMat[s]);
      dotMat[s] = utl::Abs(dotMat[s]);
      }
    utl::ProductUtlMMt(*output->GetOrientationsCartesian(), *m_FineScheme->GetOrientationsCartesian(), dotMat[numberOfShells]);
    dotMat.back() = utl::Abs(dotMat.back());

    bool has1Opt=true;
    STDVectorType dot_k(numberOfShells+1, -1);
    STDVectorType dot_k2(numberOfShells+1, -1);
    STDVectorType dotMax_k(numberOfShells+1, -1);
    STDVectorType dotMax_k2(numberOfShells+1, -1);
    double w0 = 1.0-this->m_WeightForSingleShell;
    double ws = this->m_WeightForSingleShell/numberOfShells;

    do 
      {
      STDVectorType minDot(numberOfShells+1, 100);
      int index_s=-1, index_j=-1, index_k=-1;

      Index2DVectorPointer indices = output->GetIndicesInShells();
      
      IndexVectorType index2(2,-1);
      for ( int s = 0; s < numberOfShells+1; ++s ) 
        {
        for ( int i = 0; i < N; ++i ) 
          {
          if (hasChosen[i]>0)
            continue;
          utl::Vector<double> vec = dotMat[s].GetColumn(i);
          utl::argmax2(vec, vec.Size(), index2[0], index2[1]);
          indexNear2[s][i] = index2;
          }
        }

      for ( int s = 0; s < numberOfShells; ++s ) 
        {
        for ( unsigned int k = 0; k < this->m_NumbersInShell[s]; k += 1 ) 
          {
          int ind = (*indices)[s][k];
          dot_k[s] = output->CalculateMaxDotInShell(ind, s, true);
          dot_k.back()= output->CalculateMaxDot(ind, true);

          for ( unsigned int j = 0; j < N; j += 1 ) 
            {
            if (hasChosen[j]>0)
              continue;
            // utlPrintVar(true, j, k);
            // PointType p = (*this->m_FineScheme)[j]; 
            // double dot_k2_s=-1, dot_k2=-1;

            int index2_0 = indexNear2[s][j][0];
            int index2_1 = indexNear2[s][j][1];
            if (index2_0!=k)
              dot_k2[s] = dotMat[s](index2_0,j);
            else
              dot_k2[s] = dotMat[s](index2_1,j);
            
            index2_0 = indexNear2[numberOfShells][j][0];
            index2_1 = indexNear2[numberOfShells][j][1];
            if (index2_0!=ind)
              dot_k2.back() = dotMat[numberOfShells](index2_0,j);
            else
              dot_k2.back() = dotMat[numberOfShells](index2_1,j);

            // utlPrintVar(utl::IsLogDebug(), j, s, k, dot_k[s], dot_k2[s], dot_k.back(), dot_k2.back());
            // if ( w0*std::acos(dot_k2.back()) +  ws*std::acos(dot_k2[s]) > w0*std::acos(dot_k.back()) +  ws*std::acos(dot_k[s]))
            if ( (dot_k2.back() <= dot_k.back() && dot_k2[s] < dot_k[s]) || (dot_k2.back() < dot_k.back() && dot_k2[s] <= dot_k[s]))
              {
              if ((minDot.back()>=dot_k2.back() && minDot[s]>dot_k2[s]) || (minDot.back()>dot_k2.back() && minDot[s]>=dot_k2[s]))
                {
                minDot.back()=dot_k2.back();
                minDot[s]=dot_k2[s];
                index_s = s;
                index_j = j;
                index_k = k;
                if (utl::IsLogDebug())
                  {
                  dotMax_k = dot_k;
                  dotMax_k2 = dot_k2;
                  }
                }
              }
            // (*output)[k] = p1;
            }
          }
        }

      if (index_j>=0 && index_k>=0 && index_s>=0)
        {
        utlPrintVar(utl::IsLogDebug(), index_j, index_s, index_k, std::acos(dotMax_k[index_s]), std::acos(dotMax_k2[index_s]), std::acos(dotMax_k.back()), std::acos(dotMax_k2.back()));
        hasChosen[index_j] = 1;
        int ind = (*indices)[index_s][index_k];
        (*output)[ind] = (*this->m_FineScheme)[index_j]; 
        // dotVec[index_k] = dotMax_k2;

        // for ( int k = 0; k < numberOfSamples; ++k ) 
        //   {
        //   dotVec[k] = output->CalculateMaxDot(k, true);
        //   }


        PointType p1 = (*output)[ind]; 
        for ( int i = 0; i < N; ++i ) 
          {
          if (hasChosen[i]>0)
            continue;
          PointType p2 = (*this->m_FineScheme)[i]; 

          double dotTmp = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
          if (dotTmp<0)
            dotTmp = -dotTmp;
          dotMat[numberOfShells](ind, i) = dotTmp;
          dotMat[index_s](index_k, i) = dotTmp;
          }


        }
      else
        has1Opt=false;

      } while ( has1Opt );

    }

}

}


#endif 


