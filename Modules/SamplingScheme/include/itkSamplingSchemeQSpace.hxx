/**
 *       @file  itkSamplingSchemeQSpace.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-25-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpace_hxx
#define __itkSamplingSchemeQSpace_hxx

#include "itkSamplingSchemeQSpace.h"

namespace itk
{

template <class TPixelType>
SamplingSchemeQSpace<TPixelType>
::SamplingSchemeQSpace() : Superclass(), 
  m_BVector(new STDVectorType()) 
{
  m_BThresholdSingleShell = 1e-4;
}

template <class TPixelType>
typename LightObject::Pointer
SamplingSchemeQSpace<TPixelType>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_BThresholdSingleShell = m_BThresholdSingleShell;
  rval->m_BVector = m_BVector;

  return loPtr;
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::PrintSelf(std::ostream& os, Indent indent) const
{  
  Superclass::PrintSelf(os, indent);
  utl::PrintVector(*m_BVector, "m_BVector", " ", os<<indent);
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::SetBVector(const STDVectorPointer bVec)
{
  if (m_BVector != bVec)
    {
    m_BVector = bVec;
    ConvertBVectorToQVector();
    if (m_BThresholdSingleShell>0)
      GroupBValues();
    this->Modified();
    }
}

template <class TPixelType>
typename SamplingSchemeQSpace<TPixelType>::STDVectorPointer
SamplingSchemeQSpace<TPixelType>
::GetBVectorInShell(unsigned int shellIndex)
{
  unsigned int num = this->GetNumberOfShells();
  STDVectorPointer bVec(new STDVectorType());
  if (shellIndex<num && m_BVector->size()>num)
    {
    IndexVectorType vecTmp = (*this->m_IndicesInShells)[shellIndex];
    for ( unsigned int i = 0; i < vecTmp.size(); i += 1 ) 
      {
      bVec->push_back( (*m_BVector)[vecTmp[i] ] );
      }
    }
  return bVec;
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::SetSamplingScheme3D(typename Superclass::Pointer scheme)
{
  this->clear();
  for ( int i = 0; i < scheme->size(); i += 1 ) 
    this->push_back((*scheme)[i]);

  this->m_DeltaSmall = scheme->GetDeltaSmall();
  this->m_DeltaBig = scheme->GetDeltaBig();
  this->m_Tau = scheme->GetTau();
  this->m_RadiusThresholdSingleShell = scheme->GetRadiusThresholdSingleShell();

  *this->m_RadiusVector = *scheme->GetRadiusVector();
  *this->m_IndicesInShells = *scheme->GetIndicesInShells();
}


template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::Clear()
{
  Superclass::Clear();
  m_BVector=STDVectorPointer(new STDVectorType());
  this->Modified();
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::ConvertBVectorToQVector()
{
  // utlShowPosition(true);
  // std::cout << "this->m_Tau = " << this->m_Tau << std::endl << std::flush;
  this->m_RadiusVector=STDVectorPointer(new STDVectorType(m_BVector->size()));
  double firstTerm = 4.0*M_PI*M_PI*this->m_Tau;
  for ( int i = 0; i < m_BVector->size(); i += 1 ) 
    (*this->m_RadiusVector)[i] = std::sqrt((*m_BVector)[i]/(firstTerm)); 
  this->Modified();
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::ConvertQVectorToBVector()
{
  this->m_BVector=STDVectorPointer(new STDVectorType(this->m_RadiusVector->size()));
  double firstTerm = 4.0*M_PI*M_PI*this->m_Tau;
  for ( int i = 0; i < m_BVector->size(); i += 1 ) 
    (*m_BVector)[i] = firstTerm * (*this->m_RadiusVector)[i] * (*this->m_RadiusVector)[i];
  this->Modified();
}

template <class TPixelType>
std::vector<typename SamplingSchemeQSpace<TPixelType>::STDVectorType>
SamplingSchemeQSpace<TPixelType>
::GroupBValues()
{
  // utlGlobalException(this->m_RadiusVector->size()==0 && m_BVector->size()==0, "need to set bVector or qVector");
  // if (this->m_RadiusVector->size()==0 && m_BVector->size()>0)
  //   ConvertBVectorToQVector();

  std::vector<STDVectorType> bVectors; 
  if (this->m_IndicesInShells->size()==0)
    {
    this->m_IndicesInShells = Index2DVectorPointer(new Index2DVectorType()); 
    STDVectorType bMax, bMin;
    for ( int i = 0; i < m_BVector->size(); i += 1 ) 
      {
      double b = (*m_BVector)[i];
      int j=0;
      for ( j = 0; j < bVectors.size(); j += 1 ) 
        {
        if (b>=bMin[j]-m_BThresholdSingleShell && b<=bMax[j]+m_BThresholdSingleShell)
          {
          bVectors[j].push_back(b);
          (*this->m_IndicesInShells)[j].push_back(i);
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
        IndexVectorType bIndexTemp;
        bIndexTemp.push_back(i);
        this->m_IndicesInShells->push_back(bIndexTemp);
        bMin.push_back(b);
        bMax.push_back(b);
        }
      }

    this->Modified();
    for ( int j = 0; j < bVectors.size(); j += 1 ) 
      {
      utlSAGlobalException(bMax[j]-bMin[j]>2*m_BThresholdSingleShell)
        (j)(bMax[j])(bMin[j])(m_BThresholdSingleShell).msg("the range of b values is larger than 100, which can not be in the same shell");
      }
    }
  else
    {
    for ( int i = 0; i < this->m_IndicesInShells->size(); i += 1 ) 
      {
      STDVectorType bVecTemp;
      IndexVectorType indexTemp = (*this->m_IndicesInShells)[i];
      for ( int j = 0; j < indexTemp.size(); j += 1 ) 
        {
        bVecTemp.push_back((*m_BVector)[ indexTemp[j] ]);
        }
      bVectors.push_back(bVecTemp);
      }
    }
  return bVectors;
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::CorrectBValues()
{
  utlGlobalException(m_BVector->size()==0, "no b values");
  std::vector<STDVectorType> bVectors = GroupBValues(); 

  STDVectorPointer bVec (new STDVectorType(m_BVector->size()));
  for ( int j = 0; j < bVectors.size(); j += 1 ) 
    {
    double bMean=0;
    for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
      bMean += bVectors[j][k];
    bMean /= bVectors[j].size();
    for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
      {
      (*bVec)[(*this->m_IndicesInShells)[j][k] ] = bMean;
      }
    }
  m_BVector = bVec;
  ConvertBVectorToQVector();
  this->Modified();
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::CorrectRadiusValues()
{
  utlGlobalException(this->m_RadiusVector->size()==0 && m_BVector->size()==0, "need to set bVector or qVector");
  if (this->m_RadiusVector->size()==0 && m_BVector->size()>0)
    ConvertBVectorToQVector();
  Superclass::CorrectRadiusValues();
  ConvertQVectorToBVector();
}

template <class TPixelType>
void
SamplingSchemeQSpace<TPixelType>
::RemoveSamplesNotIndexed()
{
  if (this->size()==0)
    return;
  utlGlobalException(this->m_IndicesInShells->size()==0, "need to set m_IndicesInShells first");

  if (m_BVector->size()>0)
    {
    utlGlobalException(m_BVector->size()!=this->size(), "different size of bVec and gradients");
    STDVectorType bVec = *m_BVector;
    m_BVector->clear();
    for ( int i = 0; i < this->GetNumberOfShells(); ++i ) 
      {
      IndexVectorType indices = (*this->m_IndicesInShells)[i];
      for ( int j = 0; j < indices.size(); ++j ) 
        {
        m_BVector->push_back(bVec[indices[j]]);
        }
      }
    }

  Superclass::RemoveSamplesNotIndexed();
}

}


#endif 




