/**
 *       @file  itkSamplingSchemeQSpaceWriter.hxx
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



#ifndef __itkSamplingSchemeQSpaceWriter_hxx
#define __itkSamplingSchemeQSpaceWriter_hxx

#include "itkSamplingSchemeQSpaceWriter.h"
#include "utl.h"

namespace itk
{

template <class TSamplingType>
SamplingSchemeQSpaceWriter<TSamplingType>
::SamplingSchemeQSpaceWriter()
{
  m_BFile="";
  m_OrientationFile="";
  m_Sampling = NULL;
  
  m_SaveSingleShell = false;
  m_SaveAllShellsInOneFile = true;
}


template <class TSamplingType>
void 
SamplingSchemeQSpaceWriter<TSamplingType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "m_BFile: " << m_BFile << std::endl;
  os << indent << "m_OrientationFile: " << m_OrientationFile << std::endl;
  os << indent << "m_SaveSingleShell: " << m_SaveSingleShell << std::endl;
  os << indent << "m_SaveAllShellsInOneFile: " << m_SaveAllShellsInOneFile << std::endl;
  m_Sampling->Print(os<< indent << "m_Sampling\n");
}

template <class TSamplingType>
void 
SamplingSchemeQSpaceWriter<TSamplingType>
::Update()
{
  utlGlobalException(!m_Sampling, "need to set m_Sampling");
  utlGlobalException(!m_SaveSingleShell && !m_SaveAllShellsInOneFile, "need to set m_SaveSingleShell or m_SaveAllShellsInOneFile");
  utlGlobalException(m_OrientationFile=="", "need to set m_OrientationFile");

  if (m_SaveAllShellsInOneFile)
    {
    typename SamplingType::MatrixPointer orientationMatrix = m_Sampling->GetOrientationsCartesian();
    std::cout << "save all orientations to " << m_OrientationFile << std::endl << std::flush;
    orientationMatrix->Save(m_OrientationFile);
    if (m_BFile!="")
      {
      typename SamplingType::STDVectorPointer bVector = m_Sampling->GetBVector();
      std::cout << "save all bValues to " << m_BFile << std::endl << std::flush;
      utl::SaveVector(*bVector, m_BFile);
      }
    }

  if (m_SaveSingleShell)
    {
    unsigned int numberOfShells = m_Sampling->GetNumberOfShells();
    utlGlobalException(numberOfShells==0, "need to set IndicesInShells in m_Sampling");

    typename SamplingType::MatrixPointer orientationMatrix;
    typename SamplingType::STDVectorPointer bVector;
    std::string realBFile, realOrientationFile, bExt, orientationExt, bNoExt, orientationNoExt, indexStr;
    if (m_BFile!="")
      utl::GetFileExtension(m_BFile, bExt, bNoExt);
    utl::GetFileExtension(m_OrientationFile, orientationExt, orientationNoExt);
    
    for ( int i = 0; i < numberOfShells; i += 1 ) 
      {
      bVector = m_Sampling->GetBVectorInShell(i); 
      if (bVector->size()>0)
        {
        // m_Sampling has m_BVector
        double sum=0;
        for ( unsigned int j = 0; j < bVector->size(); j += 1 ) 
          sum += bVector->operator[](j);
        sum /= (double)bVector->size();
        int meanInt = (int)(sum);
        indexStr = utl::ConvertNumberToString(meanInt);
        }
      else
        {
        // m_Sampling has no m_BVector
        indexStr = "shell" + utl::ConvertNumberToString(i+1);
        }
      if (m_BFile!="")
        {
        realBFile = bNoExt + "_" + indexStr + "." + bExt;
        utl::SaveVector(*bVector, realBFile);
        std::cout << "save b values in shell " << i+1  << " to " << realOrientationFile << std::endl << std::flush;
        }
      orientationMatrix = m_Sampling->GetOrientationsCartesianInShell(i);
      realOrientationFile = orientationNoExt + "_" + indexStr + "." + orientationExt;
      std::cout << "save orientations in shell " << i+1  << " to " << realOrientationFile << std::endl << std::flush;
      orientationMatrix->Save(realOrientationFile);
      }
    }
}

}


#endif 

