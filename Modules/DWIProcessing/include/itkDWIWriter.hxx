/**
 *       @file  itkDWIWriter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-25-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkDWIWriter_hxx
#define __itkDWIWriter_hxx

#include "utl.h"


namespace itk
{

template <class TPixelType, unsigned int VImageDimension>
DWIWriter<TPixelType, VImageDimension>
::DWIWriter() 
{
  m_UseRelativePath = true;
  m_OutputEachShell = false;

  m_MaskImage = B0ImageType::New();
  m_B0Image = MaskImageType::New();

  m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIWriter<TPixelType, VImageDimension>
::SetInput(const DWIImageType *input)
{
  // ProcessObject is not const_correct so this cast is required here.
  this->ProcessObject::SetNthInput( 0, const_cast< DWIImageType * >( input ) );
}

template <class TPixelType, unsigned int VImageDimension>
const typename DWIWriter<TPixelType, VImageDimension>::DWIImageType*
DWIWriter<TPixelType, VImageDimension>
::GetInput()
{
  if (this->GetNumberOfInputs() < 1)
    {
    return 0;
    }
  
  return static_cast<DWIImageType*>
    (this->ProcessObject::GetInput(0));
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIWriter<TPixelType, VImageDimension>
::WriteToConfigurationFile(const std::string file)
{
  std::ofstream  out;                                // create ofstream object
  
  out.open ( file.c_str() );           // open ofstream
  utlException (!out, "\nERROR : failed to open output file " << file );
  
  typename DWIImageType::Pointer dwi = const_cast<DWIImageType*>(this->GetInput());
  std::string folderPath = "", fileNoPath;
  if (!m_UseRelativePath)
    {
    utl::GetPath(file, folderPath, fileNoPath);
    }

  MatrixPointer orientationsCartesian = m_SamplingSchemeQSpace->GetOrientationsCartesian();
  if (m_OutputEachShell)
    {
    std::string BFile_ext, BFile_noext;
    utl::GetFileExtension(m_BFile, BFile_ext, BFile_noext);
    std::string DWIFile_ext, DWIFile_noext, DWIFile_final;
    utl::GetFileExtension(m_DWIFile, DWIFile_ext, DWIFile_noext);
    std::string OrientationFile_ext, OrientationFile_noext, OrientationFile_final;
    utl::GetFileExtension(m_OrientationFile, OrientationFile_ext, OrientationFile_noext);

    std::vector<STDVectorType> bVectors = m_SamplingSchemeQSpace->GroupBValues(); 
    typename SamplingSchemeQSpaceType::Index2DVectorPointer bIndices = m_SamplingSchemeQSpace->GetIndicesInShells();

    for ( int j = 0; j < bVectors.size(); j += 1 ) 
      {
      double bMean=0;
      for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
        bMean += bVectors[j][k];
      bMean /= bVectors[j].size();

      char outTemp[1024];
      sprintf(outTemp, "_b%d.", utl::RoundNumber(bMean));
      DWIFile_final = folderPath+ DWIFile_noext + outTemp + DWIFile_ext;
      OrientationFile_final = folderPath + OrientationFile_noext + outTemp + OrientationFile_ext;
      out << bMean << " " << OrientationFile_final << " " << DWIFile_final << std::endl;
      std::cout << bMean << " " << OrientationFile_final << " " << DWIFile_final << ", " << bVectors[j].size() << " dwis" << std::endl;

      MatrixType gradTemp(bVectors[j].size(),3);
      for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
        {
        gradTemp(k,0) = (*orientationsCartesian)((*bIndices)[j][k],0);
        gradTemp(k,1) = (*orientationsCartesian)((*bIndices)[j][k],1);
        gradTemp(k,2) = (*orientationsCartesian)((*bIndices)[j][k],2);
        }
      gradTemp.Save(OrientationFile_final);
      
      typename DWIImageType::Pointer dwiTempImage = DWIImageType::New();
      itk::CopyImageInformation<DWIImageType, DWIImageType>( dwi, dwiTempImage);
      dwiTempImage->SetNumberOfComponentsPerPixel(bVectors[j].size());
      dwiTempImage->Allocate();
      ImageRegionIteratorWithIndex<DWIImageType> dwiIt(dwi, dwi->GetLargestPossibleRegion());
      ImageRegionIteratorWithIndex<DWIImageType> dwiTempIt(dwiTempImage, dwiTempImage->GetLargestPossibleRegion());
      ImageRegionIteratorWithIndex<MaskImageType> maskIt;
      PixelType dwiPixel, dwiTempPixel, zeroPixel;
      dwiTempPixel.SetSize(bVectors[j].size());
      zeroPixel.SetSize(bVectors[j].size());
      zeroPixel.Fill(0.0);
      if (!IsImageEmpty(m_MaskImage))
        maskIt = ImageRegionIteratorWithIndex<MaskImageType>(this->m_MaskImage, m_MaskImage->GetLargestPossibleRegion());
      for (dwiTempIt.GoToBegin(), dwiIt.GoToBegin(), maskIt.GoToBegin(); 
        !dwiIt.IsAtEnd(); 
        ++dwiTempIt, ++dwiIt, ++maskIt) 
        {
        if (!IsImageEmpty(m_MaskImage) && maskIt.Get()<=1e-10)
          {
          dwiTempIt.Set(zeroPixel);
          continue;
          }
        dwiPixel = dwiIt.Get();
        for ( int k = 0; k < bVectors[j].size(); k += 1 ) 
          {
          dwiTempPixel[k] = dwiPixel[ (*bIndices)[j][k] ];
          }
        dwiTempIt.Set(dwiTempPixel);
        }
      itk::SaveImage(dwiTempImage, DWIFile_final);
      }

    }
  else
    {
    STDVectorPointer bVector = m_SamplingSchemeQSpace->GetBVector();
    std::cout << "Write b vector to " << folderPath + m_BFile << std::endl << std::flush;
    utl::SaveVector(*bVector, folderPath + m_BFile);
    std::cout << "Write orientation file to " << folderPath + m_OrientationFile << std::endl << std::flush;
    orientationsCartesian->Save(folderPath + m_OrientationFile);
    itk::SaveImage(dwi, folderPath + m_DWIFile, "Write DWI image to");
    out << m_BFile << " " << m_OrientationFile << " " << m_DWIFile << std::endl;
    }
  
  out.close();
  return;
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIWriter<TPixelType, VImageDimension>
::GenerateData()
{
  if (m_ConfigurationFile!="")
    WriteToConfigurationFile(m_ConfigurationFile);
  else
    utlGlobalException(true, "unknown format! ");
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIWriter<TPixelType, VImageDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "m_SamplingSchemeQSpace = " << m_SamplingSchemeQSpace << std::endl << std::flush;
  PrintVar4(true, m_ConfigurationFile, m_DWIFile, m_BFile, m_OrientationFile, os<<indent);
  PrintVar2(true, m_UseRelativePath, m_OutputEachShell, os<<indent);
}

}


#endif 
