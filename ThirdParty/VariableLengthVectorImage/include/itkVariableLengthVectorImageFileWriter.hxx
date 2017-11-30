/*=========================================================================
 
 Program:   Variable Length Vector Image File Writer
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef __itkVariableLengthVectorImageFileWriter_hxx
#define __itkVariableLengthVectorImageFileWriter_hxx

#include <fstream>
#include "itkVariableLengthVectorImageFileWriter.h"
#include "itksys/SystemTools.hxx"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


namespace itk
{

//---------------------------------------------------------
template <class TInputImage>
VariableLengthVectorImageFileWriter<TInputImage>
::VariableLengthVectorImageFileWriter()
{
  m_FileName = "";
  m_UseCompression = true;
//  m_UseInputMetaDataDictionary = true;
}

//---------------------------------------------------------
template <class TInputImage>
VariableLengthVectorImageFileWriter<TInputImage>
::~VariableLengthVectorImageFileWriter()
{
}

//---------------------------------------------------------
template <class TInputImage>
void 
VariableLengthVectorImageFileWriter<TInputImage>
::SetInput(const InputImageType *input)
{
  // ProcessObject is not const_correct so this cast is required here.
  this->ProcessObject::SetNthInput(0, 
                                   const_cast<TInputImage *>(input ) );
}


//---------------------------------------------------------
template <class TInputImage>
const typename VariableLengthVectorImageFileWriter<TInputImage>::InputImageType *
VariableLengthVectorImageFileWriter<TInputImage>
::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
    {
    return 0;
    }
  
  return static_cast<TInputImage*>
    (this->ProcessObject::GetInput(0));
}
  
//---------------------------------------------------------
template <class TInputImage>
const typename VariableLengthVectorImageFileWriter<TInputImage>::InputImageType *
VariableLengthVectorImageFileWriter<TInputImage>
::GetInput(unsigned int idx)
{
  return static_cast<TInputImage*> (this->ProcessObject::GetInput(idx));
}

//---------------------------------------------------------
template <class TInputImage>
void 
VariableLengthVectorImageFileWriter<TInputImage>
::Write()
{
  const InputImageType * input = this->GetInput();

  itkDebugMacro( <<"Writing an image file" );

  // Make sure input is available
  if ( input == 0 )
    {
    itkExceptionMacro( << "No input to writer!" );
    }

  if ( m_FileName == "" )
    {
    itkExceptionMacro( <<"No filename specified" );
    }

  // write the data
  this->GenerateData();
}


//---------------------------------------------------------
template <class TInputImage>
void 
VariableLengthVectorImageFileWriter<TInputImage>
::GenerateData(void)
{
  itkDebugMacro ( << "VariableLengthVectorImageFileWriter::GenerateData() \n" );

  if (m_FileName.length()<4 || m_FileName.compare(m_FileName.length() - 4, 4, ".vlv")!=0)
    {
    itkExceptionMacro( << "the file " << m_FileName  << " should be a text file ending with '.vlv'.");
    }

  // Setup - Input Image
  InputImageType * input = const_cast<InputImageType*>(this->GetInput());
  
  // Iterator
  ImageRegionConstIterator<InputImageType>
    inputIt( input, input->GetLargestPossibleRegion() );

  // Compute some image statistics
  inputIt.GoToBegin();
  unsigned long totalLength = 0;

  while( !inputIt.IsAtEnd() )
    {
    totalLength += inputIt.Get().Size();
    ++inputIt;
    }

  // Setup output Images
  m_LengthImage = LengthImageType::New();
  m_LengthImage->CopyInformation(input);
  m_LengthImage->SetRegions(input->GetLargestPossibleRegion());
  m_LengthImage->Allocate();

  typename ValueImageType::IndexType startIndex;
  typename ValueImageType::RegionType region;
  typename ValueImageType::SizeType size;
  typename ValueImageType::SpacingType spacing;
  
  startIndex.Fill(0);
  spacing.Fill(1);
  size[0] = totalLength;
  region.SetIndex(startIndex);
  region.SetSize(size);

  m_ValueImage = ValueImageType::New();
  m_ValueImage->SetSpacing(spacing);
  m_ValueImage->SetRegions(region);
  m_ValueImage->Allocate();

  // More iterators
  ImageRegionIterator<LengthImageType>
    lengthIt( m_LengthImage, m_LengthImage->GetRequestedRegion() );
  ImageRegionIterator<ValueImageType>
    valueIt( m_ValueImage, m_ValueImage->GetRequestedRegion() );

  lengthIt.GoToBegin();
  valueIt.GoToBegin();

  // Populate data
  inputIt.GoToBegin();
  while ( !inputIt.IsAtEnd() )
    {
    InputImagePixelType inputPixel = inputIt.Get();
    LengthType length = static_cast<LengthType>(inputPixel.Size());
    lengthIt.Set(length);
    for (LengthType k=0; k<length; k++)
      {
      valueIt.Set(inputPixel[k]);
      ++valueIt;
      }
    ++inputIt;
    ++lengthIt;
    }
  
  // Process file names
  std::string baseFileName = "";
  std::string fileNameExtension;
  std::string dataFileNameExtension = "nrrd";

  std::string fileName = itksys::SystemTools::GetFilenameName( m_FileName );
  std::string pathName = itksys::SystemTools::GetFilenamePath( m_FileName );

  std::string::size_type idx;
  idx = fileName.find_last_of('.');

  if (idx != std::string::npos)
    {
    fileNameExtension = fileName.substr(idx + 1);
    if (fileNameExtension != "vlv")
      {
      std::cout << "Renaming extension to .vlv" << std::endl;
      fileNameExtension = "vlv";
      }
    baseFileName = fileName.substr(0, idx);
    }
  std::string lengthFileName = baseFileName + "_length." + dataFileNameExtension;
  std::string valueFileName = baseFileName + "_value." + dataFileNameExtension;
  std::string headerFileName = baseFileName + "." + fileNameExtension;
  
  // Write files
  m_LengthImageFileWriter = LengthImageFileWriterType::New();
  m_ValueImageFileWriter = ValueImageFileWriterType::New();

  std::string lengthPathName = lengthFileName;
  std::string valuePathName = valueFileName;
  std::string headerPathName = headerFileName;
  
  if ( pathName != "" )
    {
    lengthPathName = pathName + "/" + lengthFileName;
    valuePathName = pathName + "/" + valueFileName;
    headerPathName = pathName + "/" + headerFileName;
    }

  m_LengthImageFileWriter->SetFileName(lengthPathName);
  m_LengthImageFileWriter->SetInput(m_LengthImage);
  m_ValueImageFileWriter->SetFileName(valuePathName);
  m_ValueImageFileWriter->SetInput(m_ValueImage);

  m_LengthImageFileWriter->SetUseCompression(this->m_UseCompression);
  m_ValueImageFileWriter->SetUseCompression(this->m_UseCompression);
//  m_LengthImageFileWriter->SetUseInputMetaDataDictionary(this->m_UseInputMetaDataDictionary);
//  m_ValueImageFileWriter->SetUseInputMetaDataDictionary(this->m_UseInputMetaDataDictionary);  
  
  m_LengthImageFileWriter->Update();
  m_ValueImageFileWriter->Update();

  // Write header
  std::ofstream outfile;
//  std::string HeaderFileName = GetHeaderFileName();
//  std::cout << HeaderFileName << std::endl;
  
  outfile.open(headerPathName.c_str(), std::fstream::out);

  InputImageRegionType outputRegion = input->GetLargestPossibleRegion();
  InputImageSizeType outputSize = outputRegion.GetSize();
  InputImageSpacingType outputSpacing = input->GetSpacing();
  InputImagePointType outputOrigin = input->GetOrigin();
  InputImageDirectionType outputDirection = input->GetDirection();

  outfile << "NDims = " << input->GetImageDimension() << std::endl;
  outfile << "DimSize = ";
//  outfile << input->GetNumberOfComponentsPerPixel() << " ";
  for (unsigned int d=0; d<input->GetImageDimension(); d++)
    {
    outfile << outputSize[d] << " ";
    }
  outfile << std::endl;
  
  outfile << "ElementSpacing = ";
//  outfile << "1" << " ";
  for (unsigned int d=0; d<input->GetImageDimension(); d++)
    {
    outfile << outputSpacing[d] << " ";
    }
  outfile << std::endl;

  outfile << "Offset = ";
//  outfile << "0" << " ";
  for (unsigned int d=0; d<input->GetImageDimension(); d++)
    {
    outfile << outputOrigin[d] << " ";
    }
  outfile << std::endl;

  outfile << "TransformMatrix = ";
  for (unsigned int d1=0; d1<input->GetImageDimension(); d1++)
    {
    for (unsigned int d2=0; d2<input->GetImageDimension(); d2++)
      outfile << outputDirection(d1,d2) << " ";
    }
  outfile << std::endl;

  outfile << "LengthElementDataFile = " << lengthFileName << std::endl;
  outfile << "ValueElementDataFile = " << valueFileName << std::endl;
  
  outfile.close();
}


//---------------------------------------------------------
template <class TInputImage>
void 
VariableLengthVectorImageFileWriter<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (m_FileName.data() ? m_FileName.data() : "(none)") << std::endl;

  if (m_UseCompression)
    {
    os << indent << "Compression: On\n";
    }
  else
    {
    os << indent << "Compression: Off\n";
    }

//  if (m_UseInputMetaDataDictionary)
//    {
//    os << indent << "UseInputMetaDataDictionary: On\n";
//    }
//  else
//    {
//    os << indent << "UseInputMetaDataDictionary: Off\n";
//    }

}

} // end namespace itk

#endif
