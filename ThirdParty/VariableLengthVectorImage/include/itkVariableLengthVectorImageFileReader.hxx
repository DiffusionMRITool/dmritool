/*=========================================================================
 
 Program:   Variable Length Vector Image File Reader
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef __itkVariableLengthVectorImageFileReader_hxx
#define __itkVariableLengthVectorImageFileReader_hxx

#include "itkVariableLengthVectorImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itksys/SystemTools.hxx"

namespace itk
{

template <class TOutputImage>
VariableLengthVectorImageFileReader<TOutputImage>
::VariableLengthVectorImageFileReader()
{
  m_FileName = "";
  m_UseStreaming = true;
  m_ImageIO = 0;
}

template <class TOutputImage>
VariableLengthVectorImageFileReader<TOutputImage>
::~VariableLengthVectorImageFileReader()
{
}

template <class TOutputImage>
void VariableLengthVectorImageFileReader<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  if (m_ImageIO)
    {
    os << indent << "ImageIO: \n";
    m_ImageIO->Print(os, indent.GetNextIndent());
    }
  else
    {
    os << indent << "ImageIO: (null)" << "\n";
    }
  
  os << indent << "m_FileName: " << m_FileName << "\n";
  os << indent << "m_UseStreaming: " << m_UseStreaming << "\n";
}

template <class TOutputImage>
void VariableLengthVectorImageFileReader<TOutputImage>
::GenerateData()
{
  itkDebugMacro ( << "VariableLengthVectorImageFileReader::GenerateData() \n" );

  if (m_FileName.length()<4 || m_FileName.compare(m_FileName.length() - 4, 4, ".vlv")!=0)
    {
    itkExceptionMacro( << "the file " << m_FileName  << " should be a text file ending with '.vlv'.");
    }

  std::string lengthFileName;
  std::string valueFileName;
  char tempLine[256];
  
  // Read Header
  std::ifstream infile;
  infile.open(m_FileName.c_str(), std::fstream::in);

  if ( !infile.is_open() )
    {
    itkExceptionMacro( << "Cannot open file: " << m_FileName );
    }
      
  std::string fileName = itksys::SystemTools::GetFilenameName( m_FileName );
  std::string pathName = itksys::SystemTools::GetFilenamePath( m_FileName );

  std::string line, extractedLine;

  OutputImageRegionType outputRegion;
  OutputImageSizeType outputSize;
  OutputImageSpacingType outputSpacing;
  OutputImageIndexType outputStartIndex;
  OutputImagePointType outputOrigin;
  OutputImageDirectionType outputDirection;

  unsigned int ndims = 0;
//  unsigned int numberOfComponentsPerPixel = 0;
  
  char * pch;
  unsigned int idx;
  
  while (!infile.eof())
    {
    if (infile.getline(tempLine, 256))
      {
      line = std::string(tempLine);

      if (line.find("NDims") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);
        pch = strtok(const_cast<char*>(extractedLine.c_str())," ");
        if (pch != NULL)
          {
          ndims = atoi(pch);
          }
        }

      if (line.find("DimSize") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);
        idx = 0;
        pch = strtok (const_cast<char*>(extractedLine.c_str())," ");
        while (pch != NULL)
          {
          outputSize[idx] = atoi(pch);
          idx++;
          pch = strtok (NULL, " ");
          }
        }

      if (line.find("ElementSpacing") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);
        idx = 0;
        pch = strtok (const_cast<char*>(extractedLine.c_str())," ");
        while (pch != NULL)
          {
          outputSpacing[idx] = atof(pch);
          idx++;
          pch = strtok (NULL, " ");
          }
        }

      if (line.find("Offset") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);
        idx = 0;
        pch = strtok (const_cast<char*>(extractedLine.c_str())," ");
        while (pch != NULL)
          {
          outputOrigin[idx] = atof(pch);
          idx++;
          pch = strtok (NULL, " ");
          }
        }

      if (line.find("TransformMatrix") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);
        idx = 0;
        pch = strtok (const_cast<char*>(extractedLine.c_str())," ");
        outputDirection.SetIdentity();
        while (pch != NULL)
          {
          outputDirection( vcl_floor(idx / ImageDimension),
                           idx % ImageDimension ) = atof(pch);
          idx++;
          pch = strtok (NULL, " ");
          }
        }

      if (line.find("LengthElementDataFile") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);

        const size_t beginStr = extractedLine.find_first_not_of(" \t");
        if (beginStr == std::string::npos)
          {
          // no content
          std::cerr << "No length file specified!";
          }
        const size_t endStr = extractedLine.find_last_not_of(" \t");
        const size_t range = endStr - beginStr + 1;
        lengthFileName = extractedLine.substr(beginStr, range);
        if ( pathName != "" )
          {
          lengthFileName = pathName + "/" + lengthFileName;
          }
        }

      if (line.find("ValueElementDataFile") != std::string::npos)
         {
         extractedLine = line.substr(line.find("=") + 1);

         const size_t beginStr = extractedLine.find_first_not_of(" \t");
         if (beginStr == std::string::npos)
           {
           // no content
           std::cerr << "No value file specified!";
           }
         const size_t endStr = extractedLine.find_last_not_of(" \t");
         const size_t range = endStr - beginStr + 1;
         valueFileName = extractedLine.substr(beginStr, range);
         if ( pathName != "" )
           {
           valueFileName = pathName + "/" + valueFileName;
           }
         }
      }
    else
      {
      break;
      }

    }
  infile.close();
   
  // Setup output Image
  OutputImageType * output = this->GetOutput();
  outputStartIndex.Fill(0);
  outputRegion.SetIndex(outputStartIndex);
  outputRegion.SetSize(outputSize);
  output->SetRegions(outputRegion);
  output->SetSpacing(outputSpacing);
//  outputOrigin.Fill(0);
  output->SetOrigin(outputOrigin);
//  outputDirection.SetIdentity();
  output->SetDirection(outputDirection);
  output->Allocate();

  // OutputImagePixelType outputPixel;
  // outputPixel.Clear();
  
  // output->FillBuffer(outputPixel);

  m_LengthImageFileReader = LengthImageFileReaderType::New();
  m_ValueImageFileReader = ValueImageFileReaderType::New();
  
  m_LengthImageFileReader->SetFileName(lengthFileName);
  m_ValueImageFileReader->SetFileName(valueFileName);
  
  m_LengthImageFileReader->Update();
  m_ValueImageFileReader->Update();
  
  m_LengthImage = m_LengthImageFileReader->GetOutput();
  m_ValueImage = m_ValueImageFileReader->GetOutput();
  
  m_LengthImage->DisconnectPipeline();
  m_ValueImage->DisconnectPipeline();

  m_ImageIO = m_ValueImageFileReader->GetImageIO();

  ImageRegionIterator<OutputImageType>
    outputIt( output, output->GetLargestPossibleRegion() );
  ImageRegionIterator<LengthImageType>
    lengthIt( m_LengthImage, m_LengthImage->GetRequestedRegion() );
  ImageRegionIterator<ValueImageType>
    valueIt( m_ValueImage, m_ValueImage->GetRequestedRegion() );

  // Populate data
  outputIt.GoToBegin();
  lengthIt.GoToBegin();
  valueIt.GoToBegin();
  while ( !outputIt.IsAtEnd() )
    {
    LengthType length = static_cast<LengthType>(lengthIt.Get());
    OutputImagePixelType outputPixel;
    outputPixel.SetSize(length);
    for (LengthType k=0; k<length; k++)
      {
      outputPixel[k] = valueIt.Get();
      ++valueIt;
      }
    outputIt.Set(outputPixel);

    ++outputIt;
    ++lengthIt;
    }
}


} //namespace ITK

#endif
