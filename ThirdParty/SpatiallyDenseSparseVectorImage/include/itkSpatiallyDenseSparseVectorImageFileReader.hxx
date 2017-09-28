/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image File Reader

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImageFileReader_hxx
#define __itkSpatiallyDenseSparseVectorImageFileReader_hxx

#include "itkSpatiallyDenseSparseVectorImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itksys/SystemTools.hxx"


namespace itk
{

template <class TOutputImage>
SpatiallyDenseSparseVectorImageFileReader<TOutputImage>
::SpatiallyDenseSparseVectorImageFileReader()
{
  m_FileName = "";
  m_UseStreaming = true;
  m_ImageIO = 0;
}

template <class TOutputImage>
SpatiallyDenseSparseVectorImageFileReader<TOutputImage>
::~SpatiallyDenseSparseVectorImageFileReader()
{
}

template <class TOutputImage>
void SpatiallyDenseSparseVectorImageFileReader<TOutputImage>
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
void SpatiallyDenseSparseVectorImageFileReader<TOutputImage>
::GenerateData()
{
  itkDebugMacro ( << "SpatiallyDenseSparseVectorImageFileReader::GenerateData() \n" );

  std::string keyFileName;
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
  unsigned int numberOfComponentsPerPixel = 0;

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
          ndims--;
          }
        }

      if (line.find("DimSize") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);
        idx = 0;
        pch = strtok (const_cast<char*>(extractedLine.c_str())," ");
        while (pch != NULL)
          {
          if ( idx == 0 )
            {
            numberOfComponentsPerPixel = atoi(pch);
            }
          else
            {
            outputSize[idx-1] = atoi(pch);
            }
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
          if ( idx > 0 )
            {
            outputSpacing[idx-1] = atof(pch);
            }
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
          if ( idx > 0 )
            {
            outputOrigin[idx-1] = atof(pch);
            }
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
          unsigned int row = vcl_floor(static_cast<double>(idx)/static_cast<double>(ImageDimension + 1));
          unsigned int col = idx % (ImageDimension + 1);
          if ( (row != 0) && (col != 0) )
            {
            outputDirection(row - 1, col - 1) = atof(pch);
            }
          idx++;
          pch = strtok (NULL, " ");
          }
        }

      if (line.find("KeyElementDataFile") != std::string::npos)
        {
        extractedLine = line.substr(line.find("=") + 1);

        const size_t beginStr = extractedLine.find_first_not_of(" \t");
        if (beginStr == std::string::npos)
        {
          // no content
          std::cerr << "No key file specified!";
        }
        const size_t endStr = extractedLine.find_last_not_of(" \t");
        const size_t range = endStr - beginStr + 1;
        keyFileName = extractedLine.substr(beginStr, range);
        if ( pathName != "" )
          {
          keyFileName = pathName + "/" + keyFileName;
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

  // Setup - Output Image
  OutputImageType * output = this->GetOutput();
  outputStartIndex.Fill(0);
  outputRegion.SetIndex(outputStartIndex);
  outputRegion.SetSize(outputSize);
  output->SetRegions(outputRegion);
  output->SetSpacing(outputSpacing);
  // outputOrigin.Fill(0);
  output->SetOrigin(outputOrigin);
  // outputDirection.SetIdentity();
  output->SetDirection(outputDirection);
  output->SetNumberOfComponentsPerPixel(numberOfComponentsPerPixel);
  output->Allocate();

  OutputImagePixelType outputPixel;
  outputPixel.SetSize(numberOfComponentsPerPixel);
  outputPixel.Fill(0);

  output->FillBuffer(outputPixel);

  OutputImagePixelContainerType * container = output->GetPixelContainer();

  m_KeyImageFileReader = KeyImageFileReaderType::New();
  m_ValueImageFileReader = ValueImageFileReaderType::New();

  m_KeyImageFileReader->SetFileName(keyFileName);
  m_ValueImageFileReader->SetFileName(valueFileName);

  m_KeyImageFileReader->Update();
  m_ValueImageFileReader->Update();

  m_KeyImage = m_KeyImageFileReader->GetOutput();
  m_ValueImage = m_ValueImageFileReader->GetOutput();

  m_ImageIO = m_ValueImageFileReader->GetImageIO();

  // Iterators
  ImageRegionIterator<KeyImageType>
    keyImageIterator( m_KeyImage, m_KeyImage->GetRequestedRegion() );
  ImageRegionIterator<ValueImageType>
    valueImageIterator( m_ValueImage, m_ValueImage->GetRequestedRegion() );

  keyImageIterator.GoToBegin();
  valueImageIterator.GoToBegin();

  // Populate Data
  SizeValueType offset;
  SizeValueType elementIndex;
  typename KeyImageType::PixelType key;
  while ( !keyImageIterator.IsAtEnd() )
    {
    key = keyImageIterator.Get();
    offset = (int)key / numberOfComponentsPerPixel;
    elementIndex = key % numberOfComponentsPerPixel;
    (*container)[offset][elementIndex] = valueImageIterator.Get();

    ++keyImageIterator;
    ++valueImageIterator;
    }

}


} //namespace ITK

#endif
