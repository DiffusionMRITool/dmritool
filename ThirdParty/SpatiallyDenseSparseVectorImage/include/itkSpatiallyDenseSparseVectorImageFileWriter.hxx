/*=========================================================================
 
 Program:   Spatially Dense Sparse Vector Image File Writer
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImageFileWriter_hxx
#define __itkSpatiallyDenseSparseVectorImageFileWriter_hxx

#include <fstream>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSpatiallyDenseSparseVectorImageFileWriter.h"
#include "itksys/SystemTools.hxx"

namespace itk
{

//---------------------------------------------------------
template <class TInputImage>
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
::SpatiallyDenseSparseVectorImageFileWriter()
{
  m_FileName = "";
  m_UseCompression = true;
//  m_UseInputMetaDataDictionary = true;
}

//---------------------------------------------------------
template <class TInputImage>
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
::~SpatiallyDenseSparseVectorImageFileWriter()
{
}

//---------------------------------------------------------
template <class TInputImage>
void 
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
::SetInput(const InputImageType *input)
{
  // ProcessObject is not const_correct so this cast is required here.
  this->ProcessObject::SetNthInput(0, 
                                   const_cast<TInputImage *>(input ) );
}


//---------------------------------------------------------
template <class TInputImage>
const typename SpatiallyDenseSparseVectorImageFileWriter<TInputImage>::InputImageType *
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
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
const typename SpatiallyDenseSparseVectorImageFileWriter<TInputImage>::InputImageType *
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
::GetInput(unsigned int idx)
{
  return static_cast<TInputImage*> (this->ProcessObject::GetInput(idx));
}

//---------------------------------------------------------
template <class TInputImage>
void 
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
::Write()
{
  const InputImageType * input = this->GetInput();

  itkDebugMacro( <<"Writing an image file" );

  // Make sure input is available
  if ( input == 0 )
    {
    itkExceptionMacro(<< "No input to writer!");
    }

  if ( m_FileName == "" )
    {
    itkExceptionMacro(<<"No filename specified");
    }

  // write the data
  this->GenerateData();
}


//---------------------------------------------------------
template <class TInputImage>
void 
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
::GenerateData(void)
{
  itkDebugMacro ( << "SpatiallyDenseSparseVectorImageFileWriter::GenerateData() \n" );

  // Setup - Input Image
  InputImageType * input = const_cast<InputImageType*>(this->GetInput());
  
  // Total number of elements
  ImageRegionIteratorWithIndex<InputImageType>
    inputImageIterator( input, input->GetRequestedRegion() );

  inputImageIterator.GoToBegin();

  SizeValueType totalElements = 0;
  while ( !inputImageIterator.IsAtEnd() )
    {
    totalElements += (input->GetInternalPixel(inputImageIterator.GetIndex())).GetSize();
    ++inputImageIterator;
    }

  // Vector length
  unsigned int numberOfComponentsPerPixel = input->GetNumberOfComponentsPerPixel();

  // Setup - Output Images
  m_KeyImage = KeyImageType::New();
  m_ValueImage = ValueImageType::New();

  typename KeyImageType::IndexType startIndex;
  typename KeyImageType::RegionType region;
  typename KeyImageType::SizeType size;
  typename KeyImageType::SpacingType spacing;
  
  startIndex.Fill(0);
  
  if ( totalElements > 0 )
    {
    size[0] = totalElements;
    }
  else
    {
    size[0] = 1; // Output at least one voxel.
    }
    
  spacing.Fill(1);
  
  region.SetIndex(startIndex);
  region.SetSize(size);

  m_KeyImage->SetSpacing(spacing);
  m_KeyImage->SetRegions(region);
  m_ValueImage->SetSpacing(spacing);
  m_ValueImage->SetRegions(region);

  m_KeyImage->Allocate();
  m_ValueImage->Allocate();
  m_KeyImage->FillBuffer(static_cast<InputImageKeyType>(0));
  m_ValueImage->FillBuffer(static_cast<InputImageValueType>(0));

  // Iterators
  ImageRegionIteratorWithIndex<KeyImageType>
    keyImageIterator( m_KeyImage, m_KeyImage->GetRequestedRegion() );
  ImageRegionIteratorWithIndex<ValueImageType>
    valueImageIterator( m_ValueImage, m_ValueImage->GetRequestedRegion() );
  
  keyImageIterator.GoToBegin();
  valueImageIterator.GoToBegin();
  inputImageIterator.GoToBegin();

  InputImageIndexType index;
  SizeValueType offset;

  if ( totalElements > 0 )
    {
    while ( !inputImageIterator.IsAtEnd() )
      {
      index = inputImageIterator.GetIndex();

      const InputImagePixelPixelMapType *internalData =
          (input->GetInternalPixel(index)).GetDataPointer();
      InputImagePixelMapConstIteratorType iterator = internalData->begin();

      offset = input->ComputeOffset(index) * numberOfComponentsPerPixel;

      while ( iterator != internalData->end() )
        {
        keyImageIterator.Set(static_cast<InputImageKeyType>(iterator->first + offset));
        valueImageIterator.Set(iterator->second);
        ++iterator;

        ++keyImageIterator;
        ++valueImageIterator;
        }
      ++inputImageIterator;
      }
    }
  else
    {
    keyImageIterator.Set(static_cast<InputImageKeyType>(0));
    valueImageIterator.Set(0);
    }

  // Process File Names
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
    if (fileNameExtension != "spr")
      {
      std::cout << "Renaming extension to .spr" << std::endl;
      fileNameExtension = "spr";
      }
    baseFileName = fileName.substr(0, idx);
    }
  std::string keyFileName = baseFileName + "_key." + dataFileNameExtension;
  std::string valueFileName = baseFileName + "_value." + dataFileNameExtension;
  std::string headerFileName = baseFileName + "." + fileNameExtension;
  
  // Write Files
  m_KeyImageFileWriter = KeyImageFileWriterType::New();
  m_ValueImageFileWriter = ValueImageFileWriterType::New();

  std::string keyPathName = keyFileName;
  std::string valuePathName = valueFileName;
  std::string headerPathName = headerFileName;
  
  if ( pathName != "" )
    {
    keyPathName = pathName + "/" + keyFileName;
    valuePathName = pathName + "/" + valueFileName;
    headerPathName = pathName + "/" + headerFileName;
    }

  m_KeyImageFileWriter->SetFileName(keyPathName);
  m_KeyImageFileWriter->SetInput(m_KeyImage);
  m_ValueImageFileWriter->SetFileName(valuePathName);
  m_ValueImageFileWriter->SetInput(m_ValueImage);

  m_KeyImageFileWriter->SetUseCompression(this->m_UseCompression);
  m_ValueImageFileWriter->SetUseCompression(this->m_UseCompression);
//  m_KeyImageFileWriter->SetUseInputMetaDataDictionary(this->m_UseInputMetaDataDictionary);
//  m_ValueImageFileWriter->SetUseInputMetaDataDictionary(this->m_UseInputMetaDataDictionary);  
  
  m_KeyImageFileWriter->Update();
  m_ValueImageFileWriter->Update();

  // Write Header
  std::ofstream outfile;
//  std::string HeaderFileName = GetHeaderFileName();
//  std::cout << HeaderFileName << std::endl;
  
  outfile.open(headerPathName.c_str(), std::fstream::out);

  InputImageRegionType outputRegion = input->GetLargestPossibleRegion();
  InputImageSizeType outputSize = outputRegion.GetSize();
  InputImageSpacingType outputSpacing = input->GetSpacing();
  InputImagePointType outputOrigin = input->GetOrigin();
  InputImageDirectionType outputDirection = input->GetDirection();

  outfile << "NDims = " << input->GetImageDimension() + 1 << std::endl;
  outfile << "DimSize = ";
  outfile << input->GetNumberOfComponentsPerPixel() << " ";
  for (unsigned int d=0; d<input->GetImageDimension(); d++)
    {
    outfile << outputSize[d] << " ";
    }
  outfile << std::endl;
  
  outfile << "ElementSpacing = ";
  outfile << "1" << " ";
  for (unsigned int d=0; d<input->GetImageDimension(); d++)
    {
    outfile << outputSpacing[d] << " ";
    }
  outfile << std::endl;
  
  outfile << "Offset = ";
  outfile << "0" << " ";
  for (unsigned int d=0; d<input->GetImageDimension(); d++)
    {
    outfile << outputOrigin[d] << " ";
    }
  outfile << std::endl;
  
  outfile << "TransformMatrix = ";
  for (unsigned int d1=0; d1 < input->GetImageDimension() + 1; d1++)
    {
    for (unsigned int d2=0; d2 < input->GetImageDimension() + 1; d2++)
      {
      if ((d1 == 0) && (d2 == 0))
        {
        outfile << "1" << " ";
        }
      else if ((d1 == 0) && (d2 != 0))
        {
        outfile << "0" << " ";
        }
      else if ((d2 == 0) && (d1 != 0))
        {
        outfile << "0" << " ";
        }
      else
        {
        outfile << outputDirection(d1-1,d2-1) << " ";
        }
      }
    }
  outfile << std::endl;

  outfile << "KeyElementDataFile = " << keyFileName << std::endl;
  outfile << "ValueElementDataFile = " << valueFileName << std::endl;

  outfile.close();
}


//---------------------------------------------------------
template <class TInputImage>
void 
SpatiallyDenseSparseVectorImageFileWriter<TInputImage>
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
