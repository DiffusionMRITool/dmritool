/*=========================================================================
 
 Program:   Cast Vector Image File Writer
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef __itkCastVectorImageFileWriter_hxx
#define __itkCastVectorImageFileWriter_hxx

#include "itkImageFileWriter.h"
#include "itkDataObject.h"
#include "itkObjectFactoryBase.h"
#include "itkImageIOFactory.h"
#include "itkCommand.h"
#include "vnl/vnl_vector.h"
#include "itkVectorImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"

namespace itk
{

//---------------------------------------------------------
template <class TInputImage>
CastVectorImageFileWriter<TInputImage>
::CastVectorImageFileWriter():
  m_PasteIORegion(TInputImage::ImageDimension)
{
  m_UseCompression = false;
  m_UseInputMetaDataDictionary = true;
  m_FactorySpecifiedImageIO = false;
  m_UserSpecifiedIORegion = false;
  m_UserSpecifiedImageIO = false;
  m_NumberOfStreamDivisions = 1;
  m_ComponentType = ImageIOBase::DOUBLE;
}

//---------------------------------------------------------
template <class TInputImage>
CastVectorImageFileWriter<TInputImage>
::~CastVectorImageFileWriter()
{
}

//---------------------------------------------------------
template <class TInputImage>
void 
CastVectorImageFileWriter<TInputImage>
::SetInput(const InputImageType *input)
{
  // ProcessObject is not const_correct so this cast is required here.
  this->ProcessObject::SetNthInput(0, 
                                   const_cast<TInputImage *>(input ) );
}


//---------------------------------------------------------
template <class TInputImage>
const typename CastVectorImageFileWriter<TInputImage>::InputImageType *
CastVectorImageFileWriter<TInputImage>
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
const typename CastVectorImageFileWriter<TInputImage>::InputImageType *
CastVectorImageFileWriter<TInputImage>
::GetInput(unsigned int idx)
{
  return static_cast<TInputImage*> (this->ProcessObject::GetInput(idx));
}

//---------------------------------------------------------
template <class TInputImage>
void 
CastVectorImageFileWriter<TInputImage>
::SetIORegion (const ImageIORegion& region) 
{
  itkDebugMacro("setting IORegion to " << region );
  if ( m_PasteIORegion != region)
    {
    m_PasteIORegion = region;
    this->Modified();
    m_UserSpecifiedIORegion = true;
    }
} 

//---------------------------------------------------------
template <class TInputImage>
void 
CastVectorImageFileWriter<TInputImage>
::Write()
{
  itkDebugMacro( <<"Writing an image file" );
  this->GenerateData();   
}


//---------------------------------------------------------
template <class TInputImage>
void 
CastVectorImageFileWriter<TInputImage>
::GenerateData(void)
{
  const InputImageType * inputPtr = this->GetInput();  

  itkDebugMacro(<<"Writing file: " << m_FileName);

//  if( strcmp( input->GetNameOfClass(), "VectorImage" ) == 0 ) 
  
  // Writer typedefs
  typedef VectorImage<unsigned char, TInputImage::ImageDimension>    UCharImageType;
  typedef VectorImage<char, TInputImage::ImageDimension>             CharImageType;
  typedef VectorImage<unsigned short, TInputImage::ImageDimension>   UShortImageType;
  typedef VectorImage<short, TInputImage::ImageDimension>            ShortImageType;
  typedef VectorImage<unsigned int, TInputImage::ImageDimension>     UIntImageType;
  typedef VectorImage<int, TInputImage::ImageDimension>              IntImageType;
  typedef VectorImage<unsigned long, TInputImage::ImageDimension>    ULongImageType;
  typedef VectorImage<long, TInputImage::ImageDimension>             LongImageType;
  typedef VectorImage<float, TInputImage::ImageDimension>            FloatImageType;
  typedef VectorImage<double, TInputImage::ImageDimension>           DoubleImageType;

  typedef ImageFileWriter<UCharImageType>    UCharWriterType;
  typedef ImageFileWriter<CharImageType>     CharWriterType;
  typedef ImageFileWriter<UShortImageType>   UShortWriterType;
  typedef ImageFileWriter<ShortImageType>    ShortWriterType;
  typedef ImageFileWriter<UIntImageType>     UIntWriterType;
  typedef ImageFileWriter<IntImageType>      IntWriterType;
  typedef ImageFileWriter<ULongImageType>    ULongWriterType;
  typedef ImageFileWriter<LongImageType>     LongWriterType;
  typedef ImageFileWriter<FloatImageType>    FloatWriterType;
  typedef ImageFileWriter<DoubleImageType>   DoubleWriterType;

  typedef CastImageFilter<TInputImage, UCharImageType>   UCharCastFilterType;
  typedef CastImageFilter<TInputImage, CharImageType>    CharCastFilterType;
  typedef CastImageFilter<TInputImage, UShortImageType>  UShortCastFilterType;
  typedef CastImageFilter<TInputImage, ShortImageType>   ShortCastFilterType;
  typedef CastImageFilter<TInputImage, UIntImageType>    UIntCastFilterType;
  typedef CastImageFilter<TInputImage, IntImageType>     IntCastFilterType;
  typedef CastImageFilter<TInputImage, ULongImageType>   ULongCastFilterType;
  typedef CastImageFilter<TInputImage, LongImageType>    LongCastFilterType;
  typedef CastImageFilter<TInputImage, FloatImageType>   FloatCastFilterType;
  typedef CastImageFilter<TInputImage, DoubleImageType>  DoubleCastFilterType;
   
  typename UCharWriterType::Pointer ucharWriter; 
  typename CharWriterType::Pointer charWriter;
  typename UShortWriterType::Pointer ushortWriter;
  typename ShortWriterType::Pointer shortWriter;
  typename UIntWriterType::Pointer uintWriter;
  typename IntWriterType::Pointer intWriter;
  typename ULongWriterType::Pointer ulongWriter;
  typename LongWriterType::Pointer longWriter;
  typename FloatWriterType::Pointer floatWriter;
  typename DoubleWriterType::Pointer doubleWriter;
  
  typename UCharCastFilterType::Pointer ucharCastFilter;
  typename CharCastFilterType::Pointer charCastFilter;
  typename UShortCastFilterType::Pointer ushortCastFilter;
  typename ShortCastFilterType::Pointer shortCastFilter;
  typename UIntCastFilterType::Pointer uintCastFilter;
  typename IntCastFilterType::Pointer intCastFilter;
  typename ULongCastFilterType::Pointer ulongCastFilter;
  typename LongCastFilterType::Pointer longCastFilter;
  typename FloatCastFilterType::Pointer floatCastFilter;
  typename DoubleCastFilterType::Pointer doubleCastFilter;
  
  switch ( m_ComponentType )
    {
    case ImageIOBase::UCHAR:
      ucharWriter = UCharWriterType::New();
      ucharCastFilter = UCharCastFilterType::New();
      ucharWriter->SetFileName( m_FileName );
      ucharCastFilter->SetInput( inputPtr );
      ucharWriter->SetInput( ucharCastFilter->GetOutput() );
      ucharWriter->SetImageIO( this->GetImageIO() );
      ucharWriter->SetIORegion( this->GetIORegion() );
      ucharWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      ucharWriter->SetUseCompression( this->GetUseCompression() );
      ucharWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      ucharWriter->Update();
      break;
    case ImageIOBase::CHAR:
      charWriter = CharWriterType::New();
      charCastFilter = CharCastFilterType::New();
      charWriter->SetFileName( m_FileName );
      charCastFilter->SetInput( inputPtr );
      charWriter->SetInput( charCastFilter->GetOutput() );
      charWriter->SetImageIO( this->GetImageIO() );
      charWriter->SetIORegion( this->GetIORegion() );
      charWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      charWriter->SetUseCompression( this->GetUseCompression() );
      charWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      charWriter->Update();
      break;
    case ImageIOBase::USHORT:
      ushortWriter = UShortWriterType::New();
      ushortCastFilter = UShortCastFilterType::New();
      ushortWriter->SetFileName( m_FileName );
      ushortCastFilter->SetInput( inputPtr );
      ushortWriter->SetInput( ushortCastFilter->GetOutput() );
      ushortWriter->SetImageIO( this->GetImageIO() );
      ushortWriter->SetIORegion( this->GetIORegion() );
      ushortWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      ushortWriter->SetUseCompression( this->GetUseCompression() );
      ushortWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      ushortWriter->Update();
      break;
    case ImageIOBase::SHORT:
      shortWriter = ShortWriterType::New();
      shortCastFilter = ShortCastFilterType::New();
      shortWriter->SetFileName( m_FileName );
      shortCastFilter->SetInput( inputPtr );
      shortWriter->SetInput( shortCastFilter->GetOutput() );
      shortWriter->SetImageIO( this->GetImageIO() );
      shortWriter->SetIORegion( this->GetIORegion() );
      shortWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      shortWriter->SetUseCompression( this->GetUseCompression() );
      shortWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      shortWriter->Update();
      break;
    case ImageIOBase::UINT:
      uintWriter = UIntWriterType::New();
      uintCastFilter = UIntCastFilterType::New();
      uintWriter->SetFileName( m_FileName );
      uintCastFilter->SetInput( inputPtr );
      uintWriter->SetInput( uintCastFilter->GetOutput() );
      uintWriter->SetImageIO( this->GetImageIO() );
      uintWriter->SetIORegion( this->GetIORegion() );
      uintWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      uintWriter->SetUseCompression( this->GetUseCompression() );
      uintWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      uintWriter->Update();
      break;
    case ImageIOBase::INT:
      intWriter = IntWriterType::New();
      intCastFilter = IntCastFilterType::New();
      intWriter->SetFileName( m_FileName );
      intCastFilter->SetInput( inputPtr );
      intWriter->SetInput( intCastFilter->GetOutput() );
      intWriter->SetImageIO( this->GetImageIO() );
      intWriter->SetIORegion( this->GetIORegion() );
      intWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      intWriter->SetUseCompression( this->GetUseCompression() );
      intWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      intWriter->Update();
      break;
    case ImageIOBase::ULONG:
      ulongWriter = ULongWriterType::New();
      ulongCastFilter = ULongCastFilterType::New();
      ulongWriter->SetFileName( m_FileName );
      ulongCastFilter->SetInput( inputPtr );
      ulongWriter->SetInput( ulongCastFilter->GetOutput() );
      ulongWriter->SetImageIO( this->GetImageIO() );
      ulongWriter->SetIORegion( this->GetIORegion() );
      ulongWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      ulongWriter->SetUseCompression( this->GetUseCompression() );
      ulongWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      ulongWriter->Update();
      break;
    case ImageIOBase::LONG:
      longWriter = LongWriterType::New();
      longCastFilter = LongCastFilterType::New();
      longWriter->SetFileName( m_FileName );
      longCastFilter->SetInput( inputPtr );
      longWriter->SetInput( longCastFilter->GetOutput() );
      longWriter->SetImageIO( this->GetImageIO() );
      longWriter->SetIORegion( this->GetIORegion() );
      longWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      longWriter->SetUseCompression( this->GetUseCompression() );
      longWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      longWriter->Update();
      break;
    case ImageIOBase::FLOAT:
      floatWriter = FloatWriterType::New();
      floatCastFilter = FloatCastFilterType::New();
      floatWriter->SetFileName( m_FileName );
      floatCastFilter->SetInput( inputPtr );
      floatWriter->SetInput( floatCastFilter->GetOutput() );
      floatWriter->SetImageIO( this->GetImageIO() );
      floatWriter->SetIORegion( this->GetIORegion() );
      floatWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      floatWriter->SetUseCompression( this->GetUseCompression() );
      floatWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      floatWriter->Update();
      break;
    case ImageIOBase::DOUBLE:
      doubleWriter = DoubleWriterType::New();
      doubleCastFilter = DoubleCastFilterType::New();
      doubleWriter->SetFileName( m_FileName );
      doubleCastFilter->SetInput( inputPtr );
      doubleWriter->SetInput( doubleCastFilter->GetOutput() );
      doubleWriter->SetImageIO( this->GetImageIO() );
      doubleWriter->SetIORegion( this->GetIORegion() );
      doubleWriter->SetNumberOfStreamDivisions( this->GetNumberOfStreamDivisions() );
      doubleWriter->SetUseCompression( this->GetUseCompression() );
      doubleWriter->SetUseInputMetaDataDictionary( this->GetUseInputMetaDataDictionary() );
      doubleWriter->Update();
      break;
    default:
      break;
    }

}


//---------------------------------------------------------
template <class TInputImage>
void 
CastVectorImageFileWriter<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (m_FileName.data() ? m_FileName.data() : "(none)") << std::endl;

  os << indent << "Image IO: ";
  if ( m_ImageIO.IsNull() )
    {
    os << "(none)\n";
    }
  else
    {
    os << m_ImageIO << "\n";
    }

  os << indent << "IO Region: " << m_PasteIORegion << "\n";
  os << indent << "Number of Stream Divisions: " << m_NumberOfStreamDivisions << "\n";

  if (m_UseCompression)
    {
    os << indent << "Compression: On\n";
    }
  else
    {
    os << indent << "Compression: Off\n";
    }

  if (m_UseInputMetaDataDictionary)
    {
    os << indent << "UseInputMetaDataDictionary: On\n";
    }
  else
    {
    os << indent << "UseInputMetaDataDictionary: Off\n";
    }

  if (m_FactorySpecifiedImageIO)
    {
    os << indent << "FactorySpecifiedmageIO: On\n";
    }
  else
    {
    os << indent << "FactorySpecifiedmageIO: Off\n";
    }

}

} // end namespace itk

#endif
