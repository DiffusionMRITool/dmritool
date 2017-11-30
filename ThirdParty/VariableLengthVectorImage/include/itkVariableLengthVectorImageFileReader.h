/*=========================================================================
 
 Program:   Variable Length Vector Image File Reader
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef __itkVariableLengthVectorImageFileReader_h
#define __itkVariableLengthVectorImageFileReader_h

#include "itkExceptionObject.h"
#include "itkImage.h"
#include "itkImageSource.h"
#include "itkImageIOBase.h"
#include "itkImageFileReader.h"
#include "itksys/SystemTools.hxx"

namespace itk
{

/** \class VariableLengthVectorImageFileReaderException
 *
 * \brief Base exception class for IO conflicts.
 *
 */
class VariableLengthVectorImageFileReaderException : public ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro( VariableLengthVectorImageFileReaderException, ExceptionObject );

  /** Constructor. */
  VariableLengthVectorImageFileReaderException(const char *file, unsigned int line,
                           const char* message = "Error in IO",
                           const char* loc = "Unknown") :
    ExceptionObject(file, line, message, loc)
  {
  }

  /** Constructor. */
  VariableLengthVectorImageFileReaderException(const std::string& file, unsigned int line,
                           const char* message = "Error in IO",
                           const char* loc = "Unknown") : 
    ExceptionObject(file, line, message, loc)
  {
  }
};


/** \class VariableLengthVectorImageFileReader
 * \brief Reads variable length vector image data.
 *
 */
template <class TOutputImage>
class ITKIOImageBase_HIDDEN VariableLengthVectorImageFileReader : public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef VariableLengthVectorImageFileReader   Self;
  typedef ImageSource<TOutputImage>             Superclass;
  typedef SmartPointer<Self>                    Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VariableLengthVectorImageFileReader, ImageSource);

  /** The dimensionality of the input and output images. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::OutputImageDimension);

  /** Some convenient typedefs. */
  typedef TOutputImage  OutputImageType;
  typedef typename OutputImageType::Pointer                       OutputImagePointer;
  typedef typename OutputImageType::SizeType                      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType                   OutputImageSpacingType;
  typedef typename OutputImageType::IndexType                     OutputImageIndexType;
  typedef typename OutputImageType::PointType                     OutputImagePointType;
  typedef typename OutputImageType::DirectionType                 OutputImageDirectionType;
  typedef typename OutputImageType::RegionType                    OutputImageRegionType;
  typedef typename OutputImageType::PixelType                     OutputImagePixelType;
  typedef typename OutputImageType::PixelContainer                OutputImagePixelContainerType;
  typedef typename OutputImageType::PixelContainer::Pointer       OutputImagePixelContainerPointerType;  
//  typedef typename OutputImageType::PixelContainer::PixelMapType  OutputImagePixelMapType;
//  typedef typename OutputImagePixelMapType::iterator              OutputImagePixelMapIteratorType;
  
  /** Specify the files to read. This is forwarded to the IO instance. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);
  
  /** Set the stream On or Off */
  itkSetMacro(UseStreaming,bool);
  itkGetConstReferenceMacro(UseStreaming,bool);
  itkBooleanMacro(UseStreaming);
  
  /** Set/Get the ImageIO helper class */
  itkGetObjectMacro(ImageIO,ImageIOBase);

protected:
  VariableLengthVectorImageFileReader();
  ~VariableLengthVectorImageFileReader();
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  
  /** Convert a block of pixels from one type to another. */
//  void DoConvertBuffer(void* buffer, size_t numberOfPixels);

  /** Test whether the given filename exist and it is readable, this
    * is intended to be called before attempting to use  ImageIO
    * classes for actually reading the file. If the file doesn't exist
    * or it is not readable, and exception with an approriate message
    * will be thrown. */
//  void TestFileExistanceAndReadability();

  typedef unsigned int LengthType;
  typedef Image<LengthType, ImageDimension> LengthImageType;
  typedef Image<typename OutputImagePixelType::ValueType, 1> ValueImageType;
      
  typename LengthImageType::Pointer m_LengthImage;
  typename ValueImageType::Pointer m_ValueImage;

  typedef itk::ImageFileReader<LengthImageType> LengthImageFileReaderType;
  typedef itk::ImageFileReader<ValueImageType> ValueImageFileReaderType;
  
  typename LengthImageFileReaderType::Pointer m_LengthImageFileReader;
  typename ValueImageFileReaderType::Pointer m_ValueImageFileReader;

  /** Does the real work. */
  virtual void GenerateData() ITK_OVERRIDE;

  std::string m_FileName;
  bool m_UseStreaming;
  ImageIOBase::Pointer m_ImageIO;
  
private:
  VariableLengthVectorImageFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  std::string m_ExceptionMessage;  

  // The region that the ImageIO class will return when we ask to
  // produce the requested region.
//  ImageIORegion m_ActualIORegion; 
};


} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVariableLengthVectorImageFileReader.hxx"
#endif

#endif // __itkSparseVariableLengthVectorImageFileReader_h
