/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image File Reader

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImageFileReader_h
#define __itkSpatiallyDenseSparseVectorImageFileReader_h

#include "itkExceptionObject.h"
#include "itkImage.h"
#include "itkImageSource.h"
#include "itkImageIOBase.h"
#include "itkImageFileReader.h"
#include "itksys/SystemTools.hxx"
#include "itkMacro.h"


namespace itk
{

/** \class ImageFileReaderException
 *
 * \brief Base exception class for IO conflicts.
 *
 * \ingroup ITKSparseVectorImage
 */
class SpatiallyDenseSparseVectorImageFileReaderException : public ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro( SpatiallyDenseSparseVectorImageFileReaderException, ExceptionObject );

  /** Constructor. */
  SpatiallyDenseSparseVectorImageFileReaderException(const char *file, unsigned int line,
                           const char* message = "Error in IO",
                           const char* loc = "Unknown") :
    ExceptionObject(file, line, message, loc)
  {
  }

  /** Constructor. */
  SpatiallyDenseSparseVectorImageFileReaderException(const std::string &file, unsigned int line,
                           const char* message = "Error in IO",
                           const char* loc = "Unknown") :
    ExceptionObject(file, line, message, loc)
  {
  }
};


/** \class SpatiallyDenseSparseVectorImageFileReader
 * \brief Reads sparse image data from key and value files.
 *
 * \ingroup ITKSpatiallyDenseSparseVectorImage
 *
 * \author Pew-Thian Yap (ptyap@med.unc.edu)
 */
template <class TOutputImage>
class ITKIOImageBase_HIDDEN SpatiallyDenseSparseVectorImageFileReader : public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SpatiallyDenseSparseVectorImageFileReader   Self;
  typedef ImageSource<TOutputImage>                   Superclass;
  typedef SmartPointer<Self>                          Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SpatiallyDenseSparseVectorImageFileReader, ImageSource);

  /** Dimension of image. */
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
  typedef typename OutputImageType::InternalPixelType             OutputImageInternalPixelType;
  typedef typename OutputImageType::ValueType                     OutputImageValueType;
  typedef typename OutputImageType::KeyType                       OutputImageKeyType;
  typedef typename OutputImageType::PixelContainer                OutputImagePixelContainerType;
  typedef typename OutputImageType::PixelContainerConstPointer    OutputImagePixelContainerConstPointerType;
  typedef typename OutputImageType::PixelMapType                  OutputImagePixelPixelMapType;
  typedef typename OutputImageType::PixelMapConstIterator         OutputImagePixelMapConstIteratorType;

  typedef Image<OutputImageKeyType, 1> KeyImageType;
  typedef Image<OutputImageValueType, 1> ValueImageType;

  typedef ImageFileReader<KeyImageType> KeyImageFileReaderType;
  typedef ImageFileReader<ValueImageType> ValueImageFileReaderType;

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
  SpatiallyDenseSparseVectorImageFileReader();
  ~SpatiallyDenseSparseVectorImageFileReader();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Convert a block of pixels from one type to another. */
//  void DoConvertBuffer(void* buffer, size_t numberOfPixels);

  /** Test whether the given filename exist and it is readable, this
    * is intended to be called before attempting to use  ImageIO
    * classes for actually reading the file. If the file doesn't exist
    * or it is not readable, and exception with an approriate message
    * will be thrown. */
//  void TestFileExistanceAndReadability();

  typename KeyImageType::Pointer m_KeyImage;
  typename ValueImageType::Pointer m_ValueImage;

  typename KeyImageFileReaderType::Pointer m_KeyImageFileReader;
  typename ValueImageFileReaderType::Pointer m_ValueImageFileReader;

  /** Does the real work. */
  virtual void GenerateData();

  std::string m_FileName;
  bool m_UseStreaming;
  ImageIOBase::Pointer m_ImageIO;

private:
  SpatiallyDenseSparseVectorImageFileReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string m_ExceptionMessage;

  // The region that the ImageIO class will return when we ask to
  // produce the requested region.
//  ImageIORegion m_ActualIORegion;
};


} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatiallyDenseSparseVectorImageFileReader.hxx"
#endif

#endif // __itkSparseSpatiallyDenseSparseVectorImageFileReader_h
