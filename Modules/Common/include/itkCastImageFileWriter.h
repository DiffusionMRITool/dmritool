/*=========================================================================
 
 Program:   Cast Image File Writer
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef __itkCastImageFileWriter_h
#define __itkCastImageFileWriter_h

#include "itkProcessObject.h"
#include "itkImageIOBase.h"
#include "itkExceptionObject.h"
#include "itkSize.h"
#include "itkImageIORegion.h"

namespace itk
{

/** \class CastImageFileWriter
 * \brief Writes image data, after casting, to a single file.
 *
 * \ingroup IOFilters 
 */
template <class TInputImage>
class ITK_EXPORT CastImageFileWriter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef CastImageFileWriter       Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CastImageFileWriter,ProcessObject);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename InputImageType::Pointer       InputImagePointer;
  typedef typename InputImageType::RegionType    InputImageRegionType; 
  typedef typename InputImageType::PixelType     InputImagePixelType;
  
  typedef typename ImageIOBase::IOComponentType  IOComponentType;
  
  /** Set/Get the image input of this writer.  */
  using Superclass::SetInput;
  virtual void SetInput(const InputImageType *input);
  const InputImageType * GetInput(void);
  const InputImageType * GetInput(unsigned int idx);
  
  /** Specify the name of the output file to write. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);
  
  /** Set/Get the ImageIO helper class. Usually this is created via the object
   * factory mechanism that determines whether a particular ImageIO can
   * write a certain file. This method provides a way to get the ImageIO 
   * instance that is created, or one can be manually set where the
   * IO factory mechanism may not work (for example, raw image files or
   * image files with non-standard filename suffix's.
   * If the user specifies the ImageIO, we assume she makes the
   * correct choice and will allow a file to be created regardless of
   * the file extension. If the factory has set the ImageIO, the
   * extension must be supported by the specified ImageIO. */
  void SetImageIO (ImageIOBase* io)
    {
    if (this->m_ImageIO != io)
      {
      this->Modified();
      this->m_ImageIO = io;
      }
    m_FactorySpecifiedImageIO = false;
    }
  itkGetObjectMacro(ImageIO,ImageIOBase);
  
  /** A special version of the Update() method for writers.  It
   * invokes start and end events and handles releasing data. It
   * eventually calls GenerateData() which does the actual writing.
   * Note: the write method will write data specified by the
   * IORegion. If not set, then then the whole image is written.  Note
   * that the region will be cropped to fit the input image's
   * LargestPossibleRegion. */
  virtual void Write(void);

  /** Specify the region to write. If left NULL, then the whole image
   * is written. */
  void SetIORegion(const ImageIORegion & region);
  const ImageIORegion &GetIORegion(void) const 
    {
    return m_PasteIORegion;
    }


  /** Set/Get the number of pieces to divide the input.  The upstream pipeline
   * will try to be executed this many times. */
  itkSetMacro(NumberOfStreamDivisions,unsigned int);
  itkGetConstReferenceMacro(NumberOfStreamDivisions,unsigned int);


  /** Aliased to the Write() method to be consistent with the rest of the
   * pipeline. */
  virtual void Update() ITK_OVERRIDE
    {
    this->Write();
    }

  /** Set the compression On or Off */
  itkSetMacro(UseCompression,bool);
  itkGetConstReferenceMacro(UseCompression,bool);
  itkBooleanMacro(UseCompression);

  /** By default the MetaDataDictionary is taken from the input image and 
   *  passed to the ImageIO. In some cases, however, a user may prefer to 
   *  introduce her/his own MetaDataDictionary. This is often the case of
   *  the ImageSeriesWriter. This flag defined whether the MetaDataDictionary 
   *  to use will be the one from the input image or the one already set in
   *  the ImageIO object. */
  itkSetMacro(UseInputMetaDataDictionary,bool);
  itkGetConstReferenceMacro(UseInputMetaDataDictionary,bool);
  itkBooleanMacro(UseInputMetaDataDictionary);

  /** Set/Get component type */
  itkSetEnumMacro(ComponentType, IOComponentType);
  itkGetEnumMacro(ComponentType, IOComponentType);
  
protected:
  CastImageFileWriter();
  ~CastImageFileWriter();
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  /** Does the real work. */
  void GenerateData(void) ITK_OVERRIDE;
  
private:
  CastImageFileWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string        m_FileName;
  
  ImageIOBase::Pointer m_ImageIO;
  bool                 m_UserSpecifiedImageIO; // track whether the ImageIO
                                                // is user specified
  
  ImageIORegion m_PasteIORegion;
  unsigned int  m_NumberOfStreamDivisions;
  bool          m_UserSpecifiedIORegion;    // track whether the region
                                            // is user specified
  bool          m_FactorySpecifiedImageIO;  //track whether the factory
                                            //  mechanism set the ImageIO
  bool m_UseCompression;
  bool m_UseInputMetaDataDictionary;        // whether to use the
                                            // MetaDataDictionary from the
                                            // input or not.
                                            
  IOComponentType m_ComponentType;
};

  
} // end namespace itk
  
#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkCastImageFileWriter_hxx)
#include "itkCastImageFileWriter.hxx"
#endif

#endif // __itkCastImageFileWriter_h
