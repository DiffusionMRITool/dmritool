/**
 *       @file  itkDWIWriter.h
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

#ifndef __itkDWIWriter_h
#define __itkDWIWriter_h



#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageSource.h"
#include "vnl/vnl_matrix.h"
#include "utlITK.h"
#include <tr1/memory>

#include "itkSamplingSchemeQSpace.h"

namespace itk
{
/**
 *   \class   DWIWriter
 *   \brief   write gradient file, b values and DWI files (with optional index file)
 *
 *   The data format has two types: 
 *
 * ******in data.txt******
 * b1.txt  grad1.txt  dwi1.nii.gz 
 * ***********************
 * All shells are stored in one file.
 * Or
 * ******in data.txt******
 * b1.txt  grad1.txt  dwi1.nii.gz 
 * b2.txt  grad2.txt  dwi2.nii.gz
 * b3.txt  grad3.txt  dwi3.nii.gz
 * ***********************
 * Each line is for a single shell data if it is set to store each shell. 
 *
 * dwi.hdr can be stored as VectorImage 
 * b.txt and grad.txt, dwi.nii.gz should have the same dimension
 *
 *
 *   \ingroup DWIProcessing
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template <class TPixelType, unsigned int VImageDimension = 3 >
class ITK_EXPORT DWIWriter : public ImageSource< VectorImage<TPixelType, VImageDimension> >
{

public:
  /** Standard class typedefs. */
  typedef DWIWriter             Self;
  typedef ImageSource< VectorImage<TPixelType, VImageDimension> > Superclass;
  typedef SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DWIWriter, ImageSource);
  
  typedef Image<TPixelType, VImageDimension>       B0ImageType;
  typedef VectorImage<TPixelType, VImageDimension> DWIImageType;
  typedef B0ImageType                              MaskImageType;
  
  typedef utl::NDArray<double,2>                     MatrixType;
  typedef utl_shared_ptr<MatrixType>                 MatrixPointer;
  typedef std::vector<double>                        STDVectorType;
  typedef utl_shared_ptr<STDVectorType >             STDVectorPointer;
  
  typedef typename DWIImageType::SizeType          SizeType;
  typedef typename DWIImageType::IndexType         IndexType;
  typedef typename DWIImageType::RegionType        RegionType;
  typedef typename DWIImageType::PixelType         PixelType;
  typedef typename DWIImageType::DirectionType     DirectionType;
  typedef typename DWIImageType::SpacingType       SpacingType;
  
  typedef typename B0ImageType::SizeType          B0SizeType;
  typedef typename B0ImageType::IndexType         B0IndexType;
  typedef typename B0ImageType::RegionType        B0RegionType;
  typedef typename B0ImageType::PixelType         B0PixelType;
  typedef typename B0ImageType::DirectionType     B0DirectionType;
  typedef typename B0ImageType::SpacingType       B0SpacingType;
  
  typedef SamplingSchemeQSpace<double>                     SamplingSchemeQSpaceType;
  typedef typename SamplingSchemeQSpaceType::Pointer       SamplingSchemeQSpacePointer;
  
  itkSetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkGetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  
  /** Set/Get the mask image. */
  itkSetObjectMacro(MaskImage, MaskImageType);
  itkGetConstObjectMacro(MaskImage, MaskImageType);
  
  /** Set/Get the mask image. */
  itkSetObjectMacro(B0Image, B0ImageType);
  itkGetObjectMacro(B0Image, B0ImageType);

  itkSetMacro(ConfigurationFile, std::string);
  itkGetMacro(ConfigurationFile, std::string);
  itkSetMacro(DWIFile, std::string);
  itkGetMacro(DWIFile, std::string);
  itkSetMacro(OrientationFile, std::string);
  itkGetMacro(OrientationFile, std::string);
  itkSetMacro(BFile, std::string);
  itkGetMacro(BFile, std::string);
  
  /** If set true, relative paths are used in configuration file  */
  itkSetMacro(UseRelativePath, bool);
  itkGetMacro(UseRelativePath, bool);
  itkBooleanMacro(UseRelativePath);
  
  /** Output each shell  */
  itkSetMacro(OutputEachShell, bool);
  itkGetMacro(OutputEachShell, bool);
  itkBooleanMacro(OutputEachShell);
  
  /** Set/Get the image input of this writer.  */
  using Superclass::SetInput;
  void SetInput(const DWIImageType *input);
  const DWIImageType * GetInput();
  
  /** 
   *   The data format has two types: 
   *
   * ******in data.txt******
   * b1.txt  grad1.txt  dwi1.nii.gz 
   * ***********************
   * All shells are stored in one file.
   * Or
   * ******in data.txt******
   * b1.txt  grad1.txt  dwi1.nii.gz 
   * b2.txt  grad2.txt  dwi2.nii.gz
   * b3.txt  grad3.txt  dwi3.nii.gz
   * ***********************
   * Each line is for a single shell data if it is set to store each shell. 
   * */
  void WriteToConfigurationFile(const std::string file);
  
  
protected:
  DWIWriter();
  virtual ~DWIWriter(){}
 
  void PrintSelf(std::ostream& os, Indent indent) const;
 
  /** Does the real work. */
  virtual void GenerateData();
 
private:
  DWIWriter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  
  SamplingSchemeQSpacePointer m_SamplingSchemeQSpace;
  
  bool m_UseRelativePath;
  bool m_OutputEachShell;

  std::string m_ConfigurationFile;
  std::string m_DWIFile;
  std::string m_BFile;
  std::string m_OrientationFile;

  typename MaskImageType::Pointer m_MaskImage;
  typename B0ImageType::Pointer m_B0Image;
};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDWIWriter_hxx)
#include "itkDWIWriter.hxx"
#endif

#endif 
