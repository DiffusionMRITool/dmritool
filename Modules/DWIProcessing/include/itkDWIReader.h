/**
 *       @file  itkDWIReader.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-11-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkDWIReader_h
#define __itkDWIReader_h

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageSource.h"
#include "utlITK.h"
#include <tr1/memory>

#include "itkSamplingSchemeQSpace.h"

namespace itk
{
/**
 *   \class   DWIReader
 *   \brief   Load gradient file, b values and DWI files (with optional index file)
 *
 *   The data format has two types: 
 *
 * ******in data.txt*******************
 * 650   grad1.txt  dwi1.nii index1.txt
 * 1500  grad2.txt  dwi2.nii index2.txt
 * 3000  grad3.txt  dwi3.nii index3.txt
 * ************************************
 * Each line is for a single shell data. 
 * If the dimension N in grad.txt is smaller than the dimension M in dwi.nii, the first M-N dimension in dwi.hpp is for b0 image
 * index.txt (optional) shows the index requested for reading.  
 * dwi.hdr can be stored as VectorImage or Image 
 *
 * ******in data.txt*********************
 * b1.txt  grad1.txt  dwi1.nii index1.txt
 * b2.txt  grad2.txt  dwi2.nii index2.txt
 * b3.txt  grad3.txt  dwi3.nii index3.txt
 * **************************************
 * b.txt and grad.txt should have the same dimension
 *
 *
 *   \ingroup DWIProcessing
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template <class TPixelType, unsigned int VImageDimension = 3 >
class ITK_EXPORT DWIReader : public ImageSource< VectorImage<TPixelType, VImageDimension> >
{

public:
  /** Standard class typedefs. */
  typedef DWIReader             Self;
  typedef ImageSource< VectorImage<TPixelType, VImageDimension> > Superclass;
  typedef SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DWIReader, ImageSource);
  
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
  
  typedef typename B0ImageType::SizeType           B0SizeType;
  typedef typename B0ImageType::IndexType          B0IndexType;
  typedef typename B0ImageType::RegionType         B0RegionType;
  typedef typename B0ImageType::PixelType          B0PixelType;
  typedef typename B0ImageType::DirectionType      B0DirectionType;
  typedef typename B0ImageType::SpacingType        B0SpacingType;
  
  typedef SamplingSchemeQSpace<double>                         SamplingSchemeQSpaceType;
  typedef typename SamplingSchemeQSpaceType::Pointer           SamplingSchemeQSpacePointer;

  /** Set/Get the mask image. */
  itkSetObjectMacro(MaskImage, MaskImageType);
  itkGetConstObjectMacro(MaskImage, MaskImageType);
  
  /** Set/Get the mask image. */
  itkSetObjectMacro(B0Image, B0ImageType);
  itkGetObjectMacro(B0Image, B0ImageType);
  
  itkSetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkGetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);

  typename DWIImageType::Pointer GetDWIImage()
    {
    this->GetOutput();
    }
  
  itkSetMacro(ConfigurationFile, std::string);
  itkGetMacro(ConfigurationFile, std::string);
  
  /** If set on, output DWI Image will be normalized using B0Image  */
  itkSetMacro(NormalizeDWI, bool);
  itkGetMacro(NormalizeDWI, bool);
  itkBooleanMacro(NormalizeDWI);
  
  // [>* If set on, the first dimension of output DWI Image will be B0Image  <]
  // itkSetMacro(DWIWithB0, bool);
  // itkGetMacro(DWIWithB0, bool);
  // itkBooleanMacro(DWIWithB0);
  
  /** If set on, convert 4D DWI image into 3D VectorImage */
  itkSetMacro(IsInput4DImage, bool);
  itkGetMacro(IsInput4DImage, bool);
  itkBooleanMacro(IsInput4DImage);
  
  /** If set on, correct values  */
  itkSetMacro(CorrectDWIValues, bool);
  itkGetMacro(CorrectDWIValues, bool);
  itkBooleanMacro(CorrectDWIValues);
  
  /** if it is true, print some warnings if the dwi values are abnormal (i.e. zero, negative, larger than b0 value)  */
  itkSetMacro(ShowWarnings, bool);
  itkGetMacro(ShowWarnings, bool);
  itkBooleanMacro(ShowWarnings);
  
  /** 
   *   The data format has two types: 
   *
   * ******in data.txt*******************
   * 650   grad1.txt  dwi1.nii index1.txt
   * 1500  grad2.txt  dwi2.nii index2.txt
   * 3000  grad3.txt  dwi3.nii index3.txt
   * ************************************
   * Each line is for a single shell data. 
   * If the dimension N in grad.txt is smaller than the dimension M in dwi.nii, the first M-N dimension in dwi.hpp is for b0 image
   * index.txt (optional) shows the index requested for reading.  
   * dwi.hdr can be stored as VectorImage or Image 
   *
   * ******in data.txt*********************
   * b1.txt  grad1.txt  dwi1.nii index1.txt
   * b2.txt  grad2.txt  dwi2.nii index2.txt
   * b3.txt  grad3.txt  dwi3.nii index3.txt
   * **************************************
   * b.txt and grad.txt should have the same dimension
   * */
  void ReadFromConfigurationFile(const std::string file);
  
  // void ReadBGradDWI(const std::string bFile, const std::string gradFile, const std::string imageFile, const std::string indexFile="");
  
  static bool DetermineIsInput4DImage(const std::string dataStr);
  
  /** normalize DWI images using b0 image  */
  void NormalizeDWI();

  /** correct DWI values if they are negative or they are larger than values in b0 image  */
  void CorrectDWI();
  
protected:
  DWIReader();
  virtual ~DWIReader(){}
 
  void PrintSelf(std::ostream& os, Indent indent) const;
 
  /** Does the real work. */
  virtual void GenerateData();
 
private:
  DWIReader(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  SamplingSchemeQSpacePointer m_SamplingSchemeQSpace;

  bool m_NormalizeDWI;
  // bool m_DWIWithB0;
  bool m_IsInput4DImage;
  bool m_CorrectDWIValues;

  std::string m_ConfigurationFile;
  
  bool m_ShowWarnings;

  typename MaskImageType::Pointer m_MaskImage;
  typename B0ImageType::Pointer m_B0Image;
};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDWIReader_hxx)
#include "itkDWIReader.hxx"
#endif

#endif 

