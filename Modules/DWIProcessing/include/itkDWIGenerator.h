/*=========================================================================

 Program:   DWI Generator

 Copyright (c) Pew-Thian Yap, Jian Cheng. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
 
#ifndef __itkDWIGenerator_h
#define __itkDWIGenerator_h


#include "itkDWIGeneratorBase.h"

namespace itk
{
  
/** \class DWIGenerator
 *  \brief Generate DWI data based on provided parameter file.
 *
 */
template <class TOutputImage, class TScalarImage=Image<float,3> >
class ITK_EXPORT DWIGenerator : public DWIGeneratorBase<TOutputImage, TScalarImage>
{
public:
  /** Standard class typedefs. */
  typedef DWIGenerator                Self;
  typedef DWIGeneratorBase<TOutputImage, TScalarImage>   Superclass;
  typedef SmartPointer<Self>          Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DWIGenerator, DWIGeneratorBase);

  /** Some convenient typedefs. */
  typedef typename Superclass::OutputImageType                  OutputImageType;
  typedef typename Superclass::OutputImagePointer               OutputImagePointer;
  typedef typename Superclass::OutputImageSizeType              OutputImageSizeType;
  typedef typename Superclass::OutputImageSpacingType           OutputImageSpacingType;
  typedef typename Superclass::OutputImageIndexType             OutputImageIndexType;
  typedef typename Superclass::OutputImagePointType             OutputImagePointType;
  typedef typename Superclass::OutputImageDirectionType         OutputImageDirectionType;
  typedef typename Superclass::OutputImageRegionType            OutputImageRegionType;
  typedef typename Superclass::OutputImagePixelType             OutputImagePixelType;
  typedef typename Superclass::OutputImageInternalPixelType     OutputImageInternalPixelType;

  /** Output image dimension */
  itkStaticConstMacro(OutputImageDimension, unsigned int, OutputImageType::ImageDimension);

  /** B0 Image */
  typedef typename Superclass::ScalarImageType           ScalarImageType;
  typedef typename Superclass::ScalarImagePointer        ScalarImagePointer;

  /** Orientation Matrice Type */
  typedef typename Superclass::PrecisionType         PrecisionType;
  typedef typename Superclass::MatrixType            MatrixType;
  typedef typename Superclass::MatrixPointer         MatrixPointer;
  typedef typename Superclass::STDVectorType         STDVectorType;
  typedef typename Superclass::STDVectorPointer      STDVectorPointer;
    
  /** Some convenient typedefs for diffusion parameters. */  
  typedef typename Superclass::DiffusionParameterValuesType    DiffusionParameterValuesType;
  typedef typename Superclass::DiffusionParameterContainerType DiffusionParameterContainerType;

  typedef typename Superclass::CylinderModelType      CylinderModelType;
  typedef typename Superclass::CylinderModelPointer   CylinderModelPointer;
  
  itkSetNDebugMacro( BackgroundDiffusionParameterValues, STDVectorType );
  itkGetMacro( BackgroundDiffusionParameterValues, STDVectorType );
    

  /** Specify the files to read. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);
  
protected:
  DWIGenerator();
  ~DWIGenerator();

  typename LightObject::Pointer InternalClone() const;
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Does the real work. */
  void GenerateData();
  
  std::string m_FileName;
  
  STDVectorType m_BackgroundDiffusionParameterValues;

private:
  DWIGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
};


} //namespace ITK

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDWIGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDWIGenerator_hxx)
#include "itkDWIGenerator.hxx"
#endif

#endif // __itkSparseDWIGenerator_h
