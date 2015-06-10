/**
 *       @file  itkDiffusionModelEstimationInSphericalCoordinateImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "03-10-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkDiffusionModelEstimationInSphericalCoordinateImageFilter_h
#define __itkDiffusionModelEstimationInSphericalCoordinateImageFilter_h


#include <tr1/memory>
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkDiffusionModelEstimationImageFilter.h"

#include "vnl/vnl_matrix.h"
#include "utlITK.h"

namespace itk
{

/**
 *   \class   DiffusionModelEstimationInSphericalCoordinateImageFilter
 *   \brief   base filter for estimation of diffusion models in spherical coordinates
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT DiffusionModelEstimationInSphericalCoordinateImageFilter
  : public DiffusionModelEstimationImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DiffusionModelEstimationInSphericalCoordinateImageFilter         Self;
  typedef DiffusionModelEstimationImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( DiffusionModelEstimationInSphericalCoordinateImageFilter, DiffusionModelEstimationImageFilter );
  
  /** Convenient Typedefs. */
  typedef typename Superclass::InputImageType           InputImageType;
  typedef typename Superclass::InputImagePointer        InputImagePointer;
  typedef typename Superclass::InputImageConstPointer   InputImageConstPointer;
  typedef typename Superclass::InputImageIndexType      InputImageIndexType;
  typedef typename Superclass::InputImageSizeType       InputImageSizeType;
  typedef typename Superclass::InputImageSpacingType    InputImageSpacingType;
  typedef typename Superclass::InputImagePixelType      InputImagePixelType;
  typedef typename Superclass::InputImageRegionType     InputImageRegionType;
  
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::OutputImageIndexType     OutputImageIndexType;
  typedef typename Superclass::OutputImageSizeType      OutputImageSizeType;
  typedef typename Superclass::OutputImageSpacingType   OutputImageSpacingType;
  typedef typename Superclass::OutputImagePixelType     OutputImagePixelType;
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
  
  typedef typename Superclass::MaskImageType            MaskImageType;
  
  typedef typename Superclass::MatrixType         MatrixType;
  typedef typename Superclass::VectorType         VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType      STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;
  
  itkSetMacro(SHRank, int);
  itkGetMacro(SHRank, int);
  itkSetMacro(RadialRank, int);
  itkGetMacro(RadialRank, int);
  
  // void SetBasisSHMatrix(const MatrixPointer& mat);
  itkGetMacro(BasisSHMatrix, MatrixPointer);
  // void SetBasisRadialMatrix(const MatrixPointer& mat);
  itkGetMacro(BasisRadialMatrix, MatrixPointer);

  /** compute SH matrix, grad is in spherical format  */
  void ComputeSHMatrix();

  /** compute radial matrix  */
  virtual void ComputeRadialMatrix () 
    {}

protected:
  DiffusionModelEstimationInSphericalCoordinateImageFilter(); 

  virtual ~DiffusionModelEstimationInSphericalCoordinateImageFilter() {};
  
  virtual void VerifyInputParameters() const;

  void PrintSelf(std::ostream& os, Indent indent) const;
  typename LightObject::Pointer InternalClone() const;

  /** rank for spherical part (SH basis)  */
  int m_SHRank;
  /** rank for radial part  */
  int m_RadialRank;
  
  MatrixPointer  m_BasisSHMatrix;
  MatrixPointer  m_BasisRadialMatrix;

private:
  DiffusionModelEstimationInSphericalCoordinateImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

  
};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDiffusionModelEstimationInSphericalCoordinateImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDiffusionModelEstimationInSphericalCoordinateImageFilter_hxx)
#include "itkDiffusionModelEstimationInSphericalCoordinateImageFilter.hxx"
#endif


#endif 

