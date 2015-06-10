/**
 *       @file  itkDiffusionModelEstimationImageFilter.h
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

#ifndef __itkDiffusionModelEstimationImageFilter_h
#define __itkDiffusionModelEstimationImageFilter_h


#include <tr1/memory>
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkMaskedImageToImageFilter.h"

#include "vnl/vnl_matrix.h"
#include "utlITK.h"
#include "utlMatrix.h"

#include "itkSamplingSchemeQSpace.h"

namespace itk
{

/**
 *   \class   DiffusionModelEstimationImageFilter
 *   \brief   base filter for estimation of diffusion models
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT DiffusionModelEstimationImageFilter
  : public MaskedImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DiffusionModelEstimationImageFilter         Self;
  typedef MaskedImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( DiffusionModelEstimationImageFilter, MaskedImageToImageFilter );
  
  /** Convenient Typedefs. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::IndexType      InputImageIndexType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SpacingType    InputImageSpacingType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  
  typedef TOutputImage                            OutputImageType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename OutputImageType::IndexType     OutputImageIndexType;
  typedef typename OutputImageType::SizeType      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType   OutputImageSpacingType;
  typedef typename OutputImageType::PixelType     OutputImagePixelType;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  
  typedef typename Superclass::MaskImageType      MaskImageType;

  typedef utl::NDArray<double,2>                  MatrixType;
  typedef utl::NDArray<double,1>                  VectorType;
  typedef utl_shared_ptr<MatrixType>              MatrixPointer;
  typedef utl_shared_ptr<VectorType>              VectorPointer;
  typedef std::vector<double>                     STDVectorType;
  typedef utl_shared_ptr<STDVectorType >          STDVectorPointer;
  
  typedef SamplingSchemeQSpace<double>                       SamplingSchemeQSpaceType;
  typedef typename SamplingSchemeQSpaceType::Pointer         SamplingSchemeQSpacePointer;
  
  typedef SamplingScheme3D<double>                           SamplingSchemeRSpaceType;
  typedef typename SamplingSchemeRSpaceType::Pointer         SamplingSchemeRSpacePointer;
  
  itkSetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkGetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);

  // itkSetMacro(BasisMatrix, MatrixPointer);
  itkGetMacro(BasisMatrix, MatrixPointer);
  
  itkSetMacro(MD0, double);
  itkGetMacro(MD0, double);

  /** compute basis matrix which incorporates radial matrix and spherical matrix  */
  virtual void ComputeBasisMatrix ()
    {}
  
  virtual void ComputeRegularizationWeight ( )
    {}

protected:
  DiffusionModelEstimationImageFilter(); 

  virtual ~DiffusionModelEstimationImageFilter() {};
  
  virtual void VerifyInputParameters() const;

  void PrintSelf(std::ostream& os, Indent indent) const;
  typename LightObject::Pointer InternalClone() const;
  
  /** Sampling Scheme in q-space  */
  SamplingSchemeQSpacePointer m_SamplingSchemeQSpace;
  
  /** typical MD value for typical scale */
  double m_MD0;

  MatrixPointer  m_BasisMatrix;

  /** used for regulaization  */
  VectorPointer  m_RegularizationWeight;

private:
  DiffusionModelEstimationImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

  
};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDiffusionModelEstimationImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDiffusionModelEstimationImageFilter_hxx)
#include "itkDiffusionModelEstimationImageFilter.hxx"
#endif


#endif 

