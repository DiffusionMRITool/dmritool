/**
 *       @file  itkODFFromSPFImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-08-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkODFFromSPFImageFilter_h
#define __itkODFFromSPFImageFilter_h

#include "itkFeaturesFromSPFImageFilter.h"


namespace itk
{  

/**
 *   \class   ODFFromSPFImageFilter
 *   \brief   calculate ODFs from SPF coefficients
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT ODFFromSPFImageFilter :
public FeaturesFromSPFImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ODFFromSPFImageFilter         Self;
  typedef FeaturesFromSPFImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( ODFFromSPFImageFilter, FeaturesFromSPFImageFilter );
  
  typedef TInputImage                                      InputImageType;
  typedef typename InputImageType::Pointer                 InputImagePointer;

  typedef typename Superclass::ScalarImageType             ScalarImageType;
  typedef typename Superclass::ScalarImagePointer          ScalarImagePointer;
  typedef typename Superclass::MatrixType                  MatrixType;
  typedef typename Superclass::MatrixPointer               MatrixPointer;
  typedef typename Superclass::VectorType                  VectorType;
  typedef typename Superclass::VectorPointer               VectorPointer;
  typedef typename Superclass::BasisType                   BasisType;
  typedef typename Superclass::STDVectorType               STDVectorType;
  typedef typename Superclass::STDVectorPointer            STDVectorPointer;
  
  typedef typename Superclass::MaskImageType               MaskImageType;
  
  /** for odf  */
  itkSetMacro(ODFOrder, int);
  itkGetMacro(ODFOrder, int);
  itkSetMacro(BMax, double);
  itkGetMacro(BMax, double);
  
  void ComputeSPFToFeatureTransform() ITK_OVERRIDE;
  

protected:
  ODFFromSPFImageFilter() : Superclass()
  {
  m_BMax=-1;
  m_ODFOrder=2;
  }

  void VerifyInputParameters() const ITK_OVERRIDE;

  virtual ~ODFFromSPFImageFilter() {};
  
  void GenerateOutputInformation() ITK_OVERRIDE;
  
  void BeforeThreadedGenerateData () ITK_OVERRIDE;
  void ThreadedGenerateData(const typename TOutputImage::RegionType& outputRegionForThread,ThreadIdType threadId ) ITK_OVERRIDE;
  
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

private:
  ODFFromSPFImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented
  
  VectorType m_P;
  VectorType m_L;
  
  /** for ODF  */
  int m_ODFOrder;
  /** maximal b value for the disk interal for ODF estimation  */
  double m_BMax;
};



} // end namespace itk


#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkODFFromSPFImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkODFFromSPFImageFilter_hxx)
#include "itkODFFromSPFImageFilter.hxx"
#endif


#endif 

