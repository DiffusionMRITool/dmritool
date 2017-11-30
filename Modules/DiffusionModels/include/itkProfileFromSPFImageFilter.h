/**
 *       @file  itkProfileFromSPFImageFilter.h
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

#ifndef __itkProfileFromSPFImageFilter_h
#define __itkProfileFromSPFImageFilter_h

#include "itkFeaturesFromSPFImageFilter.h"


namespace itk
{  

/**
 *   \class   ProfileFromSPFImageFilter
 *   \brief   Calculate DWI/EAP profile from SPF coefficients
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT ProfileFromSPFImageFilter :
public FeaturesFromSPFImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ProfileFromSPFImageFilter         Self;
  typedef FeaturesFromSPFImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( ProfileFromSPFImageFilter, FeaturesFromSPFImageFilter );
  
  itkTypedefMaskedImageToImageMacro(Superclass);

  typedef typename Superclass::MatrixType                  MatrixType;
  typedef typename Superclass::MatrixPointer               MatrixPointer;
  typedef typename Superclass::VectorType                  VectorType;
  typedef typename Superclass::VectorPointer               VectorPointer;
  typedef typename Superclass::BasisType                   BasisType;
  typedef typename Superclass::STDVectorType               STDVectorType;
  typedef typename Superclass::STDVectorPointer            STDVectorPointer;

  
  itkSetMacro(Radius, double);
  itkGetMacro(Radius, double);
  
  itkSetMacro(RadiusVector, STDVectorPointer);
  itkGetMacro(RadiusVector, STDVectorPointer);
  
  void ComputeSPFToFeatureTransform() ITK_OVERRIDE;

protected:
  ProfileFromSPFImageFilter() : Superclass(), 
  m_RadiusVector(new STDVectorType())
    {
    m_Radius=-1;
    }
  virtual ~ProfileFromSPFImageFilter() {};
  
  void VerifyInputParameters() const ITK_OVERRIDE;
  void GenerateOutputInformation() ITK_OVERRIDE;

  void BeforeThreadedGenerateData () ITK_OVERRIDE;
  void ThreadedGenerateData(const typename TOutputImage::RegionType& outputRegionForThread,ThreadIdType threadId ) ITK_OVERRIDE;

  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

private:
  ProfileFromSPFImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

  /** b value or r value, which is dependent on m_IsInQSpace  */
  double m_Radius;

  /** If it is set, the output will be samples, not SH coefficients  */
  STDVectorPointer m_RadiusVector;
};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkProfileFromSPFImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkProfileFromSPFImageFilter_hxx)
#include "itkProfileFromSPFImageFilter.hxx"
#endif




#endif 

