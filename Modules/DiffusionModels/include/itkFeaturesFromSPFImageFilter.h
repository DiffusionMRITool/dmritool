/**
 *       @file  itkFeaturesFromSPFImageFilter.h
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


#ifndef __itkFeaturesFromSPFImageFilter_h
#define __itkFeaturesFromSPFImageFilter_h

#include "itkObject.h"
#include "itkMaskedImageToImageFilter.h"
#include "itkSphericalPolarFourierEstimationImageFilter.h"

namespace itk
{  

/**
 *   \class   FeaturesFromSPFImageFilter
 *   \brief Compute some features (DWI/EAP profile, ODFs, scalar indices) from SPF coefficients.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 *  \ingroup DiffusionModels
 *  \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT FeaturesFromSPFImageFilter :
  public MaskedImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FeaturesFromSPFImageFilter         Self;
  typedef MaskedImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( FeaturesFromSPFImageFilter, MaskedImageToImageFilter );
  
  typedef TInputImage                                      InputImageType;
  typedef typename InputImageType::Pointer                 InputImagePointer;
  
  typedef Image<typename TInputImage::InternalPixelType,3> ScalarImageType;
  typedef typename ScalarImageType::Pointer                ScalarImagePointer;
  typedef utl::NDArray<double,2>                           MatrixType;
  typedef utl_shared_ptr<MatrixType>                       MatrixPointer;
  typedef utl::NDArray<double,1>                           VectorType;
  typedef utl_shared_ptr<VectorType>                       VectorPointer;
  typedef std::vector<double>                              STDVectorType;
  typedef utl_shared_ptr<STDVectorType >                   STDVectorPointer;
  
  typedef typename Superclass::MaskImageType               MaskImageType;
    
  typedef SphericalPolarFourierEstimationImageFilter<TInputImage, TOutputImage> SPFIFilterBaseType;
  
  typedef enum 
    {
    SPF=0, 
    DSPF
    } BasisType;
  
  itkSetMacro(SHRank, int);
  itkGetMacro(SHRank, int);
  itkSetMacro(RadialRank, int);
  itkGetMacro(RadialRank, int);
  
  itkSetMacro(Tau, double);
  itkGetMacro(Tau, double);
  itkSetMacro(MD0, double);
  itkGetMacro(MD0, double);

  itkSetMacro(BasisType, BasisType);
  itkGetMacro(BasisType, BasisType);
  
  itkSetMacro(Orientations, MatrixPointer);
  itkGetMacro(Orientations, MatrixPointer);

  double ComputeScale(const bool setScale=true);
  
  virtual void SetBasisScale(const double scale);
  itkGetMacro(BasisScale, double);
  
  /** Set/Get the scale image, which is normally determined by MDImage.  */
  void SetScaleImage(const ScalarImagePointer& scaleImage);
  // itkSetObjectMacro(ScaleImage, ScalarImageType);
  itkGetObjectMacro(ScaleImage, ScalarImageType);
  itkGetConstObjectMacro(ScaleImage, ScalarImageType);
  
  itkSetMacro(IsInQSpace,bool);
  itkGetMacro(IsInQSpace,bool);
  itkBooleanMacro(IsInQSpace);

  /** for dwi/eap profile  */
  itkSetMacro(IsFourier,bool);
  itkGetMacro(IsFourier,bool);
  itkBooleanMacro(IsFourier);
  
  
  itkGetMacro(SPFToFeatureTransform, MatrixPointer);
  
  virtual void ComputeSPFToFeatureTransform() {}
  

protected:
  FeaturesFromSPFImageFilter();
  virtual ~FeaturesFromSPFImageFilter() {};

  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  
  virtual void VerifyInputParameters() const ITK_OVERRIDE;

  void BeforeThreadedGenerateData () ITK_OVERRIDE;
  // void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId );
  
  void SetSPFIEstimator();

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  
  double m_BasisScale;
  double m_MD0;
  double m_Tau;
  int m_SHRank;
  int m_RadialRank;
  typename ScalarImageType::Pointer m_ScaleImage;

  /** the basis is in qspace or in the fourier space  */
  bool m_IsInQSpace;
  
  /** for profile  */
  bool m_IsFourier;
  
  MatrixPointer m_SPFToFeatureTransform;
  BasisType m_BasisType;
  
  /** Orientaitons in spherical format. If it is set, the output will be samples, not SH coefficients  */
  MatrixPointer m_Orientations;
    
  typename SPFIFilterBaseType::Pointer m_SPFIEstimator;

private:
  FeaturesFromSPFImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkFeaturesFromSPFImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkFeaturesFromSPFImageFilter_hxx)
#include "itkFeaturesFromSPFImageFilter.hxx"
#endif




#endif 

