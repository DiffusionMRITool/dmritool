/**
 *       @file  itkScalarMapFromSPFImageFilter.h
 *      @brief  
 *     Created  "03-18-2014
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkScalarMapFromSPFImageFilter_h
#define __itkScalarMapFromSPFImageFilter_h

#include "itkFeaturesFromSPFImageFilter.h"


namespace itk
{  

/**
 *   \class   ScalarMapFromSPFImageFilter
 *   \brief   calculate ODFs from SPF coefficients
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage=Image<double,3> >
class ITK_EXPORT ScalarMapFromSPFImageFilter :
public FeaturesFromSPFImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ScalarMapFromSPFImageFilter         Self;
  typedef FeaturesFromSPFImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( ScalarMapFromSPFImageFilter, FeaturesFromSPFImageFilter );
  
  typedef typename Superclass::ScalarImageType             ScalarImageType;
  typedef typename Superclass::ScalarImagePointer          ScalarImagePointer;
  typedef typename Superclass::MatrixType                  MatrixType;
  typedef typename Superclass::MatrixPointer               MatrixPointer;
  typedef typename Superclass::VectorType                  VectorType;
  typedef typename Superclass::BasisType                   BasisType;
  typedef typename Superclass::STDVectorType               STDVectorType;
  typedef typename Superclass::STDVectorPointer            STDVectorPointer;
  
  typedef typename Superclass::MaskImageType               MaskImageType;

  typedef enum 
    {
    /** return to origin probability  */
    RTO=0,  
    /** Mean Squared Displacement  */
    MSD,
    /** generalized propagator FA for 3D EAP  */
    PFA
    } MapType;
  
  itkSetMacro(MapType, MapType);
  itkGetMacro(MapType, MapType);

protected:
  ScalarMapFromSPFImageFilter() : Superclass()
  {
  }

  virtual ~ScalarMapFromSPFImageFilter() {};
  
  void VerifyInputParameters() const ITK_OVERRIDE{}
  
  void GenerateOutputInformation() ITK_OVERRIDE;
  
  void BeforeThreadedGenerateData () ITK_OVERRIDE;
  void ThreadedGenerateData(const typename TOutputImage::RegionType& outputRegionForThread,ThreadIdType threadId ) ITK_OVERRIDE;
  
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  VectorType m_SumWeight;
  MapType m_MapType;

private:
  ScalarMapFromSPFImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented
  
};



} // end namespace itk


#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkScalarMapFromSPFImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkScalarMapFromSPFImageFilter_hxx)
#include "itkScalarMapFromSPFImageFilter.hxx"
#endif


#endif 
