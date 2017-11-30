/**
 *       @file  itkSphericalPolarFourierImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-26-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSphericalPolarFourierImageFilter_h
#define __itkSphericalPolarFourierImageFilter_h

#include "itkSphericalPolarFourierEstimationImageFilter.h"
#include "itkSphericalPolarFourierGenerator.h"


namespace itk
{

/**
 *   \class   SphericalPolarFourierImageFilter
 *   \brief   estimate the coefficients in SPF model
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT SphericalPolarFourierImageFilter :
public SphericalPolarFourierEstimationImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SphericalPolarFourierImageFilter         Self;
  typedef SphericalPolarFourierEstimationImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( SphericalPolarFourierImageFilter, SphericalPolarFourierEstimationImageFilter );
  
  itkTypedefMaskedImageToImageMacro(Superclass);
  
  typedef typename Superclass::MatrixType       MatrixType;
  typedef typename Superclass::VectorType       VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType    STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;

  typedef typename Superclass::L2SolverType     L2SolverType;
  typedef typename Superclass::L1FISTASolverType     L1FISTASolverType;
  typedef typename Superclass::EstimationType   EstimationType;
  
  typedef SphericalPolarFourierRadialGenerator<double> SPFGenerator;
  
  typedef typename Superclass::SamplingSchemeQSpacePointer      SamplingSchemeQSpacePointer;
  typedef typename Superclass::SamplingSchemeQSpaceType         SamplingSchemeQSpaceType;
  
  itkGetMacro(Gn0, STDVectorPointer);
  
  std::vector<int> GetIndexNLM(const int index) const ITK_OVERRIDE;
  int GetIndexJ(const int n, const int l, const int m) const ITK_OVERRIDE;
  
  std::vector<int> DimToRank(const int dimm) const ITK_OVERRIDE;
  int RankToDim(const bool is_radial=false, const int radialRank=-1, const int shRank=-1) const ITK_OVERRIDE; 
  
  double ComputeScale(const bool setScale=true) ITK_OVERRIDE;
  
  void ComputeRadialMatrix () ITK_OVERRIDE;
  void ComputeBasisMatrix () ITK_OVERRIDE;
  void ComputeRegularizationWeight ( ) ITK_OVERRIDE;
  
  void ComputeRadialVectorForE0InBasis ( );
  void ComputeRadialVectorForE0InDWI ( );
  
  void SetBasisScale(const double scale) ITK_OVERRIDE;

protected:
  SphericalPolarFourierImageFilter();
  virtual ~SphericalPolarFourierImageFilter() {};
  
  void ComputeBasisMatrixForB0 ();
  
  // void VerifyInputParameters() const;

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  
  void BeforeThreadedGenerateData () ITK_OVERRIDE;
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId ) ITK_OVERRIDE;
  
  STDVectorPointer m_Gn0;
  VectorPointer m_G0DWI;

private:
  SphericalPolarFourierImageFilter(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

};


}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphericalPolarFourierImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSphericalPolarFourierImageFilter_hxx)
#include "itkSphericalPolarFourierImageFilter.hxx"
#endif


#endif 
