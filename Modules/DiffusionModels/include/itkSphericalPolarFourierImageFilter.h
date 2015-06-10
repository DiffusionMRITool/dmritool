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
  
  typedef typename Superclass::MatrixType       MatrixType;
  typedef typename Superclass::VectorType       VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType    STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;

  typedef typename Superclass::L2SolverType     L2SolverType;
  typedef typename Superclass::L1FISTASolverType     L1FISTASolverType;
  typedef typename Superclass::EstimationType   EstimationType;
  
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
  
  typedef typename Superclass::MaskImageType      MaskImageType;
  typedef typename Superclass::ScalarImageType    ScalarImageType;
  typedef SphericalPolarFourierRadialGenerator<double> SPFGenerator;
  
  typedef typename Superclass::SamplingSchemeQSpacePointer      SamplingSchemeQSpacePointer;
  typedef typename Superclass::SamplingSchemeQSpaceType         SamplingSchemeQSpaceType;
  
  itkGetMacro(Gn0, STDVectorPointer);
  
  std::vector<int> GetIndexNLM(const int index) const;
  int GetIndexJ(const int n, const int l, const int m) const;
  
  std::vector<int> DimToRank(const int dimm) const;
  int RankToDim(const bool is_radial=false, const int radialRank=-1, const int shRank=-1) const; 
  
  double ComputeScale(const bool setScale=true);
  
  void ComputeRadialMatrix ();
  void ComputeBasisMatrix ();
  void ComputeRegularizationWeight ( );
  
  void ComputeRadialVectorForE0InBasis ( );
  void ComputeRadialVectorForE0InDWI ( );
  
  void SetBasisScale(const double scale);

protected:
  SphericalPolarFourierImageFilter();
  virtual ~SphericalPolarFourierImageFilter() {};
  
  void ComputeBasisMatrixForB0 ();
  
  // void VerifyInputParameters() const;

  void PrintSelf(std::ostream& os, Indent indent) const;
  typename LightObject::Pointer InternalClone() const;
  
  void BeforeThreadedGenerateData ();
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId );
  
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
