/**
 *       @file  itkGeneralizedHighOrderTensorImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-22-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkGeneralizedHighOrderTensorImageFilter_h
#define __itkGeneralizedHighOrderTensorImageFilter_h

#include "itkSphericalPolarFourierEstimationImageFilter.h"


namespace itk
{

/**
 *   \class   GeneralizedHighOrderTensorImageFilter
 *   \brief   estimate the coefficients in GHOT model, which is used for md image estimation
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TInputImage, class TOutputImage=TInputImage >
class ITK_EXPORT GeneralizedHighOrderTensorImageFilter :
  public SphericalPolarFourierEstimationImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GeneralizedHighOrderTensorImageFilter         Self;
  typedef SphericalPolarFourierEstimationImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( GeneralizedHighOrderTensorImageFilter, SphericalPolarFourierEstimationImageFilter );

  itkTypedefMaskedImageToImageMacro(Superclass);

  typedef typename Superclass::MatrixType       MatrixType;
  typedef typename Superclass::VectorType       VectorType;
  typedef typename Superclass::MatrixPointer    MatrixPointer;
  typedef typename Superclass::VectorPointer    VectorPointer;
  typedef typename Superclass::STDVectorType    STDVectorType;
  typedef typename Superclass::STDVectorPointer STDVectorPointer;

  typedef typename Superclass::L2SolverType     L2SolverType;
  typedef typename Superclass::EstimationType   EstimationType;
  
  void VerifyInputParameters() const ITK_OVERRIDE;
  
  std::vector<int> DimToRank(const int dimm) const ITK_OVERRIDE;
  int RankToDim(const bool is_radial=false, const int radialRank=-1, const int shRank=-1) const ITK_OVERRIDE; 
  
  double ComputeScale(const bool setScale=true) ITK_OVERRIDE;
  
  void ComputeRadialMatrix () ITK_OVERRIDE;
  void ComputeBasisMatrix () ITK_OVERRIDE;
  void ComputeRegularizationWeight ( ) ITK_OVERRIDE;

  

protected:
  GeneralizedHighOrderTensorImageFilter();
  virtual ~GeneralizedHighOrderTensorImageFilter() {};

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  void BeforeThreadedGenerateData () ITK_OVERRIDE;

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId ) ITK_OVERRIDE;


private:
  GeneralizedHighOrderTensorImageFilter(const Self&);//purposely not implemented
  void operator=(const Self&);//purposely not implemented

};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkGeneralizedHighOrderTensorImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkGeneralizedHighOrderTensorImageFilter_hxx)
#include "itkGeneralizedHighOrderTensorImageFilter.hxx"
#endif

#endif 
