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

  typedef typename Superclass::MatrixType       MatrixType;
  typedef typename Superclass::VectorType       VectorType;
  typedef typename Superclass::MatrixPointer    MatrixPointer;
  typedef typename Superclass::VectorPointer    VectorPointer;
  typedef typename Superclass::STDVectorType    STDVectorType;
  typedef typename Superclass::STDVectorPointer STDVectorPointer;

  typedef typename Superclass::L2SolverType     L2SolverType;
  typedef typename Superclass::EstimationType   EstimationType;

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::IndexType      InputImageIndexType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SpacingType    InputImageSpacingType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  
  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename OutputImageType::IndexType     OutputImageIndexType;
  typedef typename OutputImageType::SizeType      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType   OutputImageSpacingType;
  typedef typename OutputImageType::PixelType     OutputImagePixelType;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;

  typedef typename Superclass::MaskImageType      MaskImageType;
  typedef typename Superclass::ScalarImageType    ScalarImageType;
  
  void VerifyInputParameters() const;
  
  std::vector<int> DimToRank(const int dimm) const;
  int RankToDim(const bool is_radial=false, const int radialRank=-1, const int shRank=-1) const; 
  
  double ComputeScale(const bool setScale=true);
  
  void ComputeRadialMatrix ();
  void ComputeBasisMatrix ();
  void ComputeRegularizationWeight ( );

  

protected:
  GeneralizedHighOrderTensorImageFilter();
  virtual ~GeneralizedHighOrderTensorImageFilter() {};

  void PrintSelf(std::ostream& os, Indent indent) const;

  void BeforeThreadedGenerateData ();

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId );


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
