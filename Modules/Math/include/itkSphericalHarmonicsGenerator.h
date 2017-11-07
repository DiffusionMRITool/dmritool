/**
 *       @file  itkSphericalHarmonicsGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  10-29-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */  

#ifndef __itkSphericalHarmonicsGenerator_h
#define __itkSphericalHarmonicsGenerator_h

#include <complex>
#include <itkObject.h>
#include <itkObjectFactory.h>

#include "utlDMRIStoredTables.h"


namespace itk 
{

/** \class SphericalHarmonicsGenerator
 *  \brief Generate complex and real Spherical Harmonic Basis.
 *
 *  The complex Spherical Harmonic basis is defined as
 *  \f[ 
 *    y_l^m = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}} e^{im\phi} P_l^m(\cos\theta) 
 *  \f]
 *  The real SH basis is defined as
 *  \f[ 
 *  Y_l^m(\theta,\phi) = 
 *  \left\{  \begin{array}{lcl}
 *  \sqrt{2}\mbox{Re}(y_l^{|m|}(\theta,\phi))  = \sqrt{2}(-1)^m\mbox{Re}(y_l^m(\theta,\phi)) &  \mbox{if} &  -l\leq m<0  \\  
 *  y_l^m(\theta,\phi)    &  \mbox{if}  &  m=0  \\
 *  \sqrt{2}\mbox{Im}(y_l^m(\theta,\phi))  &  \mbox{if}  &  l\geq m>0 
 *  \end{array}\right. 
 *  \f]
 *  where \f$ y_l^m(\theta,\phi)\f$ is the complex SH basis 
 *
 *  \ingroup Math
 *  \author Jian Cheng
 *
 */
template<class PreciseType = double> 
class ITK_EXPORT SphericalHarmonicsGenerator : public Object
{
public:
  /** Standard class typedefs. */
  typedef SphericalHarmonicsGenerator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SphericalHarmonicsGenerator, Object);

  /** get value of Complex SH basis  */
  static std::complex<PreciseType> ComplexSH(const int l, const int m, const PreciseType theta, const PreciseType phi);
  /** get value of Real SH basis  */
  static PreciseType RealSH(const int l, const int m, const PreciseType theta, const PreciseType phi);

  /** get value of the triple integral of Complex SH basis  */
  static PreciseType ComplexTripleIntegration(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3);
  /** get value of the triple integral of Real SH basis  */
  static PreciseType RealTripleIntegration(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3, const bool is_precalculated=true);

  /** get value of the derivative of theta for Complex SH basis  */
  static std::complex<double> ComplexDerivativeOfTheta ( const int l, const int m, const double theta, const double phi);
  /** get value of the derivative of phi for Complex SH basis  */
  static std::complex<double> ComplexDerivativeOfPhi ( const int l, const int m, const double theta, const double phi);
  
  /** get value of the derivative of theta for Real SH basis  */
  static double RealDerivativeOfTheta ( const int l, const int m, const double theta, const double phi);
  /** get value of the derivative of phi for Real SH basis  */
  static double RealDerivativeOfPhi ( const int l, const int m, const double theta, const double phi);
  
protected:
  SphericalHarmonicsGenerator();
  virtual ~SphericalHarmonicsGenerator();

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private :
  SphericalHarmonicsGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


}

namespace utl 
{
/** Initialize SH3IntegralTable with 3 dimensional sizes from 3 SH rank arguments. 
 * 
 * \param useExactSize If it is false, SH3IntegralTable may have larger size than the input three ranks based on pre-computed table. If it is true, the size is determined by the input SH ranks. 
 * */
inline void
InitializeSHTripleIntegrationTable(const int rank0=-1, const int rank1=-1, const int rank2=-1, const bool useExactSize=false)
{
  int dim[3];
  dim[0]=rank0>0? utl::RankToDimSH(rank0) : 0;
  dim[1]=rank1>0? utl::RankToDimSH(rank1) : 0;
  dim[2]=rank2>0? utl::RankToDimSH(rank2) : 0;

  const unsigned* shapeOld = SH3IntegralTable->GetShape();

  // no useExactSize, SH3IntegralTable is larger than what is requested
  if (!useExactSize && SH3IntegralTable->Size()>0 && dim[0]<=shapeOld[0] && dim[1]<=shapeOld[1] && dim[2]<=shapeOld[2])
    return;

  // SH3IntegralTable is the same as what is requested
  if (SH3IntegralTable->Size()>0 && dim[0]==shapeOld[0] && dim[1]==shapeOld[1] && dim[2]==shapeOld[2])
    return; 

  typedef itk::Image<double,3> ImageType;
  static ImageType::Pointer image = ImageType::New();
  if (itk::IsImageEmpty(image))
    itk::ReadImage<ImageType>(utl::SH3Itegralhdr, image);
  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType pixelIndex;

  unsigned shape[3];
  for ( int i = 0; i < 3; ++i ) 
    shape[i] = useExactSize? dim[i] : utl::max<int>(dim[i], size[i]);
  
  if (SH3IntegralTable->Size()==0 
    || useExactSize 
    || (!useExactSize && SH3IntegralTable->Size()>0 && (dim[0]>SH3IntegralTable->GetShape()[0] || dim[1]>SH3IntegralTable->GetShape()[1] || dim[2]>SH3IntegralTable->GetShape()[2])) )
    {
    SH3IntegralTable->ReSize(shape);

    unsigned index[3];
    for ( int i = 0; i < shape[0]; ++i ) 
      for ( int j = 0; j < shape[1]; ++j ) 
        for ( int k = 0; k < shape[2]; ++k ) 
          {
          index[0]=i, index[1]=j, index[2]=k;
          if (i<size[0] && j<size[1] && k<size[2] )
            {
            pixelIndex[0]=i, pixelIndex[1]=j, pixelIndex[2]=k;
            (*SH3IntegralTable)(index) = image->GetPixel(pixelIndex);
            }
          else
            {
            std::vector<int> v1, v2, v3;
            v1 = utl::GetIndexSHlm(i);
            v2 = utl::GetIndexSHlm(j);
            v3 = utl::GetIndexSHlm(k);
            (*SH3IntegralTable)(index) = itk::SphericalHarmonicsGenerator<double>::RealTripleIntegration(v1[0],v1[1],v2[0],v2[1],v3[0],v3[1],false);
            }
          }
    }
}

}
  
// Define instantiation macro for this template.
#define ITK_TEMPLATE_SphericalHarmonicsGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT SphericalHarmonicsGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef SphericalHarmonicsGenerator< ITK_TEMPLATE_2 x > SphericalHarmonicsGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphericalHarmonicsGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSphericalHarmonicsGenerator_hxx)
# include "itkSphericalHarmonicsGenerator.hxx"
#endif

#endif
