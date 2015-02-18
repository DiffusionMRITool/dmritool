/**
 *       @file  itkSPFScaleFromMeanDiffusivityImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-04-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSPFScaleFromMeanDiffusivityImageFilter_h
#define __itkSPFScaleFromMeanDiffusivityImageFilter_h



#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{  
/** \class SPFScaleFromMeanDiffusivityImageFilter
 *
 * \brief Compute SPF scale from mean diffusivity.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \author Jian Cheng, jian.cheng.1983@gmail.com
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 */
namespace Functor {  
  
template< class TInput, class TOutput>
class SPFScaleFromMeanDiffusivity
{
public:
  SPFScaleFromMeanDiffusivity() : m_MD0(0.7e-3), m_Tau(ONE_OVER_4_PI_2), m_IsOriginalBasis(true) {};
  ~SPFScaleFromMeanDiffusivity() {};
  bool operator!=( const SPFScaleFromMeanDiffusivity & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const SPFScaleFromMeanDiffusivity & other ) const
    {
    return other.m_MD0==m_MD0 || other.m_Tau==m_Tau;
    }

  inline TOutput operator()( const TInput & A ) const
    {
    if (A<1e-10)
      return 0.0;
      /** NOTE: in real data, the estimated ADC m_MDImage(x,y,z) may be negative or very small value due to inappropriate orders in GHOT model. 
       * This will result in problem when calculating 1F1 for EAP profile. So set these irregular values as the typical D. 
       * Note that this is just a work around. Need to check the range of ADC map after estimating the ADC. 
       * */
    TOutput scale;
    if (m_IsOriginalBasis)
      scale = 1.0 / (8*M_PI*M_PI*this->m_Tau* (A>this->m_MD0*0.01? A : this->m_MD0) ) ;  // 714.29 (700) scale for SPF basis, dual scale is 1/(4*pi^2*scale)
    else
      scale = 2*this->m_Tau*(A>this->m_MD0*0.01? A : this->m_MD0);  // 3.5462e-5 scale for SPF basis, dual scale is 1/(4*pi^2*scale)
    return scale;
    }
  void SetMD0(const double md0) {this->m_MD0 = md0; }
  const double & GetMD0() const { return m_MD0; }
  void SetTau(const double tau) {this->m_Tau = tau; }
  const double & GetTau() const { return m_Tau; }
  void SetIsOriginalBasis(const bool isSPF) {this->m_IsOriginalBasis = isSPF; }
  const bool & GetIsOriginalBasis() const { return m_IsOriginalBasis; }
  
  double m_MD0;
  double m_Tau;
  bool m_IsOriginalBasis;
};
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT SPFScaleFromMeanDiffusivityImageFilter :
      public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::SPFScaleFromMeanDiffusivity< 
   typename TInputImage::PixelType,
   typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef SPFScaleFromMeanDiffusivityImageFilter                 Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Functor::SPFScaleFromMeanDiffusivity< 
      typename TInputImage::PixelType,
      typename TOutputImage::PixelType>   >             Superclass;

  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SPFScaleFromMeanDiffusivityImageFilter, UnaryFunctorImageFilter);

  void SetMD0(const double md0)
    {
    if( md0 != this->GetFunctor().GetMD0() )
      {
      this->GetFunctor().SetMD0(md0);
      this->Modified();
      }
    }
  const double & GetMD0() const
    {
    return this->GetFunctor().GetMD0();
    }
  
  void SetTau(const double tau)
    {
    if( tau != this->GetFunctor().GetTau() )
      {
      this->GetFunctor().SetTau(tau);
      this->Modified();
      }
    }
  const double & GetTau() const
    {
    return this->GetFunctor().GetTau();
    }
  
  void SetIsOriginalBasis(const bool isSPF)
    {
    if( isSPF != this->GetFunctor().GetIsOriginalBasis() )
      {
      this->GetFunctor().SetIsOriginalBasis(isSPF);
      this->Modified();
      }
    }
  const bool & GetIsOriginalBasis() const
    {
    return this->GetFunctor().GetIsOriginalBasis();
    }

protected:
  SPFScaleFromMeanDiffusivityImageFilter() {};
  virtual ~SPFScaleFromMeanDiffusivityImageFilter() {};
   
  void PrintSelf(std::ostream &os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "md0 = " << static_cast<double>(this->GetMD0()) 
      << ", Tau = " << static_cast<double>(this->GetTau())
       << std::endl;
    }

private:
  SPFScaleFromMeanDiffusivityImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



} // end namespace itk



#endif 
