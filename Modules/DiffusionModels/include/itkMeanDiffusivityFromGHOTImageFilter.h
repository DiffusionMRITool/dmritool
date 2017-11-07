/**
 *       @file  itkMeanDiffusivityFromGHOTImageFilter.h
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

#ifndef __itkMeanDiffusivityFromGHOTImageFilter_h
#define __itkMeanDiffusivityFromGHOTImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"
#include "utlCoreMacro.h"

namespace itk
{  
namespace Functor {  
  
/**
 *   \class  MeanDiffusivityFromGHOTCoefficients 
 *   \brief   MD value from GHOT coefficients
 *
 *   \ingroup DiffusionModels
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template< class TInput, class TOutput>
class MeanDiffusivityFromGHOTCoefficients
{
public:
  MeanDiffusivityFromGHOTCoefficients() : m_Scale(-1), m_Tau(ONE_OVER_4_PI_2) {};
  ~MeanDiffusivityFromGHOTCoefficients() {};
  bool operator!=( const MeanDiffusivityFromGHOTCoefficients & other ) const
    {
    return !(*this == other);
    }
  bool operator==( const MeanDiffusivityFromGHOTCoefficients & other ) const
    {
    return other.m_Scale==m_Scale && other.m_Tau==m_Tau;
    }

  inline TOutput operator()( const TInput & A ) const
    {
    TOutput output;
    output = (A[0]/m_Scale) / (4*M_PI*M_PI*m_Tau) * (1.0/ (std::sqrt(4*M_PI)));
    return output;
    }
  void SetScale(const double scale) {this->m_Scale = scale; }
  const double & GetScale() const { return m_Scale; }
  void SetTau(const double tau) {this->m_Tau = tau; }
  const double & GetTau() const { return m_Tau; }
  
  double m_Scale;
  double m_Tau;
};
}

/** \class MeanDiffusivityFromGHOTImageFilter
 *
 * \brief Compute mean diffusivity from GHOT coefficients.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \author Jian Cheng, jian.cheng.1983@gmail.com
 *
 * \ingroup IntensityImageFilters  Multithreaded DiffusionModels
 * \sa UnaryFunctorImageFilter
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT MeanDiffusivityFromGHOTImageFilter :
      public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::MeanDiffusivityFromGHOTCoefficients< 
   typename TInputImage::PixelType,
   typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef MeanDiffusivityFromGHOTImageFilter                 Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Functor::MeanDiffusivityFromGHOTCoefficients< 
      typename TInputImage::PixelType,
      typename TOutputImage::PixelType>   >             Superclass;

  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeanDiffusivityFromGHOTImageFilter, UnaryFunctorImageFilter);

  void SetScale(const double scale)
    {
    if( scale != this->GetFunctor().GetScale() )
      {
      this->GetFunctor().SetScale(scale);
      this->Modified();
      }
    }
  const double & GetScale() const
    {
    return this->GetFunctor().GetScale();
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

protected:
  MeanDiffusivityFromGHOTImageFilter() {};
  virtual ~MeanDiffusivityFromGHOTImageFilter() {};
   
  void PrintSelf(std::ostream &os, Indent indent) const ITK_OVERRIDE
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "Scale = " << static_cast<double>(this->GetScale()) 
      << ", Tau = " << static_cast<double>(this->GetTau())
       << std::endl;
    }

private:
  MeanDiffusivityFromGHOTImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



} // end namespace itk


#endif 
