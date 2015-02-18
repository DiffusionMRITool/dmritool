/**
 *       @file  itkSHCoefficientsToGFAImaeFilter.h
 *      @brief  In each vxoel, calculate gfa from SH coefficients
 *     Created  "02-12-2014
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkSHCoefficientsToGFAImaeFilter_h
#define __itkSHCoefficientsToGFAImaeFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include <cmath>

namespace itk
{
/** \class SHCoefficientsToGFAImaeFilter
 *
 * \brief In each vxoel, calculate gfa from SH coefficients.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 */
namespace Functor
{
template <class TInput, class TOutput=double>
class SHCoefficientsToGFA
{
public:

  SHCoefficientsToGFA()
  {
  }

  ~SHCoefficientsToGFA()
  {
  }

  bool operator!=(const SHCoefficientsToGFA & other) const
  {
    return !( *this == other );
  }

  bool operator==(const SHCoefficientsToGFA & other) const
  {
    return false;
  }

  inline TOutput operator()(const TInput & A)
  {
  double norm2 = A.GetSquaredNorm();
  if (norm2>0)
    return std::sqrt(1 - A[0]*A[0]/norm2 );
  else
    return 0;
  }


};
}

template <class TInputImage, class TOutputImage=itk::Image<double,3> >
class ITK_EXPORT SHCoefficientsToGFAImaeFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::SHCoefficientsToGFA<
                            typename TInputImage::PixelType,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef SHCoefficientsToGFAImaeFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Functor::SHCoefficientsToGFA<
        typename TInputImage::PixelType,  typename TOutputImage::PixelType> 
      >             Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SHCoefficientsToGFAImaeFilter, UnaryFunctorImageFilter );

protected:
  SHCoefficientsToGFAImaeFilter()
  {
  }

  virtual ~SHCoefficientsToGFAImaeFilter()
  {
  }

private:
  SHCoefficientsToGFAImaeFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};
} // end namespace itk


#endif 
