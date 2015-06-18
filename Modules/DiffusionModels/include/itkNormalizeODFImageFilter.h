/**
 *       @file  itkNormalizeODFImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkNormalizeODFImageFilter_h
#define __itkNormalizeODFImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
namespace Functor
{
/** \ingroup DiffusionModels  */
template <class TInput, class TOutput>
class ODFNormlizeFunctor
{
public:
  ODFNormlizeFunctor()
  {
  m_ODFType=SH;
  }

  ~ODFNormlizeFunctor()
  {
  }

  bool operator!=(const ODFNormlizeFunctor & other) const
  {
    return !( *this == other );
  }

  bool operator==(const ODFNormlizeFunctor & other) const
  {
    return other.m_ODFType == m_ODFType;
  }
  
  typedef enum 
    {
    SH=0,
    SAMPLE 
    } ODFType;

  inline TOutput operator()(const TInput & A) const
  {
    unsigned int size = A.GetSize();
    TOutput value;
    value.SetSize( size );

    double norm = 0;
    for ( unsigned int k=0; k < size; k++ )
      norm += std::fabs(A[k]);      
    if (norm>0)
      {
      double normFactor = 1.0;
      if (m_ODFType==SH && std::fabs(A[0])>1e-8)
        {
        normFactor = 1.0/std::sqrt(4*M_PI) * 1.0/A[0];
        }
      if (m_ODFType==SAMPLE)
        {
        double sum = 0;
        for ( unsigned int k=0; k < size; k++ )
          sum += A[k];
        if (std::fabs(sum)>1e-8)
          normFactor = 1.0/sum * double(size)/(4*M_PI);
        }
      for ( unsigned int k=0; k < size; k++ )
        value[k] = A[k]*normFactor;
      }
    else
      {
      for ( unsigned int k=0; k < size; k++ )
        value[k] = 0;
      }
    return value;
  }

  void SetODFType(ODFType odfType)
  {
    this->m_ODFType = odfType;
  }

  const ODFType & GetODFType() const
  {
    return m_ODFType;
  }

  ODFType m_ODFType;

};
}

/** \class NormalizeODFImageFilter
 *
 * \brief Multiply input pixels by a constant.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \ingroup IntensityImageFilters  Multithreaded DiffusionModels
 * \sa UnaryFunctorImageFilter
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT NormalizeODFImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::ODFNormlizeFunctor< typename TInputImage::PixelType,typename TOutputImage::PixelType>  >
{
public:
  /** Standard class typedefs. */
  typedef NormalizeODFImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Functor::ODFNormlizeFunctor<
        typename TInputImage::PixelType, typename TOutputImage::PixelType> >             Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( NormalizeODFImageFilter, UnaryFunctorImageFilter );
  
  typedef Functor::ODFNormlizeFunctor< typename TInputImage::PixelType,
    typename TOutputImage::PixelType>   FunctorType;

  typedef typename FunctorType::ODFType ODFType;

  /** Set the constant that will be used to multiply all the image pixels */
  void SetODFType( ODFType type )
  {
    if( type != this->GetFunctor().GetODFType() )
      {
      this->GetFunctor().SetODFType( type );
      this->Modified();
      }
  }

  const ODFType & GetODFType() const
  {
    return this->GetFunctor().GetODFType();
  }

protected:
  NormalizeODFImageFilter()
  {
  }

  virtual ~NormalizeODFImageFilter()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent ) const
  {
    Superclass::PrintSelf( os, indent );

    os << indent << "ODFType: " << this->GetODFType() << std::endl;
  }

private:
  NormalizeODFImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};
} // end namespace itk


#endif 

