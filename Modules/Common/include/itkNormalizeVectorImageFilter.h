/**
 *       @file  itkNormalizeVectorImageFilter.h
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

#ifndef __itkNormalizeVectorImageFilter_h
#define __itkNormalizeVectorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
  
/** \class NormalizeVectorImageFilter
 * \brief Pixel-wise vector normalization
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {  
  
template< class TInput, class TOutput>
class NormalizeVectorImageFunctor
{
public:
  //typedef typename TInput::RealValueType RealType;
  
  NormalizeVectorImageFunctor() {m_NormalizeType=NONE;}
  ~NormalizeVectorImageFunctor() {}
  bool operator!=( const NormalizeVectorImageFunctor & other) const
    {
    return !( *this == other );
    }
  bool operator==( const NormalizeVectorImageFunctor & other ) const
    {
    return m_NormalizeType == other.m_NormalizeType;
    }
  typedef enum 
    {
    NONE=0,
    SUM, 
    L1NORM,
    L2NORM
    } NormalizeType;

  inline TOutput operator()( const TInput & A ) const
    {
    unsigned int vectorSize = A.GetSize();
    
    TOutput output;
    output.SetSize( vectorSize );

    double norm = 0;
        
    for ( unsigned int k=0; k < vectorSize; k++ )
      norm += std::fabs(A[k]);      
    
    if ( norm > 0 )
      {
      if (m_NormalizeType==NONE)
        norm = 1.0;
      else if (m_NormalizeType==L1NORM)
        {
        }
      else if (m_NormalizeType==SUM)
        {
        norm=0;
        for ( unsigned int k=0; k < vectorSize; k++ )
          norm += A[k];
        }
      else if (m_NormalizeType==L2NORM)
        {
        norm=0;
        for ( unsigned int k=0; k < vectorSize; k++ )
          norm += A[k]*A[k];
        norm = std::sqrt(norm);
        }
      if (norm!=0)
        {
        for ( unsigned int k=0; k < vectorSize; k++ )
          output[k] = A[k] / norm;      
        }
      else
        for ( unsigned int k=0; k < vectorSize; k++ )
          output[k] = 0;      
      }
    else
      {
      for ( unsigned int k=0; k < vectorSize; k++ )
        output[k] = 0;      
      }
    
    return output;
    }

  void SetNormalizeType( NormalizeType norm)
  {
    m_NormalizeType = norm;
  }
  NormalizeType GetNormalizeType() const
  {
    return m_NormalizeType;
  }
  
private:

  NormalizeType m_NormalizeType;
}; 
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT NormalizeVectorImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
    Functor::NormalizeVectorImageFunctor< typename TInputImage::PixelType, 
    typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef NormalizeVectorImageFilter                     Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Functor::NormalizeVectorImageFunctor< typename TInputImage::PixelType, 
                    typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(NormalizeVectorImageFilter, 
               UnaryFunctorImageFilter);

  typedef Functor::NormalizeVectorImageFunctor< typename TInputImage::PixelType,
    typename TOutputImage::PixelType>   FunctorType;

  typedef typename FunctorType::NormalizeType NormalizeType;

  void SetNormalizeType(NormalizeType norm)
    {
    this->GetFunctor().SetNormalizeType(norm);
    this->Modified();
    }

  NormalizeType GetNormalizeType() const
    {
    return this->GetFunctor().GetNormalizeType();
    }

protected:
  NormalizeVectorImageFilter() 
    {
    }

  virtual ~NormalizeVectorImageFilter() {}


private:
  NormalizeVectorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk



#endif 


