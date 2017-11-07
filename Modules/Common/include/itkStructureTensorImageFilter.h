/**
 *       @file  itkStructureTensorImageFilter.h
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#ifndef itkStructureTensorImageFilter_h
#define itkStructureTensorImageFilter_h

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkGradientImageFilter.h"

namespace itk
{
/**
 * \class StructureTensorImageFilter
 *
 * \brief Computes the structure tensor.
 *
 * Implementation of the structure tensor, defined by
 *
 * \f[K_\rho (\nabla u_\sigma \otimes \nabla u_\sigma),\f]
 *
 * where \f$K_\rho\f$ denotes the gaussian kernel of standard deviation \f$\rho\f$,
 * and \f$u_\sigma := K_\sigma * u\f$.
 *
 * \ingroup ImageFilter
 */
template< class TInputImage, class TOutputImage >
class StructureTensorImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  typedef StructureTensorImageFilter          Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /// Method for creation through the object factory.
  itkNewMacro(Self);
  /// Run-time type information (and related methods).
  itkTypeMacro(StructureTensorImageFilter, Superclass);

  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  typedef itk::GradientImageFilter<InputImageType, double, double, OutputImageType> GradientFilterType;

  itkSetMacro(IntensityScale, double);
  itkGetMacro(IntensityScale, double);

protected:
  virtual void GenerateData() ITK_OVERRIDE;


  StructureTensorImageFilter();
  
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "m_IntensityScale = " << m_IntensityScale << std::endl;
    }

  /** scale the gradient.  */
  double m_IntensityScale;

private:
  StructureTensorImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                      // purposely not implemented
};

} // end namespace itk

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkStructureTensorImageFilter_hxx)
#include "itkStructureTensorImageFilter.hxx"
#endif

#endif
