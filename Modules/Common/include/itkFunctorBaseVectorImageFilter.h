/**
 *       @file  itkFunctorBaseVectorImageFilter.h
 *      @brief  
 *     Created  "09-13-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef itkFunctorBaseVectorImageFilter_h
#define itkFunctorBaseVectorImageFilter_h

#include "itkMaskedImageToImageFilter.h"
#include "utlCoreMacro.h"
#include "utlTypeinfo.h"

namespace itk
{
/** \class FunctorBaseVectorImageFilter
 * \brief Implements vector-valued generic operation on one image.
 *
 *  Calculation loops over image domain and extracts vectors along x/y/z/t-axis. 
 *
 *  Mask can be a 3D image, or 4D image.
 *
 *  m_Functor is a functor with utl::Vector as input and utl::Vector as output. 
 *  m_Functor needs to define GetOutputDimension to obtain the NumberOfComponentsPerPixel in output image. 
 *
 * \author Jian Cheng
 *
 * \ingroup ITKCommon
 *
 */
template< typename TInputImage, typename TOutputImage, typename TFunction, class TMaskImage=Image<double,3> >
class ITK_EXPORT FunctorBaseVectorImageFilter
  :public MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
{
public:
  /** Standard class typedefs. */
  typedef FunctorBaseVectorImageFilter                         Self;
  typedef MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FunctorBaseVectorImageFilter, MaskedImageToImageFilter);

  /** Some typedefs. */
  typedef TFunction FunctorType;

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::IndexType      InputImageIndexType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SpacingType    InputImageSpacingType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  
  typedef TOutputImage                            OutputImageType;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename OutputImageType::IndexType     OutputImageIndexType;
  typedef typename OutputImageType::SizeType      OutputImageSizeType;
  typedef typename OutputImageType::SpacingType   OutputImageSpacingType;
  typedef typename OutputImageType::PixelType     OutputImagePixelType;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  
  typedef typename Superclass::MaskImageType      MaskImageType;

  itkSetGetMacro(VectorAxis, int);

  /** Get the functor object.  The functor is returned by reference.
   * (Functors do not have to derive from itk::LightObject, so they do
   * not necessarily have a reference count. So we cannot return a
   * SmartPointer.) */
  FunctorType &       GetFunctor() { return m_Functor; }
  const FunctorType & GetFunctor() const { return m_Functor; }

  /** Set the functor object.  This replaces the current Functor with a
   * copy of the specified Functor. This allows the user to specify a
   * functor that has ivars set differently than the default functor.
   * This method requires an operator!=() be defined on the functor
   * (or the compiler's default implementation of operator!=() being
   * appropriate). */
  void SetFunctor(const FunctorType & functor)
  {
    if ( m_Functor != functor )
      {
      m_Functor = functor;
      this->Modified();
      }
  }

protected:
  FunctorBaseVectorImageFilter() 
    {
    this->SetNumberOfRequiredInputs(1);
    }
  virtual ~FunctorBaseVectorImageFilter() {}


protected: 
  
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE
    {
    Superclass::PrintSelf(os, indent);
    utlLogOSVar(os << indent, m_VectorAxis);
    m_Functor.Print(os <<indent << "m_Functor : "<<  utl::TypeName(m_Functor)<< " ("<< &m_Functor <<") = \n");
    }
  
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();

    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_Functor = m_Functor;
    rval->m_VectorAxis = m_VectorAxis;
    return loPtr;
    }

  /** Functor work on utl::Vector  */
  FunctorType m_Functor;

  /** vector axis. 0/1/2/3 means x/y/z/t-axis */
  int m_VectorAxis=3;

private:
  FunctorBaseVectorImageFilter(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};
} // end namespace itk

#endif
