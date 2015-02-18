/**
 *       @file  itkDiffusionModelGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-02-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkDiffusionModelGenerator_h
#define __itkDiffusionModelGenerator_h

#include "itkObject.h"
#include "itkSamplingSchemeQSpace.h"

namespace itk
{


template<class PreciseType = double> 
class ITK_EXPORT DiffusionModelGenerator : public Object
{
public:
  /** Standard class typedefs. */
  typedef DiffusionModelGenerator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(DiffusionModelGenerator, Object);
  
  typedef vnl_matrix<double>                      MatrixType;
  typedef vnl_vector<double>                      VectorType;
  typedef utl_shared_ptr<MatrixType>              MatrixPointer;
  typedef utl_shared_ptr<VectorType>              VectorPointer;
  typedef std::vector<double>                     STDVectorType;
  typedef utl_shared_ptr<STDVectorType >          STDVectorPointer;
  
  typedef SamplingSchemeQSpace<double>                       SamplingSchemeQSpaceType;
  typedef typename SamplingSchemeQSpaceType::Pointer         SamplingSchemeQSpacePointer;
  
  typedef SamplingScheme3D<double>                           SamplingSchemeRSpaceType;
  typedef typename SamplingSchemeRSpaceType::Pointer         SamplingSchemeRSpacePointer;
  
  typedef SamplingSchemeQSpaceType::PointType                PointType;
  
  itkSetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkGetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkSetObjectMacro(SamplingSchemeRSpace, SamplingSchemeRSpaceType);
  itkGetObjectMacro(SamplingSchemeRSpace, SamplingSchemeRSpaceType);
  
  itkSetMacro(ODFOrder, int);
  itkGetMacro(ODFOrder, int);
  
  itkGetMacro(DWISamples, VectorPointer);
  itkGetMacro(ODFSamples, VectorPointer);
  itkGetMacro(EAPSamples, VectorPointer);
  
  virtual void Rotate (const MatrixType& mat)
    {
    }

  virtual void ComputeDWISamples ()
    {
    }

  virtual void ComputeEAPSamples ()
    {
    }
  
  virtual void ComputeODFSamples ()
    {
    }

  virtual void VerifyInputParameters() const
    {
    utlGlobalException(m_SamplingSchemeQSpace->GetNumberOfSamples()>0 && m_SamplingSchemeRSpace->GetNumberOfSamples()>0 
      && (std::fabs(m_SamplingSchemeQSpace->GetDeltaBig()- m_SamplingSchemeRSpace->GetDeltaBig())>1e-8 || std::fabs(m_SamplingSchemeQSpace->GetDeltaSmall()- m_SamplingSchemeRSpace->GetDeltaSmall())>1e-8)
      , "inconsistent m_DeltaBig or m_DeltaSmall in m_SamplingSchemeQSpace and m_SamplingSchemeRSpace");
    }

protected:
  /** DiffusionModelGenerator constructor  */
  DiffusionModelGenerator () : Superclass(), 
    m_DWISamples(new VectorType()),
    m_ODFSamples(new VectorType()),
    m_EAPSamples(new VectorType())
    {
    m_ODFOrder = 2;
    m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
    m_SamplingSchemeRSpace = SamplingSchemeRSpaceType::New();
    }

  virtual ~DiffusionModelGenerator ()
    {
    }

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    if (m_SamplingSchemeQSpace->GetNumberOfSamples()>0)
      {
      os << indent << "m_SamplingSchemeQSpace = " << m_SamplingSchemeQSpace << std::endl;
      utl::PrintVnlVector(*m_DWISamples, "m_DWISamples", " ", os<<indent);
      }
    if (m_SamplingSchemeRSpace->GetNumberOfSamples()>0)
      {
      os << indent << "m_SamplingSchemeRSpace = " << m_SamplingSchemeRSpace << std::endl;
      PrintVar1(true, m_ODFOrder, os<<indent);
      utl::PrintVnlVector(*m_EAPSamples, "m_EAPSamples", " ", os<<indent);
      utl::PrintVnlVector(*m_ODFSamples, "m_ODFSamples", " ", os<<indent);
      }
    }

  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();

    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }

    rval->m_SamplingSchemeQSpace = m_SamplingSchemeQSpace;
    rval->m_SamplingSchemeRSpace = m_SamplingSchemeRSpace;
    rval->m_ODFOrder = m_ODFOrder;
    
    rval->m_DWISamples = m_DWISamples;
    rval->m_ODFSamples = m_ODFSamples;
    rval->m_EAPSamples = m_EAPSamples;

    return loPtr;
    }

  /** Sampling Scheme in q-space for non-negativity constraint  */
  SamplingSchemeQSpacePointer m_SamplingSchemeQSpace;
  /** Sampling Scheme in r-space for non-negativity constraint  */
  SamplingSchemeRSpacePointer m_SamplingSchemeRSpace;

  int m_ODFOrder;

  VectorPointer m_DWISamples;
  VectorPointer m_ODFSamples;
  VectorPointer m_EAPSamples;

private :
  DiffusionModelGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

#endif 
