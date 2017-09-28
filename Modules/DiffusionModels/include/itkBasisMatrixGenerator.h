/**
 *       @file  itkBasisMatrixGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-18-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkBasisMatrixGenerator_h
#define __itkBasisMatrixGenerator_h

#include "itkObject.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <tr1/memory>
#include "itkSamplingSchemeQSpace.h"

namespace itk
{

/** \class BasisMatrixGenerator
 *  \brief Generate basis matrix.
 *
 * \ingroup DiffusionModels
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */

template <typename TElement = double>
class ITK_EXPORT BasisMatrixGenerator: public Object
{
public:
  /** Standard class typedefs. */
  typedef BasisMatrixGenerator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Standard part of every itk Object. */
  itkTypeMacro(BasisMatrixGenerator, Object);
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Save the template parameters. */
  typedef utl::NDArray<double,2>                  MatrixType;
  typedef utl::NDArray<double,1>                  VectorType;
  typedef utl_shared_ptr<MatrixType>              MatrixPointer;
  typedef utl_shared_ptr<VectorType>              VectorPointer;
  typedef std::vector<double>                     STDVectorType;
  typedef utl_shared_ptr<STDVectorType >          STDVectorPointer;
  
  typedef SamplingSchemeQSpace<double>                       SamplingSchemeQSpaceType;
  typedef typename SamplingSchemeQSpaceType::Pointer         SamplingSchemeQSpacePointer;
  
  typedef SamplingScheme3D<double>                           SamplingSchemeRSpaceType;
  typedef typename SamplingSchemeRSpaceType::Pointer         SamplingSchemeRSpacePointer;
  
  itkSetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);
  itkGetObjectMacro(SamplingSchemeQSpace, SamplingSchemeQSpaceType);

  itkSetObjectMacro(SamplingSchemeRSpace, SamplingSchemeRSpaceType);
  itkGetObjectMacro(SamplingSchemeRSpace, SamplingSchemeRSpaceType);

  typedef enum 
    {
    DWI=0,  
    EAP,
    ODF
    } OutputType;
  
  
  // itkSetMacro(QBasisMatrix, MatrixPointer);
  itkGetMacro(QBasisMatrixForDWI, MatrixPointer);
  
  // itkSetMacro(RBasisMatrix, MatrixPointer);
  itkGetMacro(RBasisMatrixForEAP, MatrixPointer);
  itkGetMacro(RBasisMatrixForODF, MatrixPointer);
  
  itkSetMacro(MD0, double);
  itkGetMacro(MD0, double);
  itkSetMacro(ODFOrder, int);
  itkGetMacro(ODFOrder, int);
  
  itkSetMacro(OutputType, OutputType);
  itkGetMacro(OutputType, OutputType);
  
  virtual void ComputeBasisMatrix()
    {
    }
  
  virtual void Flip(const int flipx, const int flipy, const int flipz)
    {
    }
  
  // [>* Generate basis in Q space for DWI <]
  // virtual MatrixPointer ComputeQBasisMatrixForDWI()
  //   {return MatrixPointer();}
  // [>* Generate basis in R space for EAP<]
  // virtual MatrixPointer ComputeRBasisMatrixForEAP()
  //   {return MatrixPointer();}
  // [>* Generate basis in R space for ODF<]
  // virtual MatrixPointer ComputeRBasisMatrixForODF()
  //   {return MatrixPointer();}

  virtual int GetNumberOfSamples() const
    {
    return m_OutputType==DWI? m_SamplingSchemeQSpace->GetNumberOfSamples() : m_SamplingSchemeRSpace->GetNumberOfSamples();
    }
  virtual int GetNumberOfBasis() const
    {
    return -1;
    }
  virtual MatrixPointer GetBasisMatrix() const
    {
    MatrixPointer mat;
    if (m_OutputType==DWI) 
      mat = m_QBasisMatrixForDWI;
    else if (m_OutputType==EAP)
      mat = m_RBasisMatrixForEAP;
    else if (m_OutputType==ODF)
      mat = m_RBasisMatrixForODF;
    return mat;
    }
  
protected:
  BasisMatrixGenerator();
  virtual ~BasisMatrixGenerator() 
    {}
  
  virtual void VerifyInputParameters() const; 

  void PrintSelf(std::ostream& os, Indent indent) const;
  typename LightObject::Pointer InternalClone() const;
  
  /** Sampling Scheme in q-space  */
  SamplingSchemeQSpacePointer m_SamplingSchemeQSpace;

  /** Sampling Scheme in r-space  */
  SamplingSchemeRSpacePointer m_SamplingSchemeRSpace;

  MatrixPointer    m_QBasisMatrixForDWI;  
  
  MatrixPointer    m_RBasisMatrixForEAP;  
  MatrixPointer    m_RBasisMatrixForODF;  
  
  /** typical MD value for typical scale */
  double m_MD0;
  
  /** for ODF  */
  int m_ODFOrder;

  OutputType m_OutputType;

private:
  BasisMatrixGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_BasisMatrixGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT BasisMatrixGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef BasisMatrixGenerator< ITK_TEMPLATE_2 x > BasisMatrixGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkBasisMatrixGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkBasisMatrixGenerator_hxx)
# include "itkBasisMatrixGenerator.hxx"
#endif


#endif 
