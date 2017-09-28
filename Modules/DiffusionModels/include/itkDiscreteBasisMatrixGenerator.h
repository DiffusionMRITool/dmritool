/**
 *       @file  itkDiscreteBasisMatrixGenerator.h
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

#ifndef __itkDiscreteBasisMatrixGenerator_h
#define __itkDiscreteBasisMatrixGenerator_h


#include "itkBasisMatrixGenerator.h"

namespace itk
{

/** \class DiscreteBasisMatrixGenerator
 *  \brief Generate basis matrix.
 *
 * \ingroup DiffusionModels
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */

template <typename TElement = double>
class ITK_EXPORT DiscreteBasisMatrixGenerator
  : public BasisMatrixGenerator<TElement>
{
public:
  /** Standard class typedefs. */
  typedef DiscreteBasisMatrixGenerator Self;
  typedef BasisMatrixGenerator<TElement> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Standard part of every itk Object. */
  itkTypeMacro(DiscreteBasisMatrixGenerator, BasisMatrixGenerator);
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Save the template parameters. */
  typedef typename Superclass::MatrixType         MatrixType;
  typedef typename Superclass::VectorType         VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType      STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;
    
  itkSetMacro(BasisOrientations, MatrixPointer);
  itkGetMacro(BasisOrientations, MatrixPointer);
  
  int GetNumberOfBasis() const override
    {
    int numberOfBasis = this->m_BasisOrientations->Rows();
    if (m_UseIsotropicTerm)
      numberOfBasis += 1;
    return numberOfBasis;
    }
  
  /** Some parameters */
  itkSetMacro(UseIsotropicTerm, bool);
  itkGetMacro(UseIsotropicTerm, bool);
  itkBooleanMacro(UseIsotropicTerm);
  
  void Flip(const int flipx, const int flipy, const int flipz) override
    {
    std::vector<int> flipVec(3);
    flipVec[0]=flipx, flipVec[1]=flipy, flipVec[2]=flipz;
    MatrixPointer mat (new MatrixType());
    *mat = utl::FlipOrientations(*m_BasisOrientations, flipVec);
    m_BasisOrientations = mat;
    }
  
protected:
  DiscreteBasisMatrixGenerator();
  virtual ~DiscreteBasisMatrixGenerator() 
    {}
  
  virtual void VerifyInputParameters() const;

  void PrintSelf(std::ostream& os, Indent indent) const;
  typename LightObject::Pointer InternalClone() const;

  bool          m_UseIsotropicTerm;
  
  MatrixPointer    m_BasisOrientations;

private:
  DiscreteBasisMatrixGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_DiscreteBasisMatrixGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT DiscreteBasisMatrixGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef DiscreteBasisMatrixGenerator< ITK_TEMPLATE_2 x > DiscreteBasisMatrixGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDiscreteBasisMatrixGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDiscreteBasisMatrixGenerator_hxx)
# include "itkDiscreteBasisMatrixGenerator.hxx"
#endif


#endif 
