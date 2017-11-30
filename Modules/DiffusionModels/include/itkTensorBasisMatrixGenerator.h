/**
 *       @file  itkTensorBasisMatrixGenerator.h
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

#ifndef __itkTensorBasisMatrixGenerator_h
#define __itkTensorBasisMatrixGenerator_h

#include "itkDiscreteBasisMatrixGenerator.h"
#include "itkDiffusionTensor.h"

namespace itk
{

/** \class TensorBasisMatrixGenerator
 *  \brief Generate basis matrix.
 *
 *  \ingroup DiffusionModels
 *  \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */

template <typename TElement = double>
class ITK_EXPORT TensorBasisMatrixGenerator 
  : public DiscreteBasisMatrixGenerator<TElement>
{
public:
  /** Standard class typedefs. */
  typedef TensorBasisMatrixGenerator Self;
  typedef DiscreteBasisMatrixGenerator<TElement> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(TensorBasisMatrixGenerator, DiscreteBasisMatrixGenerator);

  /** Save the template parameters. */
  typedef typename Superclass::MatrixType         MatrixType;
  typedef typename Superclass::VectorType         VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType      STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;
  
  typedef  DiffusionTensor<double>  TensorType;
  
  /** tensor parameters */
  itkSetMacro(EigenValue1, double);  
  itkGetMacro(EigenValue1, double);  
  itkSetMacro(EigenValue2, double);  
  itkGetMacro(EigenValue2, double);  
  itkSetMacro(EigenValue3, double);  
  itkGetMacro(EigenValue3, double);  

  itkSetMacro(EigenValueISO, double);  
  itkGetMacro(EigenValueISO, double);  
  
  void ComputeBasisMatrix() ITK_OVERRIDE;

  // [>* Generate basis in Q space for DWI<]
  // virtual MatrixPointer ComputeQBasisMatrixForDWI();
  // [>* Generate basis in R space for EAP<]
  // virtual MatrixPointer ComputeRBasisMatrixForEAP();
  // [>* Generate basis in R space for ODF<]
  // virtual MatrixPointer ComputeRBasisMatrixForODF();
  
protected:
  TensorBasisMatrixGenerator();
  virtual ~TensorBasisMatrixGenerator() 
    {}
  
  virtual void VerifyInputParameters() const ITK_OVERRIDE; 

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;

  double m_EigenValue1;
  double m_EigenValue2;
  double m_EigenValue3;
  
  double m_EigenValueISO;

private:
  TensorBasisMatrixGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

// Define instantiation macro for this template.
#define ITK_TEMPLATE_TensorBasisMatrixGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT TensorBasisMatrixGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef TensorBasisMatrixGenerator< ITK_TEMPLATE_2 x > TensorBasisMatrixGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkTensorBasisMatrixGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkTensorBasisMatrixGenerator_hxx)
# include "itkTensorBasisMatrixGenerator.hxx"
#endif



#endif 

