/**
 *       @file  itkSolverBase.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-26-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSolverBase_h
#define __itkSolverBase_h

#include <itkObject.h>
#include <itkObjectFactory.h>
#include "vnl/vnl_matrix.h"
#include "utlSTDHeaders.h"
#include "utlNDArray.h"


namespace itk
{

/**
 *   \class   SolverBase
 *   \brief   Base class for some optimization solvers using primal-dual updates
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup OptimizationSolver
 */
template <class TPrecision>
class ITK_EXPORT SolverBase : public Object
{
public:
  /** Standard class typedefs. */
  typedef SolverBase           Self;
  typedef Object               Superclass;
  typedef SmartPointer<Self>   Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SolverBase, Object);

  typedef TPrecision                                  ValueType;
  typedef utl::NDArray<ValueType,2>                   MatrixType;
  typedef utl::NDArray<ValueType,1>                   VectorType;
  typedef utl_shared_ptr<MatrixType>                  MatrixPointer;
  typedef utl_shared_ptr<VectorType>                  VectorPointer;
  typedef std::vector<ValueType>                      ValueContainerType;
  typedef utl_shared_ptr<std::vector<ValueType> >     ValueContainerPointer;
  
  virtual int GetXDimension() const 
    {return 0;}
  itkGetConstReferenceMacro(x, VectorType);

  virtual void VerifyInputs() const {}
  virtual void Initialize(const VectorType& xInitial=VectorType()); 
  virtual void EndSolve() {}

  /** if x is not set, evaluate the cost function for m_x  */
  virtual ValueType EvaluateCostFunction(const VectorType& x=VectorType()) const {return ValueType(0.0);}
  virtual ValueType EvaluateCostFunction(const MatrixType& x=MatrixType()) const {return ValueType(0.0);}
  /** if x is not set, evaluate the gradients of the cost function for m_x  */
  virtual VectorType EvaluateGradientOfCostFunction(const VectorType& x) const {return VectorType(); }
  
  virtual void Solve(const VectorType& xInitial=VectorType()); 
  
  virtual void Solve(const MatrixType& xInitial=MatrixType())
    {
    this->VerifyInputs();
    }

  virtual void Clear();

protected:
  SolverBase();
  virtual ~SolverBase() {};
  
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual typename LightObject::Pointer InternalClone() const;
  
  /** Nx1 vector primal variable */
  VectorType m_x;

private:
  SolverBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

// Define instantiation macro for this template.
#define ITK_TEMPLATE_SolverBase(_, EXPORT, TypeX, TypeY)         \
  namespace itk                                                 \
  {                                                             \
  _( 1 ( class EXPORT SolverBase< ITK_TEMPLATE_1 TypeX > ) )     \
  namespace Templates                                           \
  {                                                             \
  typedef SolverBase< ITK_TEMPLATE_1 TypeX > SolverBase##TypeY; \
  }                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkSolverBase+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSolverBase_hxx)
#include "itkSolverBase.hxx"
#endif

#endif 
