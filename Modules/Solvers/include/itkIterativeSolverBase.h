/**
 *       @file  itkIterativeSolverBase.h
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

#ifndef __itkIterativeSolverBase_h
#define __itkIterativeSolverBase_h

#include "itkSolverBase.h"

namespace itk
{

/**
 *   \class   IterativeSolverBase
 *   \brief   Base class for some optimization solvers using primal-dual updates
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup OptimizationSolver
 */
template <class TPrecision>
class ITK_EXPORT IterativeSolverBase : public SolverBase<TPrecision>
{
public:
  /** Standard class typedefs. */
  typedef IterativeSolverBase           Self;
  typedef SolverBase<TPrecision>        Superclass;
  typedef SmartPointer<Self>   Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(IterativeSolverBase, SolverBase);

  typedef TPrecision ValueType;
  typedef typename Superclass::MatrixType MatrixType;
  typedef typename Superclass::VectorType VectorType;

  typedef std::vector<ValueType> ValueContainerType;
  
  typedef enum 
    {
    NONE=0,    
    STOP_MIN_CHANGE,    
    STOP_MAX_NUM_ITERATION,
    CONTINUE,
    RESTART               
    } UpdateInfomationType;
  
  // itkSetMacro(CostFunction, ValueContainerType);
  itkGetConstReferenceMacro(CostFunction, ValueContainerType);
  // itkSetMacro(DifferenceNormOfPrimalResidual, ValueContainerType);
  itkGetConstReferenceMacro(DifferenceNormOfPrimalResidual, ValueContainerType);
  // itkSetMacro(DifferenceNormOfDualResidual, ValueContainerType);
  itkGetConstReferenceMacro(DifferenceNormOfDualResidual, ValueContainerType);
  // itkSetMacro(EPSOfPrimalResidual, ValueContainerType);
  itkGetConstReferenceMacro(EPSOfPrimalResidual, ValueContainerType);
  // itkSetMacro(EPSOfDualResidual, ValueContainerType);
  itkGetConstReferenceMacro(EPSOfDualResidual, ValueContainerType);
  // itkSetMacro(DifferenceNormOfSeparateVariable, ValueContainerType);
  // itkGetConstReferenceMacro(DifferenceNormOfSeparateVariable, ValueContainerType);

  itkSetMacro(MaxNumberOfIterations, int);
  itkGetConstReferenceMacro(MaxNumberOfIterations, int);
  // itkSetMacro(NumberOfIterations, int);
  itkGetConstReferenceMacro(NumberOfIterations, int);
  itkSetMacro(MinRelativeChangeOfCostFunction, ValueType);
  itkGetConstReferenceMacro(MinRelativeChangeOfCostFunction, ValueType);
  itkSetMacro(MinRelativeChangeOfPrimalResidual, ValueType);
  itkGetConstReferenceMacro(MinRelativeChangeOfPrimalResidual, ValueType);
  itkSetMacro(MinRelativeChangeOfDualResidual, ValueType);
  itkGetConstReferenceMacro(MinRelativeChangeOfDualResidual, ValueType);
  // itkSetMacro(MinRelativeChangeOfSeparateVarible, ValueType);
  // itkGetConstReferenceMacro(MinRelativeChangeOfSeparateVariable, ValueType);

  itkSetMacro(EPSOfVaribles, ValueType);
  itkGetConstReferenceMacro(EPSOfVaribles, ValueType);
  
  itkSetMacro(NumberOfChangeLessThanThreshold, int);
  itkGetConstReferenceMacro(NumberOfChangeLessThanThreshold, int);

  void Initialize(const VectorType& xInitial=VectorType()) ITK_OVERRIDE; 
  virtual void Iterate() {}
  
  /** Update history information and monitor stop conditions  */
  virtual void HistoryUpdateAndConvergenceCheck() {}
  
  virtual void Solve(const VectorType& xInitial=VectorType()) ITK_OVERRIDE; 

  virtual void Clear() ITK_OVERRIDE;

protected:
  IterativeSolverBase();
  virtual ~IterativeSolverBase() {};
  
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  virtual typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  
  ValueContainerType m_CostFunction;
  ValueContainerType m_DifferenceNormOfPrimalResidual;
  ValueContainerType m_DifferenceNormOfDualResidual;
  // ValueContainerType m_DifferenceNormOfSeparateVariable;
  
  ValueContainerType m_EPSOfPrimalResidual;
  ValueContainerType m_EPSOfDualResidual;

  int m_MaxNumberOfIterations;
  int m_NumberOfIterations;
  ValueType m_MinRelativeChangeOfCostFunction;
  ValueType m_MinRelativeChangeOfPrimalResidual;
  ValueType m_MinRelativeChangeOfDualResidual;
  // ValueType m_MinRelativeChangeOfSeparateVariable;

  ValueType m_EPSOfVaribles;

  UpdateInfomationType m_UpdateInformation;
  int m_NumberOfChangeLessThanThreshold;

private:
  IterativeSolverBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

// Define instantiation macro for this template.
#define ITK_TEMPLATE_IterativeSolverBase(_, EXPORT, TypeX, TypeY)         \
  namespace itk                                                 \
  {                                                             \
  _( 1 ( class EXPORT IterativeSolverBase< ITK_TEMPLATE_1 TypeX > ) )     \
  namespace Templates                                           \
  {                                                             \
  typedef IterativeSolverBase< ITK_TEMPLATE_1 TypeX > IterativeSolverBase##TypeY; \
  }                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkIterativeSolverBase+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkIterativeSolverBase_hxx)
#include "itkIterativeSolverBase.hxx"
#endif

#endif 
