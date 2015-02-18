/**
 *       @file  itkIterativeSolverBase.hxx
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



#ifndef __itkIterativeSolverBase_hxx
#define __itkIterativeSolverBase_hxx

#include "itkIterativeSolverBase.h"
#include "utl.h"


namespace itk
{

//---------------------------------------------------------
template <class TPrecision>
IterativeSolverBase<TPrecision>
::IterativeSolverBase() : Superclass()
{
  m_EPSOfVaribles = 1e-4;
  m_MaxNumberOfIterations = 1000;
  m_NumberOfIterations = 0;
  m_MinRelativeChangeOfCostFunction = 1e-2;
  m_MinRelativeChangeOfPrimalResidual = 1e-2;
  m_MinRelativeChangeOfDualResidual = 1e-2;
  m_UpdateInformation = NONE;
  m_NumberOfChangeLessThanThreshold=0;
  // m_MinRelativeChangeOfSeparateVariable = 1e-2;
}

template < class TPrecision >
void
IterativeSolverBase<TPrecision>
::Solve ( const VectorType& xInitial) 
{
  this->VerifyInputs();
  Initialize(xInitial);
}

template < class TPrecision >
void
IterativeSolverBase<TPrecision>
::Initialize ( const VectorType& xInitial) 
{
  Superclass::Initialize(xInitial);
  this->m_CostFunction.clear();
  this->m_DifferenceNormOfPrimalResidual.clear();
  this->m_DifferenceNormOfDualResidual.clear();

  this->m_EPSOfPrimalResidual.clear();
  this->m_EPSOfDualResidual.clear();
  m_NumberOfIterations = 0;
  m_NumberOfChangeLessThanThreshold=0;
}

template < class TPrecision >
void
IterativeSolverBase<TPrecision>
::Clear ( ) 
{
  Superclass::Clear();
  this->m_CostFunction.clear();
  this->m_DifferenceNormOfPrimalResidual.clear();
  this->m_DifferenceNormOfDualResidual.clear();

  this->m_EPSOfPrimalResidual.clear();
  this->m_EPSOfDualResidual.clear();
  m_NumberOfIterations = 0;
  m_UpdateInformation = NONE;
}

template < class TPrecision >
typename LightObject::Pointer
IterativeSolverBase<TPrecision>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();
  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_CostFunction = m_CostFunction;
  rval->m_DifferenceNormOfPrimalResidual = m_DifferenceNormOfPrimalResidual;
  rval->m_DifferenceNormOfDualResidual = m_DifferenceNormOfDualResidual;
  rval->m_EPSOfPrimalResidual = m_EPSOfPrimalResidual;
  rval->m_EPSOfDualResidual = m_EPSOfDualResidual;
  rval->m_MaxNumberOfIterations = m_MaxNumberOfIterations;
  rval->m_NumberOfIterations = m_NumberOfIterations;
  rval->m_MinRelativeChangeOfCostFunction = m_MinRelativeChangeOfCostFunction;
  rval->m_MinRelativeChangeOfPrimalResidual = m_MinRelativeChangeOfPrimalResidual;
  rval->m_MinRelativeChangeOfDualResidual = m_MinRelativeChangeOfDualResidual;
  rval->m_EPSOfVaribles = m_EPSOfVaribles;
  rval->m_UpdateInformation = m_UpdateInformation;
  rval->m_NumberOfChangeLessThanThreshold = m_NumberOfChangeLessThanThreshold;
  return loPtr;
}

template <class TPrecision>
void
IterativeSolverBase<TPrecision>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar4(true, m_MaxNumberOfIterations, m_NumberOfIterations, m_EPSOfVaribles, m_UpdateInformation, os<<indent);
  // PrintVar4(true, m_MinRelativeChangeOfCostFunction, m_MinRelativeChangeOfPrimalResidual, 
  // m_MinRelativeChangeOfDualResidual, m_MinRelativeChangeOfSeparateVariable, os<<indent);
  PrintVar4(true, m_MinRelativeChangeOfCostFunction, m_MinRelativeChangeOfPrimalResidual, 
  m_MinRelativeChangeOfDualResidual, m_NumberOfChangeLessThanThreshold, os<<indent);
  utl::PrintVector(m_CostFunction, "m_CostFunction", " ", os<<indent);
  utl::PrintVector(m_DifferenceNormOfPrimalResidual, "m_DifferenceNormOfPrimalResidual", " ", os<<indent);
  utl::PrintVector(m_DifferenceNormOfDualResidual, "m_DifferenceNormOfDualResidual", " ", os<<indent);
  utl::PrintVector(m_EPSOfPrimalResidual, "m_EPSOfPrimalResidual", " ", os<<indent);
  utl::PrintVector(m_EPSOfDualResidual, "m_EPSOfDualResidual", " ", os<<indent);
}

} // end namespace itk

#endif

