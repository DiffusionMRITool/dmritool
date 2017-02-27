/**
 *       @file  itkL2RegularizedLeastSquaresSolver.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-25-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkL2RegularizedLeastSquaresSolver_hxx
#define __itkL2RegularizedLeastSquaresSolver_hxx

#include "itkL2RegularizedLeastSquaresSolver.h"
#include "utl.h"
#include "utlVNLBlas.h"
#include "utlVNLLapack.h"


namespace itk
{

template <class TPrecision>
L2RegularizedLeastSquaresSolver<TPrecision>
::L2RegularizedLeastSquaresSolver() : Superclass(), 
  m_A(new MatrixType()),
  m_b(new VectorType()),
  m_Lambda(new MatrixType()),
  m_LS(new MatrixType())
{
  m_ConditionNumber = -1;
  // empty matrix is symmetric
  m_IsLambdaSymmetric = true;
}

template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::SetA(const MatrixPointer& mat)
{
  itkDebugMacro("setting A to " << *mat);
  if ( *this->m_A != *mat )
    {
    // NOTE: use value copy because mat can be changed outside, while m_LS can only be changed inside.
    m_A=MatrixPointer(new MatrixType());
    *this->m_A = *mat;
    this->Modified();
    m_LS = MatrixPointer(new MatrixType());
    m_ConditionNumber=-1;
    }
}

template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::SetLambda(const MatrixPointer& mat)
{
  itkDebugMacro("setting Lambda to " << *mat);
  if ( *this->m_Lambda != *mat )
    {
    this->m_Lambda=MatrixPointer(new MatrixType());
    *this->m_Lambda = *mat;
    this->Modified();
    m_IsLambdaSymmetric = m_Lambda->IsSymmetric();
    m_LS = MatrixPointer(new MatrixType());
    m_ConditionNumber=-1;
    }
}

template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::VerifyInputs() const
{
  Superclass::VerifyInputs();
  int N = GetXDimension();
  int M = m_A->Rows();
  utlGlobalException(M<=0 || N<=0, "need to set m_A" );
  utlGlobalException(M!=m_b->Size(), "wrong size of m_b! m_A->rows=()"<<m_A->Rows() <<", m_b->Size()="<<m_b->Size());
  if (m_Lambda->Size()>0)
    {
    utlGlobalException(m_Lambda->Rows()!=m_Lambda->Columns(), "m_Lambda needs to be square");
    utlGlobalException(m_Lambda->Rows()!=N, "wrong size of m_Lambda! m_Lambda->Rows()="<<m_Lambda->Rows() << ", N="<<N);
    }
}

template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::Clear()
{
  Superclass::Clear();
  ClearA();
  Clearb();
  ClearLambda();
}

/** m_LS is re-calcutated only if m_A or m_Lambda is re-set */
template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::Initialize(const VectorType& xInitial)
{
  Superclass::Initialize(xInitial);
  if (m_LS->Size()==0)
    {
    m_LS = MatrixPointer(new MatrixType());
    // utl::ProductVnlMtM(*m_A, *m_A, *m_LS);
    utl::ProductUtlXtX(*m_A, *m_LS);
    if (m_Lambda->Size()>0)
      *m_LS += *m_Lambda;
      // utl::vAdd(m_LS->Size(),m_LS->data_block(), m_Lambda->data_block(), m_LS->data_block());
    m_ConditionNumber = m_LS->GetInfNorm();
    MatrixPointer tmp(new MatrixType());
    if (m_IsLambdaSymmetric)
      *tmp = m_LS->PInverseSymmericMatrix();
    else
      *tmp = m_LS->PInverseMatrix();
    m_ConditionNumber *= tmp->GetInfNorm();
    utl::ProductUtlMMt(*tmp, *m_A, *m_LS);
    }
}

template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::Solve(const VectorType& xInitial)
{
  VerifyInputs();
  Initialize(xInitial);
  utl::ProductUtlMv(*m_LS, *m_b, this->m_x);
}

template <class TPrecision>
typename L2RegularizedLeastSquaresSolver<TPrecision>::ValueType
L2RegularizedLeastSquaresSolver<TPrecision>
::EvaluateCostFunction(const VectorType& x) const
{
  const VectorType* xx = (x.Size()!=0? (&x) : (&this->m_x));
  VectorType tmp;
  utl::ProductUtlMv(*m_A, *xx, tmp);
  ValueType cost = utl::ToVector<double>(tmp - *m_b)->GetSquaredTwoNorm();
  if (m_Lambda->Size()>0)
    {
    utl::ProductUtlvM(*xx, *m_Lambda, tmp);
    cost += tmp.InnerProduct(*xx);
    }
  return cost;
}

template < class TPrecision >
typename LightObject::Pointer
L2RegularizedLeastSquaresSolver<TPrecision>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();
  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  *rval->m_A = *m_A;
  *rval->m_b = *m_b;
  *rval->m_Lambda = *m_Lambda;
  *rval->m_LS = *m_LS;
  rval->m_IsLambdaSymmetric = m_IsLambdaSymmetric;
  rval->m_ConditionNumber = m_ConditionNumber;
  return loPtr;
}

template <class TPrecision>
void
L2RegularizedLeastSquaresSolver<TPrecision>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  utl::PrintUtlMatrix(*m_A, "m_A", " ", os << indent);
  utl::PrintUtlVector(*m_b, "m_b", " ", os << indent);
  utl::PrintUtlMatrix(*m_Lambda, "m_Lambda", " ", os << indent);
  if (m_LS->Size()>0)
    utl::PrintUtlMatrix(*m_LS, "m_LS", " ", os << indent);
  os << indent << "m_ConditionNumber = " << m_ConditionNumber << std::endl << std::flush;
  os << indent << "m_IsLambdaSymmetric = " << m_IsLambdaSymmetric << std::endl << std::flush;
}

}

#endif 
