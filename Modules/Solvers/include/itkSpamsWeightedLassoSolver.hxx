/**
 *       @file  itkSpamsWeightedLassoSolver.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "02-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkSpamsWeightedLassoSolver_hxx
#define __itkSpamsWeightedLassoSolver_hxx

#include "itkSpamsWeightedLassoSolver.h"
#include "utlCore.h"
// #include "utlITKSpams.h"
#include "utlSpams.h"

namespace itk
{

template < class TPrecision >
SpamsWeightedLassoSolver<TPrecision>
::SpamsWeightedLassoSolver () : Superclass(), 
  m_A(new MatrixType()),
  m_B(new MatrixType()),
  m_W(new MatrixType()),
  m_X(new MatrixType()),
  m_As(new SpamsMatrixType()),
  m_Bs(new SpamsMatrixType()),
  m_Ws(new SpamsMatrixType()),
  m_Xs(new SpamsSpMatrixType())
{
  m_ConstraintType=PENALTY;
  m_Positive = false;
  m_Lambda = 1.0;
  m_NumberOfThreads = -1;
}

template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::SetA (const MatrixPointer& mat) 
{
  itkDebugMacro("setting A to " << *mat);
  if ( this->m_A != mat )
    {
    this->m_A = mat;
    m_As = SpamsMatrixPointer(new SpamsMatrixType());
    spams::UtlMatrixToMatrix(*m_A, *m_As);
    this->Modified();
    }
}

template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::SetB (const MatrixPointer& b) 
{
  itkDebugMacro("setting B to " << *b);
  if ( this->m_B != b )
    {
    m_B = b;
    m_Bs = SpamsMatrixPointer(new SpamsMatrixType());
    spams::UtlMatrixToMatrix(*m_B, *m_Bs);
    this->Modified();
    }
}

template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::Setb (const VectorPointer& b) 
{
  itkDebugMacro("setting B to " << *b);
  MatrixPointer B( new MatrixType(b->Size(),1));
  B->SetColumn(0, *b);
  if ( this->m_B != B )
    {
    m_B = B;
    m_Bs = SpamsMatrixPointer(new SpamsMatrixType());
    spams::UtlMatrixToMatrix(*m_B, *m_Bs);
    this->Modified();
    }
}

template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::SetW (const MatrixPointer& w) 
{
  itkDebugMacro("setting W to " << *w);
  if ( this->m_W != w )
    {
    this->m_W = w;
    m_Ws = SpamsMatrixPointer(new SpamsMatrixType());
    spams::UtlMatrixToMatrix(*m_W, *m_Ws);
    this->Modified();
    }
}

template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::Setw (const VectorPointer& w) 
{
  itkDebugMacro("setting W to " << *w);
  MatrixPointer W( new MatrixType(w->Size(),1) );
  W->SetColumn(0, *w);
  if ( this->m_W != W )
    {
    m_W = W;
    m_Ws = SpamsMatrixPointer(new SpamsMatrixType());
    spams::UtlMatrixToMatrix(*m_W, *m_Ws);
    this->Modified();
    }
}


template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::VerifyInputs () const
{
  Superclass::VerifyInputs();
  int N = this->GetXDimension();
  int M = this->GetXNumber();
  utlException(m_W->Columns()!=M, "wrong size of m_W!, m_W->Columns()="<<m_W->Columns()<<", M="<<M);
  utlException(m_A->Columns()!=N, "wrong size of m_A!, m_A->Columns()="<<m_A->Columns()<<", N="<<N);
  utlException(m_A->Rows()!=m_B->Rows(), "wrong rows of m_A!, m_A->Rows()="<<m_A->Rows()<<", m_B->Rows()="<<m_B->Rows());
}

template < class TPrecision >
typename LightObject::Pointer
SpamsWeightedLassoSolver<TPrecision>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();
  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_A = m_A;
  rval->m_As = m_As;
  rval->m_B = m_B;
  rval->m_Bs = m_Bs;
  rval->m_W = m_W;
  rval->m_Ws = m_Ws;

  rval->m_X = m_X;
  rval->m_Xs = m_Xs;
  // if (m_Xs->m()*m_Xs->n()>0)
  //   rval->m_Xs->copy(*m_Xs);
  // else
  //   rval->m_Xs= SpamsSpMatrixPointer(new SpamsSpMatrixType());

  rval->m_ConstraintType = m_ConstraintType;
  rval->m_Lambda = m_Lambda;
  rval->m_Positive = m_Positive;
  rval->m_NumberOfThreads = m_NumberOfThreads;
  
  return loPtr;
}


template < class TPrecision >
void
SpamsWeightedLassoSolver<TPrecision>
::Solve (const VectorType& ) 
{
  utlShowPosition(this->GetDebug());
  this->VerifyInputs();
  // this->Initialize();
  // SpamsSpMatrixType xs;
  int N = this->GetXDimension();
  int M = this->GetXNumber();
  spams::constraint_type mode;
  if (m_ConstraintType==L1CONS)
    {
    mode = spams::L1COEFFS;
    utlGlobalException(true, "TODO");
    }
  else if (m_ConstraintType==L2CONS)
    {
    mode = spams::L2ERROR;
    utlGlobalException(true, "TODO");
    }
  else if (m_ConstraintType==PENALTY)
    {
    mode = spams::PENALTY;
    // m_Bs.print("m_Bs");
    // m_As.print("m_As");
    // m_Ws.print("m_Ws");
    // std::cout << "m_Lambda = " << m_Lambda << std::endl << std::flush;
    // std::cout << "utl::min(m_A.Columns(), m_A.Rows()) = " << utl::min<int>(m_A.Columns(), m_A.Rows()) << std::endl << std::flush;
    // std::cout << "mode = " << mode << std::endl << std::flush;
    // std::cout << "m_Positive = " << m_Positive << std::endl << std::flush;
    m_Xs= SpamsSpMatrixPointer(new SpamsSpMatrixType()); 
    spams::lassoWeight<double>(*m_Bs,*m_As,*m_Ws,*m_Xs,utl::min(m_A->Columns(), m_A->Rows()), 0.5*m_Lambda, mode,m_Positive,m_NumberOfThreads);
    }
  else 
    utlGlobalException(true, "wrong m_ConstraintType");
  m_X= MatrixPointer(new MatrixType()); 
  spams::SpMatrixToUtlMatrix(*m_Xs, *m_X);
  if (M==1)
    this->m_x = m_X->GetColumn(0);
  // utl::PrintUtlMatrix(m_X, "m_X");
  // utl::PrintUtlVector(this->m_x, "m_x");
}

template < class TPrecision >
typename SpamsWeightedLassoSolver<TPrecision>::ValueType
SpamsWeightedLassoSolver<TPrecision>
::EvaluateCostFunctionInColumn (const VectorType& x, const int col) const 
{
  utlException(x.Size()==0, "need to give a vector");
  const VectorType* xx = &x;
  VectorType e = (*m_A) * (*xx)-m_B->GetColumn(col);
  ValueType func = e.GetSquaredTwoNorm() + m_Lambda* utl::ToVector<double>(m_W->GetColumn(col) % (*xx))->GetOneNorm();
  return func;
}

template < class TPrecision >
typename SpamsWeightedLassoSolver<TPrecision>::ValueType
SpamsWeightedLassoSolver<TPrecision>
::EvaluateCostFunction (const MatrixType& x) const 
{
  ValueType func=0;
  const MatrixType* xx = (x.Size()!=0? (&x) : (this->m_X.get()));
  for ( int i = 0; i < xx->Columns(); i += 1 ) 
    {
    VectorType vec;
    vec = xx->GetColumn(i);
    func += EvaluateCostFunctionInColumn(vec, i);
    }
  return func;
}


template <class TPrecision>
void
SpamsWeightedLassoSolver<TPrecision>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar3(true, m_ConstraintType, m_Positive, m_Lambda, os<<indent);
  utl::PrintUtlMatrix(*m_A, "m_A", " ", os<<indent);
  utl::PrintUtlMatrix(*m_B, "m_B", " ", os<<indent);
  utl::PrintUtlMatrix(*m_W, "m_W", " ", os<<indent);
  utl::PrintUtlMatrix(*m_X, "m_X", " ", os<<indent);

  // m_As.print("m_As");
  // m_Bs.print("m_Bs");
  // m_Ws.print("m_Ws");
  // m_Xs.print("m_Xs");
}

}

#endif 



