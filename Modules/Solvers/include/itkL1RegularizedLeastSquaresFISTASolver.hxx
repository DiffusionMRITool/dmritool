/**
 *       @file  itkL1RegularizedLeastSquaresFISTASolver.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-03-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkL1RegularizedLeastSquaresFISTASolver_hxx
#define __itkL1RegularizedLeastSquaresFISTASolver_hxx

#include "itkL1RegularizedLeastSquaresFISTASolver.h"
#include "utl.h"
#include "utlVNLBlas.h"

namespace itk
{

template < class TPrecision >
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::L1RegularizedLeastSquaresFISTASolver () : Superclass(),
  m_A(new MatrixType()),
  m_b(new VectorType()),
  m_w(new VectorType()),
  m_At(new MatrixType()),
  m_AtA(new MatrixType()),
  m_Atb(new VectorType()),
  m_xOld(new VectorType())
{
  m_Step=-1;
  m_L2Solver = L2SolverType::New();
  m_UseL2SolverForInitialization = false;
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::SetA (const MatrixPointer& mat) 
{
  itkDebugMacro("setting A to " << *mat);
  if ( *this->m_A != *mat )
    {
    m_A=MatrixPointer(new MatrixType());
    m_At=MatrixPointer(new MatrixType());
    m_AtA=MatrixPointer(new MatrixType());
    *this->m_A = *mat;
    // utl::MatrixCopy(*mat, *this->m_A, 1.0, 'N');
    this->Modified();
    this->m_A->GetTranspose(*this->m_At);
    // utl::MatrixCopy(*this->m_A, *m_At, 1.0, 'T');
    utl::ProductUtlXtX(*m_A, *m_AtA);
    m_Step = 0.5/m_AtA->GetTwoNorm();
    if (m_b->Size()>0)
      {
      utlException(m_A->Rows()!=m_b->Size(), "wrong size of m_A");
      m_Atb=VectorPointer(new VectorType());
      utl::ProductUtlMv(*m_At, *m_b, *m_Atb);
      }
    }
  if (m_UseL2SolverForInitialization)
    m_L2Solver->SetA(m_A);
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::Setb (const VectorPointer& b) 
{
  itkDebugMacro("setting b to " << *b);
  if ( *this->m_b != *b )
    {
    m_b=VectorPointer(new VectorType());
    *this->m_b = *b;
    this->Modified();
    if (m_A->Size()>0)
      {
      utlException(m_At->Columns()!=m_b->Size(), "wrong size of m_A");
      m_Atb=VectorPointer(new VectorType());
      utl::ProductUtlMv(*m_At, *m_b, *m_Atb);
      }
    }
  if (m_UseL2SolverForInitialization)
    m_L2Solver->Setb(m_b);
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::Setw (const VectorPointer& w) 
{
  itkDebugMacro("setting w to " << *w);
  if ( this->m_w != w )
    {
    this->m_w = w;
    this->Modified();
    }
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::SetwForInitialization (const VectorPointer& w) 
{
  int N = w->Size();
  MatrixPointer lambda(new MatrixType(N, N));
  lambda->Fill(0.0);
  lambda->SetDiagonal(*w);
  utlException(!m_UseL2SolverForInitialization, "need to set m_UseL2SolverForInitialization");
  m_L2Solver->SetLambda(lambda);
  this->Modified();
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::VerifyInputs () const
{
  Superclass::VerifyInputs();
  int N = this->GetXDimension();
  utlException(m_w->Size()!=N, "wrong size of m_w!, m_w->Size()="<<m_w->Size()<<", N="<<N);
  utlException(m_A->Columns()!=N, "wrong size of m_A!, m_A->Columns()="<<m_A->Columns()<<", N="<<N);
  utlException(m_A->Rows()!=m_b->Size(), "wrong rows of m_A!, m_A->Rows()="<<m_A->Rows()<<", m_b->Size()="<<m_b->Size());
}

template < class TPrecision >
typename LightObject::Pointer
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();
  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_A = m_A;
  rval->m_b = m_b;
  rval->m_w = m_w;
  rval->m_At = m_At;
  rval->m_AtA = m_AtA;
  rval->m_Atb = m_Atb;
  rval->m_Step = m_Step;
  rval->m_UseL2SolverForInitialization = m_UseL2SolverForInitialization;
  rval->m_L2Solver = m_L2Solver->Clone();
  return loPtr;
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::Initialize ( const VectorType& xInitial) 
{
  utlShowPosition(this->GetDebug());
  Superclass::Initialize(xInitial);
  if (xInitial.Size()==0)
    {
    utlException(!m_UseL2SolverForInitialization, "need to set m_UseL2SolverForInitialization");
    if (m_L2Solver->GetLambda()->Size()==0)
      SetwForInitialization(this->m_w);
    m_L2Solver->Solve();
    this->m_x = m_L2Solver->Getx();
    // m_L2Solver->Print(std::cout<<"m_L2Solver : ");
    }
  if (this->GetDebug())
    utl::PrintUtlVector(this->m_x, "m_x initialization");
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::Solve ( const VectorType& xInitial) 
{
  utlShowPosition(this->GetDebug());
  this->VerifyInputs();
  Initialize(xInitial);
  Iterate();
  // EndSolve();
}

template < class TPrecision >
typename L1RegularizedLeastSquaresFISTASolver<TPrecision>::ValueType
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::EvaluateCostFunction (const VectorType& x) const 
{
  const VectorType* xx = (x.Size()!=0? (&x) : (&this->m_x));
  VectorPointer e(new VectorType());
  utl::ProductUtlMv(*m_A, *xx, *e);
  utl::vSub(e->Size(), e->GetData(), m_b->GetData(), e->GetData());
  // *e -= *m_b;
  double eNorm = utl::cblas_nrm2(e->Size(), e->GetData(), 1);
  VectorType tmp(e->Size());
  utl::vMul(e->Size(), m_w->GetData(), (double*)xx->GetData(), tmp.GetData());
  ValueType func = eNorm*eNorm + utl::cblas_asum(e->Size(), tmp.GetData(),1);
  // ValueType func = e->squared_magnitude() + element_product(*m_w, *xx).one_norm();
  return func;
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::HistoryUpdateAndConvergenceCheck ( ) 
{
  ValueType fValue = EvaluateCostFunction(), changePercentage=0, changePercentage_x=0;
  this->m_CostFunction.push_back(fValue);
  int size = this->m_CostFunction.size();
  VectorType tmp(this->m_x.Size());
  if (size>=2)
    {
    changePercentage = (this->m_CostFunction[size-2] - this->m_CostFunction[size-1])/this->m_CostFunction[size-2];
    double xOldNorm = utl::cblas_nrm2(this->m_xOld->Size(), this->m_xOld->GetData(), 1);
    if (xOldNorm>0)
      {
      utl::vSub(this->m_x.Size(), this->m_x.GetData(), m_xOld->GetData(), tmp.GetData());
      changePercentage_x = utl::cblas_nrm2(tmp.Size(), tmp.GetData(), 1) / xOldNorm;
      }
    if (changePercentage <= this->m_MinRelativeChangeOfCostFunction && changePercentage>=0 && changePercentage_x<=this->m_MinRelativeChangeOfPrimalResidual)
      {
      this->m_UpdateInformation = Self::CONTINUE;
      this->m_NumberOfChangeLessThanThreshold++;
      }
    else
      {
      this->m_UpdateInformation = Self::CONTINUE;
      this->m_NumberOfChangeLessThanThreshold = 0;
      }
    if (this->m_NumberOfChangeLessThanThreshold == 3)
      this->m_UpdateInformation = Self::STOP_MIN_CHANGE;
    }
  else
    this->m_UpdateInformation = Self::CONTINUE;
  // if (this->GetDebug())
  //   {
  //   // utl::PrintUtlVector(this->m_x, "m_x");
  //   utlPrintVar4(true, fValue, changePercentage, changePercentage_x, this->m_MinRelativeChangeOfCostFunction);
  //   std::cout << "m_NumberOfChangeLessThanThreshold = " << this->m_NumberOfChangeLessThanThreshold << std::endl << std::flush;
  //   }
}

template < class TPrecision >
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::Iterate () 
{
  utlShowPosition(this->GetDebug());

  VectorType y=this->m_x, xg=this->m_x;
  double t = 1.0, t_new=1.0, func=-1.0;
  VectorType w_new = (*m_w) % m_Step, tmp(xg.Size());
  this->m_CostFunction.push_back(EvaluateCostFunction());

  int N = GetXDimension();
  this->m_NumberOfIterations=0;
  int num_change_less = 0;
  m_xOld->ReSize(this->m_x.Size());
  while ( this->m_NumberOfIterations <= this->m_MaxNumberOfIterations ) 
    {

    utl::cblas_copy(this->m_x.Size(), this->m_x.GetData(), 1, m_xOld->GetData(), 1);
    if (this->GetDebug())
      {
      std::cout << "iter = " << this->m_NumberOfIterations << std::endl << std::flush;
      }
  
    // FISTA
    utl::ProductUtlMv(*m_AtA, y, tmp);
    utl::vSub(tmp.Size(), tmp.GetData(), m_Atb->GetData(), tmp.GetData());
    utl::cblas_scal<double>(tmp.Size(), m_Step*2, tmp.GetData(), 1);
    utl::vSub(y.Size(), y.GetData(), tmp.GetData(), xg.GetData());
    // xg = y - m_Step*2*(tmp - *m_Atb);
    for ( int j = 0; j < N; j += 1 ) 
      {
      xg[j] = xg[j]>w_new[j] ? (xg[j]-w_new[j]) : (xg[j]<-w_new[j] ? (xg[j]+w_new[j]) : 0 ); 
      }
    t_new = 0.5+0.5*std::sqrt(1+4*t*t);

    utl::vSub(xg.Size(), xg.GetData(), this->m_x.GetData(), tmp.GetData());
    utl::cblas_scal<double>(xg.Size(), (t-1)/t_new, tmp.GetData(), 1);
    utl::vAdd(xg.Size(), xg.GetData(), tmp.GetData(), y.GetData());
    // y = xg + (t-1)/t_new * (xg-this->m_x);

    // // ISTA
    // xg = m_x - m_Step*2*(m_AtA*m_x-m_Atb);
    // for ( int j = 0; j < N; j += 1 ) 
    //   
    //   xg[j] = xg[j]>w_new[j] ? (xg[j]-w_new[j]) : (xg[j]<-w_new[j] ? (xg[j]+w_new[j]) : 0 ); 
    //   }

    // update
    // this->m_x = xg; 
    utl::cblas_copy(this->m_x.Size(), xg.GetData(), 1, this->m_x.GetData(), 1);
    t = t_new;

    this->m_NumberOfIterations++;

    HistoryUpdateAndConvergenceCheck();
    if (this->m_UpdateInformation==Self::STOP_MIN_CHANGE)
      break;
    }
}

template <class TPrecision>
void
L1RegularizedLeastSquaresFISTASolver<TPrecision>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  utl::PrintUtlMatrix(*m_A, "m_A", " ", os<<indent);
  utl::PrintUtlVector(*m_b, "m_b", " ", os<<indent);
  utl::PrintUtlVector(*m_w, "m_w", " ", os<<indent);
  PrintVar2(true, m_Step, m_UseL2SolverForInitialization, os<<indent);
  if (m_L2Solver)
    m_L2Solver->Print(os<<indent<< "m_L2Solver for initialization = ");
}

}


#endif 




