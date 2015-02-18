/**
 *       @file  itkSolverBase.hxx
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



#ifndef __itkSolverBase_hxx
#define __itkSolverBase_hxx

#include "itkSolverBase.h"
#include "utl.h"


namespace itk
{

//---------------------------------------------------------
template <class TPrecision>
SolverBase<TPrecision>
::SolverBase() : Superclass()
{
}

template < class TPrecision >
void
SolverBase<TPrecision>
::Solve ( const VectorType& xInitial) 
{
  this->VerifyInputs();
}

template < class TPrecision >
void
SolverBase<TPrecision>
::Initialize ( const VectorType& xInitial) 
{
  int N = GetXDimension();

  if (xInitial.Size()==0)
    {
    m_x.ReSize(N);
    m_x.Fill(0.0);
    }
  else
    {
    utlException(xInitial.Size()!=N, "wrong size of xInitial. xInitial.Size()="<<xInitial.Size());
    m_x = xInitial;
    }
}

template < class TPrecision >
void
SolverBase<TPrecision>
::Clear ( ) 
{
  m_x.Clear();
}

template < class TPrecision >
typename LightObject::Pointer
SolverBase<TPrecision>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_x = m_x;
  rval->SetDebug(this->GetDebug());
  return loPtr;
}

template <class TPrecision>
void
SolverBase<TPrecision>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  utl::PrintUtlVector(this->m_x, "m_x", " ", os << indent);
}

} // end namespace itk

#endif

