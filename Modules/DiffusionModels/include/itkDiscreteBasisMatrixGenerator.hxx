/**
 *       @file  itkDiscreteBasisMatrixGenerator.hxx
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

#ifndef __itkDiscreteBasisMatrixGenerator_hxx
#define __itkDiscreteBasisMatrixGenerator_hxx

#include "itkDiscreteBasisMatrixGenerator.h"

namespace itk
{

template<typename TElement>
DiscreteBasisMatrixGenerator<TElement>
::DiscreteBasisMatrixGenerator() : Superclass(),
  m_BasisOrientations(new MatrixType())
{
  m_UseIsotropicTerm = false;
}

template<typename TElement>
typename LightObject::Pointer
DiscreteBasisMatrixGenerator<TElement>
::InternalClone() const
{
  utlShowPosition(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_UseIsotropicTerm = m_UseIsotropicTerm;
  rval->m_BasisOrientations = m_BasisOrientations;
  return loPtr;
}

template<typename TElement>
void
DiscreteBasisMatrixGenerator<TElement>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar1(true, m_UseIsotropicTerm, os<<indent);
  
  utl::PrintUtlMatrix(*m_BasisOrientations, "m_BasisOrientations", " ", os<<indent);
}

template<typename TElement>
void
DiscreteBasisMatrixGenerator<TElement>
::VerifyInputParameters() const
{
  Superclass::VerifyInputParameters();
  utlGlobalException(m_BasisOrientations->Rows()==0, "no m_BasisOrientations");
}

}


#endif 



