/**
 *       @file  itkBasisMatrixGenerator.hxx
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

#ifndef __itkBasisMatrixGenerator_hxx
#define __itkBasisMatrixGenerator_hxx

#include "itkBasisMatrixGenerator.h"

namespace itk
{

template<typename TElement>
BasisMatrixGenerator<TElement>
::BasisMatrixGenerator() : Superclass(),
  m_QBasisMatrixForDWI(new MatrixType()), 
  m_RBasisMatrixForEAP(new MatrixType()), 
  m_RBasisMatrixForODF(new MatrixType()) 
{
  m_MD0 = 0.7e-3;
  m_ODFOrder = 2;
  m_OutputType = DWI;

  m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
  m_SamplingSchemeRSpace = SamplingSchemeQSpaceType::New();
}

template<typename TElement>
typename LightObject::Pointer
BasisMatrixGenerator<TElement>
::InternalClone() const
{
  utlShowPosition(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_MD0 = m_MD0;
  rval->m_ODFOrder = m_ODFOrder;
  rval->m_OutputType = m_OutputType;
  
  // NOTE: shared_ptr is thread safe, if the data is read only, thus do not need to copy the data block
  rval->m_QBasisMatrixForDWI = m_QBasisMatrixForDWI;

  rval->m_RBasisMatrixForEAP = m_RBasisMatrixForEAP;
  rval->m_RBasisMatrixForODF = m_RBasisMatrixForODF;

  rval->m_SamplingSchemeQSpace = m_SamplingSchemeQSpace;
  rval->m_SamplingSchemeRSpace = m_SamplingSchemeRSpace;

  return loPtr;
}

template<typename TElement>
void
BasisMatrixGenerator<TElement>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar2(true, m_MD0, m_ODFOrder, os<<indent);
  PrintEnum3(true, m_OutputType, DWI, EAP, ODF, os <<indent);
  
  os << indent << "m_SamplingSchemeQSpace: " << m_SamplingSchemeQSpace << std::endl;
  utl::PrintUtlMatrix(*m_QBasisMatrixForDWI, "m_QBasisMatrixForDWI", " ", os<<indent);

  os << indent << "m_SamplingSchemeRSpace: " << m_SamplingSchemeRSpace << std::endl;
  utl::PrintUtlMatrix(*m_RBasisMatrixForEAP, "m_RBasisMatrixForEAP", " ", os<<indent);
  utl::PrintUtlMatrix(*m_RBasisMatrixForODF, "m_RBasisMatrixForODF", " ", os<<indent);
}

template<typename TElement>
void
BasisMatrixGenerator<TElement>
::VerifyInputParameters() const
{
  if (m_OutputType==DWI)
    {
    MatrixPointer qOrientations = m_SamplingSchemeQSpace->GetOrientationsCartesian();
    STDVectorPointer bVector = m_SamplingSchemeQSpace->GetBVector();
    STDVectorPointer qVector = m_SamplingSchemeQSpace->GetRadiusVector();
    utlGlobalException(bVector->size()==0 && qVector->size()==0, "no qVector nor bVector for DWI output");
    utlGlobalException(qOrientations->Rows()==0, "no qOrientations for DWI output");
    utlSAGlobalException(qVector->size()>0 && qVector->size()!=qOrientations->Rows())(qVector->size())(qOrientations->Rows()).msg("sizes are different");
    utlSAGlobalException(bVector->size()>0 && bVector->size()!=qOrientations->Rows())(bVector->size())(qOrientations->Rows()).msg("sizes are different");
    }
  else
    {
    MatrixPointer rOrientations = m_SamplingSchemeRSpace->GetOrientationsCartesian();
    STDVectorPointer rVector = m_SamplingSchemeRSpace->GetRadiusVector();
    utlGlobalException(rOrientations->Rows()==0, "no rOrientations for ODF or EAP output");
    utlGlobalException(m_OutputType==EAP && rVector->size()==0, "no m_RVector for EAP output");
    utlSAGlobalException(rVector->size()>0 && rVector->size()!=rOrientations->Rows())(rVector->size())(rOrientations->Rows()).msg("sizes are different");
    }
}

}


#endif 

