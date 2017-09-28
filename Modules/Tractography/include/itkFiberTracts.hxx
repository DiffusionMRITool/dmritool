/**
 *       @file  itkFiberTracts.hxx
 *      @brief  
 *     Created  "07-24-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkFiberTracts_hxx
#define __itkFiberTracts_hxx

#include "itkFiberTracts.h"

namespace itk
{

template< class TValue >
typename LightObject::Pointer
FiberTracts< TValue >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  
  rval->m_Header = m_Header;
  rval->m_Fibers = m_Fibers;
  return loPtr;
}

template< class TValue >
void
FiberTracts< TValue >
::PrintSelf(std::ostream & os, Indent indent) const 
{
  Superclass::PrintSelf(os, indent);
  PrintFibersHeader(os, indent);

  utlSAGlobalException(m_Header->n_count!=m_Fibers->Size())
    (m_Header->n_count)(m_Fibers->Size()).msg("the number of tracts are not consistent in m_Header and m_Fibers");

  os << indent << std::endl;
  os << indent << "The number of fibers: " << GetNumberOfFibers() << std::endl;
  for ( int i = 0; i < m_Fibers->Size(); ++i ) 
    {
    FiberPointer fiber = (*m_Fibers)[i];
    os << indent << std::endl;
    os << indent << "Fiber #" << i << " :    number of points = " << fiber->GetNumberOfPoints() << std::endl;
    fiber->Print(os, indent);
    }
  os << indent<< std::endl;
}

}


#endif 
