/**
 *       @file  itkFiber.hxx
 *      @brief  
 *     Created  "07-24-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkFiber_hxx
#define __itkFiber_hxx

#include "itkFiber.h"

namespace itk
{

template< class TValue >
typename LightObject::Pointer
Fiber< TValue >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  
  rval->m_Tract = m_Tract;
  rval->m_Properties = m_Properties;
  rval->m_Scalars = m_Scalars;
  return loPtr;
}

template< class TValue >
void
Fiber< TValue >
::PrintSelf(std::ostream & os, Indent indent) const 
{
  // Superclass::PrintSelf(os, indent);

  int n_s = GetNumberOfScalarsPerPoint();
  int n_p = GetNumberOfProperties();
  PrintVar(true, os<<indent, n_s, n_p);

  // m_Tract->Print(os, indent);
  auto points = m_Tract->GetVertexList();
  os << indent << "points (" << (points) << "): [";
  for ( int i = 0; i < points->Size(); ++i ) 
    {
    os << points->GetElement(i);
    os << ((i==points->Size()-1) ? "]\n" : ", ");
    }

  if (n_p>0)
    utl::PrintVector(*m_Properties, "m_Properties", " ", os <<indent);
  if (n_s>0)
    {
    os << indent << "m_Scalars : [";
    for ( int i = 0; i < m_Scalars->size(); ++i ) 
      {
      os << (*m_Scalars)[i];
      os << ((i==m_Scalars->size()-1) ? "]\n" : ", ");
      }
    }
}

}


#endif 

