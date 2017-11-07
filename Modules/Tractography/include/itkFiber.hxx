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
#include "utlDMRI.h"

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
typename Fiber<TValue>::Pointer
Fiber< TValue >
::DeepClone() const
{
  Pointer out = Self::New();
  *out->m_Scalars = *m_Scalars;
  *out->m_Properties = *m_Properties;

  auto tract = out->GetTract();
  auto vertexList = m_Tract->GetVertexList();
  int numPoints = GetNumberOfPoints();
  VertexType vertex;
  for ( int i = 0; i < numPoints; ++i ) 
    {
    vertex = (*vertexList)[i];
    tract->AddVertex(vertex);
    }

  return out;
}

template< class TValue >
double
Fiber< TValue >
::DistanceToPoint(double x, double y, double z) const
{
  VertexType vertex;
  double d2=0, d2Min=std::numeric_limits<double>::max();
  double dx, dy, dz;
  for ( int i = 0; i < GetNumberOfPoints(); ++i ) 
    {
    vertex = GetPoint(i);
    dx = x-vertex[0];
    dy = y-vertex[1];
    dz = z-vertex[2];
    d2 = dx*dx+dy*dy+dz*dz;
    if (d2 < d2Min)
      d2Min = d2;
    }
  return std::sqrt(d2Min);
}

template< class TValue >
std::vector<double>
Fiber< TValue >
::GetPointDistanceStats() const
{
  std::vector<double> distVec;
  VertexType p0, p1;
  for ( int i = 1; i < GetNumberOfPoints(); ++i ) 
    {
    p0 = GetPoint(i-1);
    p1 = GetPoint(i);
    distVec.push_back(p0.EuclideanDistanceTo(p1));
    }
  return utl::GetContainerStats(distVec.begin(), distVec.end());
}

template< class TValue >
void
Fiber< TValue >
::RemoveScalarsByName(const std::string& name, const std::vector<std::string>& nameVec)
{
  int numPoints = GetNumberOfPoints();
  for ( int i = 0; i < numPoints; ++i ) 
    {
    utl::RemoveScalarsByName((*m_Scalars)[i], nameVec, name);
    }
}

template< class TValue >
void
Fiber< TValue >
::RemovePropertiesByName(const std::string& name, const std::vector<std::string>& nameVec)
{
  utl::RemoveScalarsByName(*m_Properties, nameVec, name);
}

template< class TValue >
void
Fiber< TValue >
::PrintSelf(std::ostream & os, Indent indent) const 
{
  // Superclass::PrintSelf(os, indent);

  int dim_s = GetDimensionOfScalarsPerPoint();
  int dim_p = GetDimensionOfProperties();
  PrintVar(true, os<<indent, dim_s, dim_p);

  // m_Tract->Print(os, indent);
  auto points = m_Tract->GetVertexList();
  os << indent << "points (" << (points) << "): (numberOfPoints: "<<  points->size() << ") [";
  for ( int i = 0; i < points->Size(); ++i ) 
    {
    os << points->GetElement(i);
    os << ((i==points->Size()-1) ? "]\n" : ", ");
    }

  if (dim_p>0)
    utl::PrintVector(*m_Properties, "m_Properties", " ", os <<indent);
  if (dim_s>0)
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

