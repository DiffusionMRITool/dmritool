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
typename FiberTracts<TValue>::Pointer
FiberTracts< TValue >
::DeepClone() const
{
  Pointer out = Self::New();
  itk::CopyTrackvisHeader(*m_Header, *out->m_Header);
  auto fibers = out->GetFibers();
  int num = GetNumberOfFibers();
  for ( int i = 0; i < num; ++i ) 
    {
    auto fiber = this->GetFiber(i)->DeepClone();
    fibers->InsertElement( fibers->Size(), fiber);
    }
  return out;
}

template< class TValue >
int
FiberTracts< TValue >
::GetNumberOfPoints() const
{
  int num=0;
  for ( int i = 0; i < GetNumberOfFibers(); ++i ) 
    num += GetFiber(i)->GetNumberOfPoints();
  return num;
}

template< class TValue >
void
FiberTracts< TValue >
::GetPointStats(std::vector<double>& numPointsStats, std::vector<double>& pointDistStats) const
{
  std::vector<double> numVec, distVec;
  for ( int i = 0; i < GetNumberOfFibers(); ++i ) 
    {
    auto fiber = GetFiber(i);
    numVec.push_back(fiber->GetNumberOfPoints());
    auto vec = fiber->GetPointDistanceStats();
    distVec.push_back(vec[2]); // mean dist
    }
  numPointsStats = utl::GetContainerStats(numVec.begin(), numVec.end());
  pointDistStats = utl::GetContainerStats(distVec.begin(), distVec.end());
}

template< class TValue >
void
FiberTracts< TValue >
::RemoveScalarsByName(const std::string& name)
{
  std::vector<std::string> scalarNames = utl::CovertChar2DArrayToStringArray(GetHeader()->scalar_name, 10);
  for ( int i = 0; i < GetNumberOfFibers(); ++i ) 
    this->GetFiber(i)->RemoveScalarsByName(name, scalarNames);

  itk::RemoveScalarName(*m_Header, name);
}

template< class TValue >
void
FiberTracts< TValue >
::RemovePropertiesByName(const std::string& name)
{
  std::vector<std::string> scalarNames = utl::CovertChar2DArrayToStringArray(GetHeader()->property_name, 10);
  for ( int i = 0; i < GetNumberOfFibers(); ++i ) 
    this->GetFiber(i)->RemovePropertiesByName(name, scalarNames);

  itk::RemovePropertyName(*m_Header, name);
}

template< class TValue >
typename FiberTracts<TValue>::Pointer
FiberTracts< TValue >
::SelectByIndicesOfTracts(const std::vector<int>& indices)
{
  Pointer result = Self::New();
  itk::CopyTrackvisHeader(*m_Header, *result->m_Header);
  result->m_Header->n_count = indices.size();

  auto fibers = result->GetFibers();
  for ( int i = 0; i < indices.size(); ++i ) 
    fibers->InsertElement(i, GetFiber(indices[i]));

  return result;
}

template< class TValue >
typename FiberTracts<TValue>::Pointer
FiberTracts< TValue >
::SelectByBallROI(double x, double y, double z, double radius)
{
  Pointer result = Self::New();
  itk::CopyTrackvisHeader(*m_Header, *result->m_Header);

  int num=0;
  auto fibers = result->GetFibers();
  for ( int i = 0; i < GetNumberOfFibers(); ++i ) 
    {
    auto fiber_i = GetFiber(i);
    double dist = fiber_i->DistanceToPoint(x,y,z);
    if (dist < radius)
      {
      fibers->InsertElement(fibers->Size(), fiber_i);
      num++;
      }
    }
  result->m_Header->n_count=num;

  return result;
}

template< class TValue >
void
FiberTracts< TValue >
::AppendFiber(const FiberPointer fiber)
{
  m_Fibers->InsertElement(m_Fibers->Size(), fiber);
  m_Header->n_count++;
}

template< class TValue >
void
FiberTracts< TValue >
::AppendFibers(const Pointer fibers)
{
  utlGlobalException(!itk::IsSameStructure(*m_Header, *fibers->m_Header), "two headers do not have the same structure. Cannot merge.");
  for ( int i = 0; i < fibers->GetNumberOfFibers(); ++i ) 
    AppendFiber(fibers->GetFiber(i));
}

template< class TValue >
void
FiberTracts< TValue >
::PrintSelf(std::ostream & os, Indent indent) const 
{
  Superclass::PrintSelf(os, indent);
  PrintFibersHeader(os, indent);

  std::vector<double> numPointsStats, pointDistStats;
  GetPointStats(numPointsStats, pointDistStats);

  utlSAGlobalException(m_Header->n_count!=m_Fibers->Size())
    (m_Header->n_count)(m_Fibers->Size()).msg("the number of tracts are not consistent in m_Header and m_Fibers");

  os << indent << std::endl;
  os << indent << "The number of fibers: " << GetNumberOfFibers() << std::endl;
  os << indent << "Total number of points: " << GetNumberOfPoints() << std::endl << std::flush;
  utl::PrintVector(numPointsStats, "stats of the number of points (min,max,mean,std): ", " ", os << indent, false);
  utl::PrintVector(pointDistStats, "stats of distances between points (min,max,mean,std): ", " ", os << indent, false);
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
