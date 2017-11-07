/**
 *       @file  itkFiberTracts.h
 *      @brief  
 *     Created  "07-22-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkFiberTracts_h
#define __itkFiberTracts_h

#include "itkFiber.h"

#include "itkTrackvisHeader.h"
#include "utlSmartAssert.h"
#include "utlCore.h"
#include "utlITKMacro.h"

namespace itk
{


template< typename TValue=double >
class FiberTracts: public DataObject
{
public:
  /** Standard class typedefs. */
  typedef FiberTracts                Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  
  itkNewMacro(Self);

  itkTypeMacro(FiberTracts, DataObject );

  typedef TValue ValueType;

  typedef TrackVisHeaderType  HeaderType;
  typedef std::shared_ptr<HeaderType>  HeaderPointer;
  
  typedef Fiber<ValueType>                   FiberType;
  typedef typename FiberType::Pointer        FiberPointer;
  typedef typename FiberType::TractType      TractType;
  typedef typename FiberType::TractPointer   TractPointer;
  typedef typename TractType::VertexType     VertexType;

  typedef VectorContainer< unsigned int, FiberPointer >  FiberContainerType;
  typedef typename FiberContainerType::Pointer FiberContainerPointer;
  
  itkSetGetMacro(Header, HeaderPointer);

  itkSetGetMacro(Fibers, FiberContainerPointer);

  int GetNumberOfFibers() const
    {
    return m_Fibers->Size();
    }

  int GetNumberOfPoints() const;

  FiberPointer GetFiber(const int index) const
    {
    return m_Fibers->GetElement(index);
    }

  typename FiberTracts<TValue>::Pointer DeepClone() const;

  void PrintFibersHeader(std::ostream & os=std::cout, Indent indent=0) const
    {
    itk::PrintTractVisHeader(*m_Header, os <<indent);
    }

  bool HasPropertyName(const std::string& name) const
    {
    return itk::HasPropertyName(*m_Header, name);
    }
  bool HasScalarName(const std::string& name) const
    {
    return itk::HasScalarName(*m_Header, name);
    }

  /** remove scalars by name  */
  void RemoveScalarsByName(const std::string& name);
  /** remove properties by name  */
  void RemovePropertiesByName(const std::string& name);

  /** Get the stats of the number of points in each tract, and the stats of distances between two points in each tract.  */
  void GetPointStats(std::vector<double>& numPointsStats, std::vector<double>& pointDistStats) const;

  /** Get fiber tracts based on indices  */
  Pointer SelectByIndicesOfTracts(const std::vector<int>& indices);
  
  /** Get fiber tracts across in a ball ROI (origin x,y,z and radius)  */
  Pointer SelectByBallROI(double x, double y, double z, double radius);

  /** Append a fiber. Should have the same number of scalars and the same number of properties.   */
  void AppendFiber(const FiberPointer fiber);

  /** Append fibers. Should have the same number of scalars and the same number of properties.   */
  void AppendFibers(const Pointer fibers);

protected:
  FiberTracts(): m_Header(new HeaderType())
    {
    }
  ~FiberTracts(){}

  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;

  HeaderPointer m_Header;

  FiberContainerPointer m_Fibers = FiberContainerType::New();


private:
  FiberTracts(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};


}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkFiberTracts_hxx)
#include "itkFiberTracts.hxx"
#endif

#endif 
