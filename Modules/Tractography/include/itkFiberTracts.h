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


template< typename TValueType=double >
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

  typedef TValueType ValueType;

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

  FiberPointer GetFiber(const int index) const
    {
    return m_Fibers->GetElement(index);
    }

  void PrintFibersHeader(std::ostream & os=std::cout, Indent indent=0) const
    {
    itk::PrintTractVisHeader(*m_Header, os <<indent);
    }

protected:
  FiberTracts(): m_Header(new HeaderType())
    {
    }
  ~FiberTracts(){}

  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  typename LightObject::Pointer InternalClone() const;

  HeaderPointer m_Header;

  FiberContainerPointer m_Fibers = FiberContainerType::New();


private:
  FiberTracts(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberTracts.hxx"
#endif

#endif 
