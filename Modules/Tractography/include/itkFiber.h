/**
 *       @file  itkFiber.h
 *      @brief  
 *     Created  "07-24-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkFiber_h
#define __itkFiber_h

#include <itkDataObject.h>
#include <itkObjectFactory.h>

#include "itkSlowPolyLineParametricPath.h"

#include "utlSmartAssert.h"
#include "utlCore.h"
#include "utlITKMacro.h"

namespace itk 
{

template< typename TValueType=double >
class Fiber: public DataObject
{
public:
  /** Standard class typedefs. */
  typedef Fiber                      Self;
  typedef DataObject                 Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;
  
  itkNewMacro(Self);

  itkTypeMacro(Fiber, DataObject );

  typedef TValueType ValueType;

  typedef SlowPolyLineParametricPath< 3 >       TractType;
  typedef typename TractType::Pointer           TractPointer;
  typedef typename TractType::VertexType        VertexType;

  typedef std::vector<ValueType>                STDVectorType;
  typedef std::shared_ptr<STDVectorType>        STDVectorPointer;
  typedef std::vector<std::vector<ValueType> >  STD2DVectorType;
  typedef std::shared_ptr<STD2DVectorType>      STD2DVectorPointer;

  int GetNumberOfPoints() const
    {
    return m_Tract->GetVertexList()->Size();
    }

  int GetNumberOfProperties() const
    {
    return m_Properties->size();
    }
  
  int GetNumberOfScalarsPerPoint() const
    {
    return m_Scalars->size()>0? (*m_Scalars)[0].size() : 0;
    }

  VertexType GetPoint(const int index) const
    {
    return m_Tract->GetVertexList()->GetElement(index);
    }
  
  Vector<double,3> GetDirection(const double pos) const
    {
    if (GetNumberOfPoints()==1)
      {
      // the fiber only has one point, then there is no direction for that point
      Vector<double,3> vec;
      vec.Fill(0.0);
      return vec;
      }
    else
      {
      // for fiber with more than 2 points
      return m_Tract->EvaluateDerivative(pos);
      }
    }

  itkSetGetMacro(Properties, STDVectorPointer);

  itkSetGetMacro(Scalars, STD2DVectorPointer);

  itkSetGetMacro(Tract, TractPointer);

protected:
  Fiber(): m_Properties(new STDVectorType()), m_Scalars(new STD2DVectorType())
    {}
  ~Fiber(){}

  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
  
  typename LightObject::Pointer InternalClone() const;

  TractPointer m_Tract = TractType::New();

  /** properties for a single tract  */
  STDVectorPointer  m_Properties;

  /** various scalars for points along a single tract  */
  STD2DVectorPointer m_Scalars;


private:
  Fiber(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiber.hxx"
#endif


#endif 
