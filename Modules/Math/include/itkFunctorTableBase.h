/**
 *       @file  itkFunctorTableBase.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkFunctorTableBase_h
#define __itkFunctorTableBase_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "utlSTDHeaders.h"
#include "utlCoreMacro.h"

namespace itk
{

/**
 *   \class   FunctorTableBase
 *   \brief   use FunctorTableBase to accelerate evaluation of functions.
 *
 *   \ingroup Math
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TFunctor, class TParameters, class TFunctorValue  >
class ITK_EXPORT FunctorTableBase
  : public Object
{
public:
  /** Standard class typedefs. */
  typedef FunctorTableBase         Self;
  typedef Object  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( FunctorTableBase, Object );
  
  typedef TParameters                                ParametersType; 
  typedef TFunctor                                   FunctorType;
  typedef TFunctorValue                              FunctorValueType;

  typedef std::vector<double>                        STDVectorType;
  typedef utl_shared_ptr<STDVectorType >             STDVectorPointer;

  
  /** Get the functor object.  The functor is returned by reference.
   * (Functors do not have to derive from itk::LightObject, so they do
   * not necessarily have a reference count. So we cannot return a
   * SmartPointer.) */
  FunctorType &       GetFunctor() { return m_Functor; }
  const FunctorType & GetFunctor() const { return m_Functor; }

  /** Set the functor object.  This replaces the current Functor with a
   * copy of the specified Functor. This allows the user to specify a
   * functor that has ivars set differently than the default functor.
   * This method requires an operator!=() be defined on the functor
   * (or the compiler's default implementation of operator!=() being
   * appropriate). */
  void SetFunctor(const FunctorType & functor)
  {
    if ( m_Functor != functor )
      {
      m_Functor = functor;
      this->Modified();
      }
  }

  virtual void Initialize()
    {
    }

  void BuildTable()
    {
    }

  unsigned long GetTableSize() const
    {
    return -1;
    }

  bool IsTableBuilt() const 
    {
    return false;
    }

  /**  \note virtual function is a little bit slower. Thus it is not efficient if a virtual function is called too many times   */
  FunctorValueType GetFunctionValue ( const ParametersType& param)
    {
    return m_Functor(param);
    }

protected:
  FunctorTableBase() : Superclass() 
    {
    }
  virtual ~FunctorTableBase() {};

  // void PrintSelf(std::ostream& os, Indent indent) const
  //   {
  //   Superclass::PrintSelf(os, indent);
  //   }
  
  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();
    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_Functor = m_Functor;
    return loPtr;
    }
  
  FunctorType m_Functor;

private:
  FunctorTableBase(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

};


}


#endif 

