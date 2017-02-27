/**
 *       @file  itkFunctorHashTable.h
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

#ifndef __itkFunctorHashTable_h
#define __itkFunctorHashTable_h

#include "itkFunctorTableBase.h"
#include "utlSTDHeaders.h"

namespace itk
{

/**
 *   \class   FunctorHashTable
 *   \brief   use FunctorHashTable to accelerate evaluation of functions.
 *
 *   \ingroup Math
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TFunctor, class TParameters, class TFunctorValue = double, class THash=utl_hash<TParameters> >
class ITK_EXPORT FunctorHashTable
  : public FunctorTableBase<TFunctor, TParameters, TFunctorValue>
{
public:
  /** Standard class typedefs. */
  typedef FunctorHashTable         Self;
  typedef FunctorTableBase<TFunctor, TParameters, TFunctorValue>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( FunctorHashTable, FunctorTableBase );

  typedef typename Superclass::ParametersType        ParametersType; 
  typedef typename Superclass::FunctorType           FunctorType;
  typedef typename Superclass::FunctorValueType      FunctorValueType;

  typedef typename Superclass::STDVectorType         STDVectorType;
  typedef typename Superclass::STDVectorPointer      STDVectorPointer;

  typedef THash                                                       HashType;
  typedef utl_unordered_map<ParametersType,FunctorValueType, HashType >         HashTableType;
  typedef utl_shared_ptr<HashTableType >                              HashTablePointer; 
  typedef typename  HashTableType::iterator                           HashTableIterator;

  
  // itkSetMacro(Hash, HashTablePointer);
  itkGetMacro(Hash, HashTablePointer);
  
  unsigned long GetTableSize() const
    {
    return m_Hash->size();
    }
  
  bool IsTableBuilt() const 
    {
    return m_Hash->size()>0;
    }

  FunctorValueType GetFunctionValue ( const ParametersType& param)
    {
    HashTableIterator iter = m_Hash->find(param);
    if ( iter!= m_Hash->end() ) 
      {
      // if (this->GetDebug())
      //   std::cout << "parameter = " << param << ", is in the table" << std::endl << std::flush;
      return iter->second;
      }
    else
      {
      // if (this->GetDebug())
      //   std::cout << "parameter = " << param << ", is not in the table" << std::endl << std::flush;
      FunctorValueType val = m_Functor(param);
      m_Hash->insert(std::pair<ParametersType, FunctorValueType> (param, val));
      // (*m_Hash)[param]=val;
      return val;
      }
    }
  
protected:
  FunctorHashTable() : Superclass(), 
    m_Hash(new HashTableType())
    {
    }
  virtual ~FunctorHashTable() {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    }
  
  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();
    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_Hash = m_Hash;
    return loPtr;
    }
  
  FunctorType m_Functor;

  HashTablePointer m_Hash;

  
private:
  FunctorHashTable(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented
};

}

#endif 

