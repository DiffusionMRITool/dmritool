/**
 *       @file  itkUnaryFunctorLookUpTable.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-04-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkUnaryFunctorLookUpTable_h
#define __itkUnaryFunctorLookUpTable_h

#include "itkFunctorTableBase.h"

namespace itk
{

/**
 *   \class   UnaryFunctorLookUpTable
 *   \brief   use UnaryFunctorLookUpTable to accelerate evaluation of functions.
 *
 *   \ingroup Math
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class TFunctor >
class ITK_EXPORT UnaryFunctorLookUpTable
  : public FunctorTableBase<TFunctor, double, double>
{
public:
  /** Standard class typedefs. */
  typedef UnaryFunctorLookUpTable         Self;
  typedef  FunctorTableBase<TFunctor, double, double> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( UnaryFunctorLookUpTable, FunctorTableBase );
  
  typedef typename Superclass::ParametersType        ParametersType; 
  typedef typename Superclass::FunctorType           FunctorType;
  typedef typename Superclass::FunctorValueType      FunctorValueType;

  typedef typename Superclass::STDVectorType         STDVectorType;
  typedef typename Superclass::STDVectorPointer      STDVectorPointer;

  itkSetMacro(VariableMax, double);
  itkGetMacro(VariableMax, double);
  
  itkSetMacro(VariableMin, double);
  itkGetMacro(VariableMin, double);
  
  itkSetMacro(NumberOfBins, int);
  itkGetMacro(NumberOfBins, int);
  
  itkGetMacro(Table, STDVectorPointer);

  void Initialize()
    {
    BuildTable();
    }
  
  void BuildTable()
    {
    utlGlobalException(m_NumberOfBins<=0, "need to compute m_NumberOfBins first");
    utlGlobalException(m_VariableMax<=m_VariableMin, "m_VariableMax should be larger than m_VariableMin");
    m_Table = STDVectorPointer(new STDVectorType(m_NumberOfBins+1));
    m_Delta = (m_VariableMax-m_VariableMin)/(double)m_NumberOfBins;
    m_DeltaInv = 1.0/m_Delta;

    for (int i=0; i<= m_NumberOfBins; i++)   
      (*m_Table)[i] = this->m_Functor( m_Delta*i + m_VariableMin );
    }

  unsigned long GetTableSize() const
    {
    return m_Table->size();
    }

  bool IsTableBuilt() const 
    {
    return m_NumberOfBins>0 && m_Table->size()==m_NumberOfBins+1;
    }

  FunctorValueType GetFunctionValue ( const ParametersType& var)
    {
    utlGlobalException(m_Table->size()==0, "need to compute m_Table first");
    if (var<=m_VariableMin)
      return (*m_Table)[0];
    if (var>=m_VariableMax)
      return m_Table->back();

    double xDouble = (var-m_VariableMin)*m_DeltaInv;
    int x = (int) std::floor(xDouble);
    return (*m_Table)[x] + ((*m_Table)[x+1]-(*m_Table)[x])*(xDouble-x);
    }

protected:
  UnaryFunctorLookUpTable() : Superclass(), 
    m_Table(new STDVectorType())
    {
    m_VariableMax = 0;
    m_VariableMin = 0;
    m_Delta = 0;
    m_DeltaInv = 0;
    m_NumberOfBins = -1;
    }
  ~UnaryFunctorLookUpTable() {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    PrintVar4(true, m_VariableMax, m_VariableMin, m_NumberOfBins, m_Delta, os<<indent);
    for (int i=0; i<= m_NumberOfBins; i++)   
      utlPrintVar3(true, i, m_Delta*i+m_VariableMin, (*m_Table)[i]);
    }
  
  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();
    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_VariableMax = m_VariableMax;
    rval->m_VariableMin = m_VariableMin;
    rval->m_Delta = m_Delta;
    rval->m_DeltaInv = m_DeltaInv;
    rval->m_NumberOfBins = m_NumberOfBins;

    rval->m_Table = m_Table;
    return loPtr;
    }

  double m_VariableMax;
  double m_VariableMin;
  double m_Delta;
  double m_DeltaInv;

  /** number of bins between m_VariableMin and m_VariableMax  */
  int m_NumberOfBins;

  /** the size of m_Table is m_NumberOfBins+1  */
  STDVectorPointer m_Table;
  

private:
  UnaryFunctorLookUpTable(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

};


}



#endif 

