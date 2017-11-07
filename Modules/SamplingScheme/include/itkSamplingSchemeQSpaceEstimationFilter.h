/**
 *       @file  itkSamplingSchemeQSpaceEstimationFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-14-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpaceEstimationFilter_h
#define __itkSamplingSchemeQSpaceEstimationFilter_h

#include "itkLightProcessObject.h"
#include "itkSamplingScheme3D.h"
#include "utlITK.h"

namespace itk
{

/**
 *   \class   SamplingSchemeQSpaceEstimationFilter
 *   \brief   base class for the filters to estimate the sampling scheme in Q-space
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup SamplingScheme
 */
template <class TSamplingType>
class ITK_EXPORT SamplingSchemeQSpaceEstimationFilter : public LightProcessObject
{

public:
  /** Standard class typedefs. */
  typedef SamplingSchemeQSpaceEstimationFilter                      Self;
  typedef LightProcessObject                                  Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;
  
  typedef TSamplingType                                       SamplingType;
  typedef typename SamplingType::Pointer                      SamplingPointer;
  typedef typename SamplingType::ConstPointer                 SamplingConstPointer;
  
  typedef typename SamplingType::ValueType                    ValueType;
  typedef typename SamplingType::MatrixType                   MatrixType;
  typedef typename SamplingType::MatrixPointer                MatrixPointer;
  typedef typename SamplingType::STDVectorType                STDVectorType;
  typedef typename SamplingType::STDVectorPointer             STDVectorPointer;
  typedef typename SamplingType::IndexVectorType              IndexVectorType;
  typedef typename SamplingType::Index2DVectorType            Index2DVectorType;
  typedef typename SamplingType::Index2DVectorPointer         Index2DVectorPointer;

  /** Standard Macros */
  itkTypeMacro(SamplingSchemeQSpaceEstimationFilter, LightProcessObject);
  itkNewMacro(Self);
  
  typedef enum 
    {
    DISTANCE=0, 
    ELECTROSTATIC
    } CriteriaType;

  itkSetMacro(CriteriaType, CriteriaType);
  itkGetMacro(CriteriaType, CriteriaType);

  itkSetMacro(ElectrostaticOrder, double);
  itkGetMacro(ElectrostaticOrder, double);
  
  itkSetNDebugMacro(NumbersInShell, IndexVectorType);
  itkGetMacro(NumbersInShell, IndexVectorType);
  
  itkSetMacro(WeightForSingleShell, double);
  itkGetMacro(WeightForSingleShell, double);
  
  itkSetObjectMacro(InitialOrientations, SamplingType);
  itkGetObjectMacro(InitialOrientations, SamplingType);
  
  itkGetObjectMacro(OutputOrientations, SamplingType);
  itkGetConstObjectMacro(OutputOrientations, SamplingType);

  virtual void GenerateData() ITK_OVERRIDE{}

  bool IsSetInitialization() const
    {
    return m_InitialOrientations && m_InitialOrientations->GetNumberOfSamples()>0;
    }

protected:
  
  virtual void Initialization(){}

  SamplingSchemeQSpaceEstimationFilter()
    {
    m_WeightForSingleShell = 0.5;
    m_ElectrostaticOrder = 2.0;
    m_CriteriaType = DISTANCE;
    m_InitialOrientations = NULL;
    m_OutputOrientations = SamplingType::New();
    }

  ~SamplingSchemeQSpaceEstimationFilter(){}
  
  double m_ElectrostaticOrder;

  /** It contains number of samples in each shell  */
  IndexVectorType m_NumbersInShell;
  
  /** weight in the cost function for single shell term. It is in \f$(0,1)\f$ */
  double m_WeightForSingleShell;

  /** Initalization of the first several points  */
  SamplingPointer  m_InitialOrientations;
  
  SamplingPointer  m_OutputOrientations;
  // Index2DVectorPointer  m_Indices;

  CriteriaType m_CriteriaType;

private:
  SamplingSchemeQSpaceEstimationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented
};

}


#endif 

