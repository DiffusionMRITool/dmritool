/**
 *       @file  itkSamplingSchemeQSpaceIncrementalEstimationFilter.h
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

#ifndef __itkSamplingSchemeQSpaceIncrementalEstimationFilter_h
#define __itkSamplingSchemeQSpaceIncrementalEstimationFilter_h

#include "itkSamplingSchemeQSpaceEstimationFilter.h"
#include "itkSamplingScheme3D.h"


namespace itk
{

/**
 *   \class   SamplingSchemeQSpaceIncrementalEstimationFilter
 *   \brief   incremental estimation of single/multi-shell orientations
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup SamplingScheme
 */
template <class TSamplingType>
class ITK_EXPORT SamplingSchemeQSpaceIncrementalEstimationFilter 
  : public SamplingSchemeQSpaceEstimationFilter<TSamplingType>
{

public:
  /** Standard class typedefs. */
  typedef SamplingSchemeQSpaceIncrementalEstimationFilter Self;
  typedef SamplingSchemeQSpaceEstimationFilter<TSamplingType>                 Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;
  
  typedef TSamplingType                                       SamplingType;
  typedef typename SamplingType::Pointer                      SamplingPointer;
  
  typedef typename Superclass::ValueType                    ValueType;
  typedef typename Superclass::MatrixType                   MatrixType;
  typedef typename Superclass::MatrixPointer                MatrixPointer;
  typedef typename Superclass::STDVectorType                STDVectorType;
  typedef typename Superclass::STDVectorPointer             STDVectorPointer;
  typedef typename Superclass::IndexVectorType              IndexVectorType;
  typedef typename Superclass::Index2DVectorType            Index2DVectorType;
  typedef typename Superclass::Index2DVectorPointer         Index2DVectorPointer;

  /** Standard Macros */
  itkTypeMacro(SamplingSchemeQSpaceIncrementalEstimationFilter, SamplingSchemeQSpaceEstimationFilter);
  itkNewMacro(Self);
  
  itkSetMacro(TessellationOrder, unsigned int);
  itkGetMacro(TessellationOrder, unsigned int);
  
  itkSetMacro(FineOrientations, MatrixPointer);
  itkGetMacro(FineOrientations, MatrixPointer);

  void GenerateData();

  // static void IndicesOfMaximalDistance(const MatrixType& mat, unsigned int& r1, unsigned int& r2);

protected:
  SamplingSchemeQSpaceIncrementalEstimationFilter();
  ~SamplingSchemeQSpaceIncrementalEstimationFilter(){}
  
  void Initialization();
  
  /** the order of tessellation for the orignal fine mesh  */
  unsigned int m_TessellationOrder;
  
  /** It is the input data or generated from m_TessellationOrder  */
  MatrixPointer    m_FineOrientations;

private:
  SamplingSchemeQSpaceIncrementalEstimationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented

};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSamplingSchemeQSpaceIncrementalEstimationFilter_hxx)
#include "itkSamplingSchemeQSpaceIncrementalEstimationFilter.hxx"
#endif

#endif 

