/**
 *       @file  itkSamplingSchemeQSpaceIMOCEstimationFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "10-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpaceIMOCEstimationFilter_h
#define __itkSamplingSchemeQSpaceIMOCEstimationFilter_h

#include "itkSamplingSchemeQSpaceEstimationFilter.h"
#include "itkSamplingScheme3D.h"

#include "itkListSample.h"
#include "itkKdTreeGenerator.h"


namespace itk
{

/**
 *   \class   SamplingSchemeQSpaceIMOCEstimationFilter
 *   \brief   Estimation of single/multi-shell orientations using Iterative Maximum Overlap Construction (IMOC)
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup SamplingScheme
 */
template <class TSamplingType>
class ITK_EXPORT SamplingSchemeQSpaceIMOCEstimationFilter 
  : public SamplingSchemeQSpaceEstimationFilter<TSamplingType>
{

public:
  /** Standard class typedefs. */
  typedef SamplingSchemeQSpaceIMOCEstimationFilter Self;
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

  typedef itk::Point<double,3> MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef TreeGeneratorType::KdTreeType TreeType;

  /** Standard Macros */
  itkTypeMacro(SamplingSchemeQSpaceIMOCEstimationFilter, SamplingSchemeQSpaceEstimationFilter);
  itkNewMacro(Self);
  
  itkSetGetMacro(TessellationOrder, unsigned int);
  
  itkSetGetMacro(FineOrientations, MatrixPointer);

  itkSetGetMacro(AngleMinChange, double);

  itkSetGetMacro(ChooseMinimalCoverageShell, bool);

  void GenerateData() ITK_OVERRIDE;

protected:
  SamplingSchemeQSpaceIMOCEstimationFilter();
  ~SamplingSchemeQSpaceIMOCEstimationFilter(){}
  
  void Initialization() ITK_OVERRIDE;

  bool IsSatisfiedSeparationAngles(const std::vector<double>& angles);
  
  /** the order of tessellation for the orignal fine mesh  */
  unsigned int m_TessellationOrder=7;
  
  /** It is generated from m_TessellationOrder  */
  MatrixPointer m_FineOrientations;
  SamplingPointer m_FineScheme = SamplingType::New();
  
  double m_MinDistanceInFineScheme=-1;

  typename TreeGeneratorType::Pointer m_TreeGenerator=nullptr;
  typename TreeType::Pointer m_KDTree= nullptr;
  typename SampleType::Pointer m_Sample;

  double m_AngleMinChange=0.0001;

  /** 
   * If true, always choose the shell with the minimal total coverage in each step.  
   * It is used only for multi-shell case where the number of samples are different in shells. 
   * It is better to always set it as true.
   * */
  bool m_ChooseMinimalCoverageShell=true;

private:
  SamplingSchemeQSpaceIMOCEstimationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented

};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSamplingSchemeQSpaceIMOCEstimationFilter_hxx)
#include "itkSamplingSchemeQSpaceIMOCEstimationFilter.hxx"
#endif


#endif 
