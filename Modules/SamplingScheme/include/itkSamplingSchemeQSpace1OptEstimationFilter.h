/**
 *       @file  itkSamplingSchemeQSpace1OptEstimationFilter.h
 *      @brief  
 *     Created  "01-16-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#ifndef __itkSamplingSchemeQSpace1OptEstimationFilter_h
#define __itkSamplingSchemeQSpace1OptEstimationFilter_h

#include "itkSamplingSchemeQSpaceEstimationFilter.h"
#include "itkSamplingScheme3D.h"

namespace itk
{

/**
 *   \class   SamplingSchemeQSpace1OptEstimationFilter
 *   \brief   incremental estimation of single/multi-shell orientations
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup SamplingScheme
 */
template <class TSamplingType>
class ITK_EXPORT SamplingSchemeQSpace1OptEstimationFilter 
  : public SamplingSchemeQSpaceEstimationFilter<TSamplingType>
{

public:
  /** Standard class typedefs. */
  typedef SamplingSchemeQSpace1OptEstimationFilter Self;
  typedef SamplingSchemeQSpaceEstimationFilter<TSamplingType>                 Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;
  
  typedef TSamplingType                                       SamplingType;
  typedef typename SamplingType::Pointer                      SamplingPointer;
  typedef typename TSamplingType::PointType                   PointType;
  
  typedef typename Superclass::ValueType                    ValueType;
  typedef typename Superclass::MatrixType                   MatrixType;
  typedef typename Superclass::MatrixPointer                MatrixPointer;
  typedef typename Superclass::STDVectorType                STDVectorType;
  typedef typename Superclass::STDVectorPointer             STDVectorPointer;
  typedef typename Superclass::IndexVectorType              IndexVectorType;
  typedef typename Superclass::Index2DVectorType            Index2DVectorType;
  typedef typename Superclass::Index2DVectorPointer         Index2DVectorPointer;

  /** Standard Macros */
  itkTypeMacro(SamplingSchemeQSpace1OptEstimationFilter, SamplingSchemeQSpaceEstimationFilter);
  itkNewMacro(Self);
  
  itkSetGetMacro(TessellationOrder, unsigned int);
  
  itkSetGetMacro(FineOrientations, MatrixPointer);

  void GenerateData() ITK_OVERRIDE;

protected:
  SamplingSchemeQSpace1OptEstimationFilter();
  ~SamplingSchemeQSpace1OptEstimationFilter(){}
  
  void Initialization() ITK_OVERRIDE;
  
  /** the order of tessellation for the orignal fine mesh  */
  unsigned int m_TessellationOrder;
  
  /** It is the input data or generated from m_TessellationOrder  */
  MatrixPointer    m_FineOrientations;
  
  SamplingPointer m_FineScheme;

private:
  SamplingSchemeQSpace1OptEstimationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented

};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSamplingSchemeQSpace1OptEstimationFilter_hxx)
#include "itkSamplingSchemeQSpace1OptEstimationFilter.hxx"
#endif

#endif 

