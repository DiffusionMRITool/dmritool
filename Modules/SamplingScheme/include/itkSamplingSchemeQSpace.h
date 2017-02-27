/**
 *       @file  itkSamplingSchemeQSpace.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-25-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpace_h
#define __itkSamplingSchemeQSpace_h

#include "itkSamplingScheme3D.h"
#include "itkObject.h"

namespace itk 
{

/**
 *   \class   SamplingSchemeQSpace
 *   \brief   this class describes sampling in Q space. 
 *
 *   The sampling in Q space can be single shell sampling or multiple shell sampling, 
 *   which is determined by m_IndicesInShells. 
 *
 *   \author  Jian Cheng  
 *   \ingroup SamplingScheme
 */
template<class TPixelType = double> 
class ITK_EXPORT SamplingSchemeQSpace : public SamplingScheme3D<TPixelType>
{
public:
  /** Standard class typedefs. */
  typedef SamplingSchemeQSpace Self;
  typedef SamplingScheme3D<TPixelType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SamplingSchemeQSpace, SamplingScheme3D);
  
  typedef typename Superclass::STDVectorType          STDVectorType;
  typedef typename Superclass::STDVectorPointer       STDVectorPointer;
  typedef typename Superclass::IndexVectorType        IndexVectorType;
  typedef typename Superclass::Index2DVectorType      Index2DVectorType;
  typedef typename Superclass::Index2DVectorPointer   Index2DVectorPointer;
  
  itkSetMacro(BThresholdSingleShell, double);
  itkGetMacro(BThresholdSingleShell, double);

  void SetBVector(const STDVectorPointer bVec);
  itkGetMacro(BVector, STDVectorPointer);
  STDVectorPointer GetBVectorInShell(unsigned int shellIndex);

  void SetSamplingScheme3D(typename Superclass::Pointer scheme3D);
  
  void Clear();

  void ConvertBVectorToQVector();
  
  void ConvertQVectorToBVector();
  
  std::vector<STDVectorType> GroupBValues();
  
  void CorrectBValues();
  void CorrectRadiusValues();
  
  /** remove samples not in m_IndicesInShells  */
  void RemoveSamplesNotIndexed();

protected:
  SamplingSchemeQSpace ();

  virtual ~SamplingSchemeQSpace ()
    {
    }
  
  typename LightObject::Pointer InternalClone() const;
 
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** b values whose distance is smallter the threshold will be considered in 
   * the same shell, and they will be replaced as their mean b value. */
  double m_BThresholdSingleShell;

  STDVectorPointer m_BVector;

private:
  SamplingSchemeQSpace(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

  // std::vector<STDVectorType> GroupRadiusValues();
  // void CorrectRadiusValues();
};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSamplingSchemeQSpace_hxx)
#include "itkSamplingSchemeQSpace.hxx"
#endif

#endif 

