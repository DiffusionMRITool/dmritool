/**
 *       @file  itkSamplingSchemeQSpaceWriter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-15-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSamplingSchemeQSpaceWriter_h
#define __itkSamplingSchemeQSpaceWriter_h

#include "itkLightProcessObject.h"


namespace itk
{
/**
 *   \class   SamplingSchemeQSpaceWriter
 *   \brief   write orientations and b values (single/multiple shells) into text files
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup SamplingScheme
 */
template <class TSamplingType>
class ITK_EXPORT SamplingSchemeQSpaceWriter : public LightProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SamplingSchemeQSpaceWriter        Self;
  typedef LightProcessObject          Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer< const Self >  ConstPointer;
  
  typedef TSamplingType                        SamplingType;
  typedef typename SamplingType::Pointer       SamplingPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SamplingSchemeQSpaceWriter, LightProcessObject);
  
  itkSetObjectMacro(Sampling, SamplingType);
  itkGetObjectMacro(Sampling, SamplingType);

  itkSetMacro(BFile, std::string);
  itkGetMacro(BFile, std::string);
  
  itkSetMacro(OrientationFile, std::string);
  itkGetMacro(OrientationFile, std::string);
  
  itkSetMacro(SaveSingleShell, bool);
  itkGetMacro(SaveSingleShell, bool);
  itkBooleanMacro(SaveSingleShell);
  
  itkSetMacro(SaveAllShellsInOneFile, bool);
  itkGetMacro(SaveAllShellsInOneFile, bool);
  itkBooleanMacro(SaveAllShellsInOneFile);
  
  
  void Update();

private:
  SamplingSchemeQSpaceWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

protected:

  SamplingSchemeQSpaceWriter();
  
  ~SamplingSchemeQSpaceWriter() {}

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  /** If it is set, save b values  */
  std::string m_BFile;

  /** This is must be set  */
  std::string m_OrientationFile;
  
  // [>* If it is set, save the index for single shell data.  <]
  // std::string m_IndexFile;

  /** input sampling scheme including orientations and b values  */
  SamplingPointer m_Sampling;

  /** If it is true, save single shell data. */
  bool m_SaveSingleShell;

  /** If it is set, save all data.  */
  bool m_SaveAllShellsInOneFile;
  
};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSamplingSchemeQSpaceWriter_hxx)
#include "itkSamplingSchemeQSpaceWriter.hxx"
#endif

#endif 

