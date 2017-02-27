/**
 *       @file  utlITKMacro.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "12-31-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlITKMacro_h
#define __utlITKMacro_h

#include "itkMacro.h"


/** @addtogroup utlHelperFunctions
@{ */


#define itkSetNoConstMacro(name, type)                      \
  virtual void Set##name (type _arg)                        \
    {                                                       \
    itkDebugMacro("setting " #name " to " << _arg);         \
    if ( this->m_##name != _arg )                           \
      {                                                     \
      this->m_##name = _arg;                                \
      this->Modified();                                     \
      }                                                     \
    }


#define itkGetPointerMacro(name, type)  \
  virtual type * Get##name ()           \
    {                                   \
    return &(m_##name);                 \
    }


#define itkGetConstPointerMacro(name, type)   \
  virtual const type * Get##name () const     \
    {                                         \
    return &(m_##name);                       \
    }


#define itkSetNDebugMacro(name, type)                      \
  virtual void Set##name (const type _arg)                 \
    {                                                      \
    if ( this->m_##name != _arg )                          \
      {                                                    \
      this->m_##name = _arg;                               \
      this->Modified();                                    \
      }                                                    \
    }


#define itkSetNDebugByConstReferenceMacro(name, type)      \
  virtual void Set##name (const type& _arg)                \
    {                                                      \
    if ( this->m_##name != _arg )                          \
      {                                                    \
      this->m_##name = _arg;                               \
      this->Modified();                                    \
      }                                                    \
    }


#define itkSetByConstReferenceMacro(name, type)            \
  virtual void Set##name (const type& _arg)                \
    {                                                      \
    itkDebugMacro("setting " #name " to " << _arg);        \
    if ( this->m_##name != _arg )                          \
      {                                                    \
      this->m_##name = _arg;                               \
      this->Modified();                                    \
      }                                                    \
    }


#define itkSetGetMacro(name, type)   \
  itkSetMacro(name, type);           \
  itkGetMacro(name, type);


#define itkSetGetBooleanMacro(name)  \
  itkSetMacro(name,bool);            \
  itkGetMacro(name,bool);            \
  itkBooleanMacro(name);


#define itkFunctorSetMacro(name, type)      \
void Set##name(type name)                   \
  {                                         \
  this->GetFunctor().Set##name(name);       \
  this->Modified();                         \
  }                                         


#define itkFunctorGetMacro(name, type)      \
type Get##name() const                      \
{                                           \
  return this->GetFunctor().Get##name();    \
}


#define itkFunctorSetGetMacro(name, type)   \
  itkFunctorSetMacro(name, type);           \
  itkFunctorGetMacro(name, type);


/** Show Position when debug and muti-thread are used. It is only used for filters derived from MaskedImageToImageFilter  */
#define itkShowPositionThreadedLogger(cond)                                                                            \
do                                                                                                                     \
{                                                                                                                      \
  if ((cond))                                                                                                          \
    {                                                                                                                  \
    if (this->m_ThreadID>=0 && this->GetDebug() && this->GetNumberOfThreads()>1)                                       \
      {                                                                                                                \
      std::string __threadIDStr = this->ThreadIDToString();                                                            \
      std::ostringstream __msg;                                                                                        \
      __msg << "\n" << __threadIDStr <<  "Work Flow Position: In File: " <<__FILE__<< ", Line: " << __LINE__ << "\n"   \
      << __threadIDStr << "Function: " << __utl_LOCATION__ << "\n" << std::flush;                                      \
      this->WriteLogger(__msg.str());                                                                                  \
      }                                                                                                                \
    else                                                                                                               \
      utlOSShowPosition(cond, std::cout);                                                                              \
    }                                                                                                                  \
} while(0)


#define itkThreadedLogger(cond, msgStr)                                                                   \
do                                                                                                        \
  {                                                                                                       \
  if ((cond))                                                                                             \
    {                                                                                                     \
    std::ostringstream __msg;                                                                             \
    __msg << (msgStr) << std::endl << std::flush;                                                         \
    this->WriteLogger(__msg.str());                                                                       \
    }                                                                                                     \
  } while ( 0 )                                                                                           


/** It should be used with itk::VectorImageChannelFilter to perform scalarImageFilter in each channel of a VectorImage, 
 * then compose it back to output (a VectorImage). */
#define itkVectorImageFilterComposeMacro(VectorImageType, ScalarImageFilterType, scalarImageFilter, input, output)        \
  do {                                                                                                                    \
  typedef typename itk::VectorImageChannelFilter<VectorImageType, VectorImageType, ScalarImageFilterType>  FilterType;    \
  typename FilterType::Pointer filter = FilterType::New();                                                                \
  utlException(std::string(input->GetNameOfClass())!="VectorImage", "input should be VectorImage");                       \
  filter->SetInput(input);                                                                                                \
  filter->SetFilter(scalarImageFilter);                                                                                   \
  filter->Update();                                                                                                       \
  output = filter->GetOutput();                                                                                           \
  } while(0)                                                                                                              \

    /** @} */

#endif 

