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

#define itkSetGetBooleanMacro(name)  \
  itkSetMacro(name,bool);            \
  itkGetMacro(name,bool);            \
  itkBooleanMacro(name);

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

#endif 

