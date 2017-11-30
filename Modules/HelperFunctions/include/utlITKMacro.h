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


enum{
  /** itk::Image<double, 4>, stored volume by volume  */
  IMAGE_ND=0,
  /** itk::VectorImage<double, 3>, stored continuously in each voxel, default image type used in dmritool  */
  IMAGE_VECTOR,
  /** itk::Image<itk::VariableLengthVector<double>, double>, a variable length vector in each voxel  */
  IMAGE_VARIABLELENGTH,
  /** itk::SpatiallyDenseSparseVectorImage<double, 3>, a sparse vector in each voxel  */
  IMAGE_SPARSE
};

#ifndef ITK_OVERRIDE

// Copied from itkMacro.h. It is used for old version of ITK without definition of the macros.
#if __cplusplus >= 201103L
// In c++11 the override keyword allows you to explicity define that a function
// is intended to override the base-class version.  This makes the code more
// managable and fixes a set of common hard-to-find bugs.
#define ITK_OVERRIDE override
// In functions that should not be implemented, use the C++11 mechanism
// to ensure that thye are purposely not implemented
#define ITK_DELETE_FUNCTION =delete
// In c++11 there is an explicit nullptr type that introduces a new keyword to
// serve as a distinguished null pointer constant: nullptr. It is of type
// nullptr_t, which is implicitly convertible and comparable to any pointer type
// or pointer-to-member type. It is not implicitly convertible or comparable to
// integral types, except for bool.
#define ITK_NULLPTR  nullptr
// In C++11 the throw-list specification has been deprecated,
// replaces with the noexcept specifier. Using this function
// specification adds the run-time check that the method does not
// throw, if it does throw then std::terminate will be called.
// Use cautiously.
#define ITK_NOEXCEPT noexcept
#define ITK_HAS_CXX11_STATIC_ASSERT
#define ITK_HAS_CXX11_RVREF
#else
#define ITK_OVERRIDE
#define ITK_DELETE_FUNCTION
#define ITK_NULLPTR  NULL
#define ITK_NOEXCEPT throw()
#endif

#endif


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


#define itkTypedefMaskedImageToImageMacro(Superclass)                                    \
  typedef typename Superclass::InputImageType           InputImageType;                  \
  typedef typename Superclass::InputImagePointer        InputImagePointer;               \
  typedef typename Superclass::InputImageConstPointer   InputImageConstPointer;          \
  typedef typename Superclass::InputImageIndexType      InputImageIndexType;             \
  typedef typename Superclass::InputImageSizeType       InputImageSizeType;              \
  typedef typename Superclass::InputImageSpacingType    InputImageSpacingType;           \
  typedef typename Superclass::InputImagePixelType      InputImagePixelType;             \
  typedef typename Superclass::InputImageRegionType     InputImageRegionType;            \
                                                                                         \
  typedef typename Superclass::OutputImageType          OutputImageType;                 \
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;              \
  typedef typename Superclass::OutputImageIndexType     OutputImageIndexType;            \
  typedef typename Superclass::OutputImageSizeType      OutputImageSizeType;             \
  typedef typename Superclass::OutputImageSpacingType   OutputImageSpacingType;          \
  typedef typename Superclass::OutputImagePixelType     OutputImagePixelType;            \
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;           \
                                                                                         \
  typedef typename Superclass::MaskImageType            MaskImageType;                   \
  typedef typename MaskImageType::Pointer               MaskImagePointer;                \
                                                                                         \
  typedef typename itk::Image<double,3>                 ScalarImageType;                 \
  typedef typename ScalarImageType::Pointer             ScalarImagePointer;              



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

