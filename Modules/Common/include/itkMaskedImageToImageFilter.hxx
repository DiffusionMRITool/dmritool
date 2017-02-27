/**
 *       @file  itkMaskedImageToImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-20-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkMaskedImageToImageFilter_hxx
#define __itkMaskedImageToImageFilter_hxx

#include "itkMaskedImageToImageFilter.h"
#include "itkStdStreamLogOutput.h"
#include "itkProgressReporter.h"

#include "utl.h"

namespace itk
{
template< class TInputImage, class TOutputImage, class TMaskImage >
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::MaskedImageToImageFilter() : Superclass() 
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
  m_MaskImage = NULL;
  m_ThreadID=-1;

  m_LogLevel= LOG_NORMAL;

  // NOTE: m_Logger makes the estimation a little bit slower once it is newed, even though it is not used. 
  // m_Logger = itk::ThreadLogger::New();
  // m_LoggerVector= LoggerVectorPointer(new LoggerVectorType())
}

template< class TInputImage, class TOutputImage, class TMaskImage >
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::~MaskedImageToImageFilter()
{}

template< class TInputImage, class TOutputImage, class TMaskImage >
void
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::SetMaskImage(const std::string file)
{
  MaskImagePointer mask = MaskImageType::New();
  itk::ReadImage<MaskImageType>(file, mask);
  this->SetMaskImage(mask);
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::VerifyMaskInformation() const
{
  utlGlobalException(!m_MaskImage, "no m_MaskImage");
  InputImagePointer input = const_cast<InputImageType *>(this->GetInput());
  utlGlobalException(!input, "no input");
  utlGlobalException((!itk::VerifyImageInformation<InputImageType, MaskImageType>(input, m_MaskImage, true)), "wrong mask information");
  // utlGlobalException((!itk::VerifyImageSize<InputImageType, MaskImageType>(input, m_MaskImage, true)), "wrong mask information");
}

template< class TInputImage, class TOutputImage, class TMaskImage >
typename LightObject::Pointer
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_MaskImage = m_MaskImage;
  rval->m_Logger = m_Logger;
  rval->m_LoggerVector = m_LoggerVector;
  rval->m_ThreadID = m_ThreadID;
  rval->m_LogLevel = m_LogLevel;
  rval->SetNumberOfThreads(this->GetNumberOfThreads());
  rval->SetDebug(this->GetDebug());
  return loPtr;
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::CreateLoggerVector()
{
  utlException(!(this->GetDebug() && this->GetNumberOfThreads()>1), "create m_LoggerVector only when m_Debug is true and multiple threads are used");
  if (!m_Logger)
    {
    itk::StdStreamLogOutput::Pointer coutput = itk::StdStreamLogOutput::New();
    coutput->SetStream(std::cout);
    m_Logger = itk::ThreadLogger::New();
    m_Logger->SetPriorityLevel(itk::LoggerBase::DEBUG);
    m_Logger->SetLevelForFlushing(itk::LoggerBase::MUSTFLUSH);
    m_Logger->AddLogOutput(coutput);
    }

  int numberofThreads = this->GetNumberOfThreads();
  m_LoggerVector = LoggerVectorPointer(new LoggerVectorType(numberofThreads));
  for (int ii = 0; ii < numberofThreads; ++ii)
    (*m_LoggerVector)[ii] = m_Logger;
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::WriteLogger(const std::string str, const LoggerBase::PriorityLevelType level) const
{
  if (m_ThreadID>=0 && this->GetDebug() && this->GetNumberOfThreads()>1)
    {
    utlSAException(this->m_LoggerVector->size()!=this->GetNumberOfThreads())
      (this->m_LoggerVector->size())(this->GetNumberOfThreads()).msg("need to set m_LoggerVector");
    (*this->m_LoggerVector)[m_ThreadID]->Write(level, str);
    (*this->m_LoggerVector)[m_ThreadID]->Flush();
    }
  else
    std::cout << str << std::flush;
}

template< class TInputImage, class TOutputImage, class TMaskImage >
std::string
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::ThreadIDToString() const
{
  if (m_ThreadID>=0 && this->GetDebug() && this->GetNumberOfThreads()>1)
    {
    std::ostringstream msg;
    msg << "<Thread " << m_ThreadID << "> : ";
    return msg.str();
    }
  else
    return "";
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::InitializeThreadedLibraries()
{
  utl::InitializeThreadedLibraries(this->GetNumberOfThreads());
}

template< class TInputImage, class TOutputImage, class TMaskImage >
void
MaskedImageToImageFilter< TInputImage, TOutputImage, TMaskImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  if (!IsImageEmpty(m_MaskImage))
    os << indent << "MaskImage = " << m_MaskImage << std::endl << std::flush;
}

}

#endif 


