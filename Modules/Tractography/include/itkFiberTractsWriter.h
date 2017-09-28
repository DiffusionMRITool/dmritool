/**
 *       @file  itkFiberTractsWriter.h
 *      @brief  
 *     Created  "07-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkFiberTractsWriter_h
#define __itkFiberTractsWriter_h

#include <itkLightProcessObject.h>

#include "itkFiberTracts.h"
#include "utlCoreMacro.h"

namespace itk
{

class ITK_EXPORT FiberTractsWriter : public LightProcessObject
{
public:
  /** Standard class typedefs. */
  typedef FiberTractsWriter       Self;
  typedef LightProcessObject   Superclass;
  typedef SmartPointer<Self>   Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FiberTractsWriter, LightProcessObject);

  typedef FiberTracts<double>                     FiberTractsType;
  typedef typename FiberTractsType::Pointer       FiberTractsPointer;
  typedef typename FiberTractsType::FiberType     FiberType;
  typedef typename FiberTractsType::FiberPointer  FiberPointer;
  typedef typename FiberType::STDVectorType       STDVectorType;
  typedef typename FiberType::VertexType          VertexType;

  itkSetGetMacro(FileName, std::string);
  itkSetGetMacro(FiberTracts, FiberTractsPointer);

  /** Does the real work. */
  virtual void Update();

  void WriteTractsTRK();

  void WriteTractsTCK();
  
  // void WriteTractsVTK();
    
protected:
  FiberTractsWriter(){}
  ~FiberTractsWriter(){}

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    PrintVar(true, os<< indent, m_FileName);
    m_FiberTracts->Print(os, indent);
    }
  
  std::string m_FileName;

  FiberTractsPointer m_FiberTracts = FiberTractsType::New();


private:
  FiberTractsWriter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberTractsWriter.hxx"
#endif


#endif 
