/**
 *       @file  itkFiberTractsReader.h
 *      @brief  
 *     Created  "08-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkFiberTractsReader_h
#define __itkFiberTractsReader_h

#include <itkLightProcessObject.h>

#include "itkFiberTracts.h"
#include "utlCoreMacro.h"

namespace itk
{

class ITK_EXPORT FiberTractsReader : public LightProcessObject
{
public:
  /** Standard class typedefs. */
  typedef FiberTractsReader       Self;
  typedef LightProcessObject   Superclass;
  typedef SmartPointer<Self>   Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FiberTractsReader, LightProcessObject);

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

  void ReadTractsTRK();

  void ReadTractsTCK();
  
  void ReadTractsVTK();
    
protected:
  FiberTractsReader(){}
  ~FiberTractsReader(){}

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    PrintVar(true, os<< indent, m_FileName);
    m_FiberTracts->Print(os, indent);
    }
  
  std::string m_FileName;

  FiberTractsPointer m_FiberTracts = FiberTractsType::New();


private:
  FiberTractsReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberTractsReader.hxx"
#endif


#endif 
