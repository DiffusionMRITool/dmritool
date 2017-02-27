/**
 *       @file  itkSphericalPolarFourierGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-09-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkSphericalPolarFourierGenerator_h
#define __itkSphericalPolarFourierGenerator_h

#include <itkObject.h>
#include <itkObjectFactory.h>
#include "itkDiffusionTensor.h"


namespace itk 
{

/**
 *   \class   SphericalPolarFourierRadialGenerator
 *   \brief   radial part of general SPF basis
 *   \ingroup Math
 *   \author  Jian Cheng  
 */
template<class PreciseType = double> 
class ITK_EXPORT SphericalPolarFourierRadialGenerator : public Object
{
public:
  /** Standard class typedefs. */
  typedef SphericalPolarFourierRadialGenerator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SphericalPolarFourierRadialGenerator, Object);

  typedef enum {
  SPF=0, //=0
  DSPF,  // =1
  SHORE // =2
  } SPFType;



  itkSetMacro(SPFType, SPFType);
  itkGetMacro(SPFType, SPFType);
  
  itkSetMacro(Scale, double );
  itkGetMacro(Scale, double );
  
  itkSetMacro(N, int);
  itkGetMacro(N, int);
  
  itkSetMacro(L, int);
  itkGetMacro(L, int);
  
  void SetNLSPF(const int n, const int l, const SPFType spf)
    {
    SetN(n);  SetL(l);  SetSPFType(spf);
    }

  PreciseType Evaluate(const PreciseType qrValue, const bool isFourier=false) const;

  PreciseType GetNormalizeFacotr(const bool isFourier=false) const;

protected:
  /** SphericalPolarFourierRadialGenerator constructor  */
  SphericalPolarFourierRadialGenerator ();

  virtual ~SphericalPolarFourierRadialGenerator ()
    {
    }

  int m_N;
  int m_L;
  SPFType m_SPFType;
  double m_Scale;

private :
  SphericalPolarFourierRadialGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

/**
 *   \class   SphericalPolarFourierGenerator
 *   \brief   general SPF basis in radial part, SH basis in spherical part
 *   \ingroup Math
 *   \author  Jian Cheng  
 */
template<class PreciseType = double> 
class ITK_EXPORT SphericalPolarFourierGenerator : public Object
{
public:
  /** Standard class typedefs. */
  typedef SphericalPolarFourierGenerator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SphericalPolarFourierGenerator, Object);

  typedef SphericalPolarFourierRadialGenerator<PreciseType>         RadialType;
  typedef typename RadialType::Pointer         RadialPointer;
  typedef typename RadialType::SPFType         SPFType;
  using RadialType::SPF;
  using RadialType::DSPF;
  using RadialType::SHORE;
  

  void SetN(const int n) 
    { 
    if (n!=m_Radial->GetN()) 
      {
      m_Radial->SetN(n); 
      this->Modified();
      }
    }
  void SetL(const int l)
    { 
    if (l!=m_Radial->GetL()) 
      {
      m_Radial->SetL(l); 
      this->Modified();
      }
    }

  int GetN() const {return m_Radial->GetN();}
  int GetL() const {return m_Radial->GetL();}
  
  itkSetMacro(M, int);
  itkGetMacro(M, int);

  void SetNLM(const int n, const int l, const int m)
    {
    SetN(n);
    SetL(l);
    SetM(m);
    }

  void SetScale(const double scale) 
    { 
    if (scale!=m_Radial->GetScale()) 
      {
      m_Radial->SetScale(scale); 
      this->Modified();
      }
    }
  double GetScale() const {return m_Radial->GetScale();}
  
  void SetSPFType(const SPFType model) 
    { 
    if (model!=m_Radial->GetSPFType()) 
      {
      m_Radial->SetSPFType(model); 
      this->Modified();
      }
    }
  typename RadialType::SPFType GetSPFType() const {return m_Radial->GetSPFType();}

  PreciseType EvaluateRadial(const PreciseType qrValue, const bool isFourier=false) const;

  PreciseType Evaluate(const PreciseType qrValue, const PreciseType theta, const PreciseType phi, const bool isFourier=false) const;
  
  PreciseType GetNormalizeFacotr(const bool isFourier=false) const
    {
    return m_Radial->GetNormalizeFacotr(isFourier);
    }

protected:
  /** SphericalPolarFourierGenerator constructor  */
  SphericalPolarFourierGenerator ()
    {
    m_Radial = RadialType::New();
    m_M = 0;
    }

  virtual ~SphericalPolarFourierGenerator ()
    {
    }
  
  void PrintSelf(std::ostream & os, Indent indent) const;

  RadialPointer m_Radial;
  int m_M;

private :
  SphericalPolarFourierGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};


}

// Define instantiation macro for this template.
#define ITK_TEMPLATE_SphericalPolarFourierRadialGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT SphericalPolarFourierRadialGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef SphericalPolarFourierRadialGenerator< ITK_TEMPLATE_2 x > SphericalPolarFourierRadialGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphericalPolarFourierRadialGenerator+-.h"
#endif

#define ITK_TEMPLATE_SphericalPolarFourierGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT SphericalPolarFourierGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef SphericalPolarFourierGenerator< ITK_TEMPLATE_2 x > SphericalPolarFourierGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphericalPolarFourierGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSphericalPolarFourierGenerator_hxx)
# include "itkSphericalPolarFourierGenerator.hxx"
#endif

#endif 
