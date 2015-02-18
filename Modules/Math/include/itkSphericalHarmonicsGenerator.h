/**
 *       @file  itkSphericalHarmonicsGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  10-29-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */  

#ifndef __itkSphericalHarmonicsGenerator_h
#define __itkSphericalHarmonicsGenerator_h

#include <complex>
#include <itkObject.h>
#include <itkObjectFactory.h>


namespace itk 
{

/** \class SphericalHarmonicsGenerator
 *  \brief Generate complex and real Spherical Harmonic Basis.
 *
 *  The complex Spherical Harmonic basis is defined as
 *  \f[ 
 *    y_l^m = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}} e^{im\phi} P_l^m(\cos\theta) 
 *  \f]
 *  The real SH basis is defined as
 *  \f[ 
 *  Y_l^m(\theta,\phi) = 
 *  \left\{  \begin{array}{lcl}
 *  \sqrt{2}\mbox{Re}(y_l^{|m|}(\theta,\phi))  = \sqrt{2}(-1)^m\mbox{Re}(y_l^m(\theta,\phi)) &  \mbox{if} &  -l\leq m<0  \\  
 *  y_l^m(\theta,\phi)    &  \mbox{if}  &  m=0  \\
 *  \sqrt{2}\mbox{Im}(y_l^m(\theta,\phi))  &  \mbox{if}  &  l\geq m>0 
 *  \end{array}\right. 
 *  \f]
 *  where \f$ y_l^m(\theta,\phi)\f$ is the complex SH basis 
 *
 *  \author Jian Cheng
 *
 */
template<class PreciseType = double> 
class ITK_EXPORT SphericalHarmonicsGenerator : public Object
{
public:
  /** Standard class typedefs. */
  typedef SphericalHarmonicsGenerator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SphericalHarmonicsGenerator, Object);

  /** get value of Complex SH basis  */
  static std::complex<PreciseType> ComplexSH(const int l, const int m, const PreciseType theta, const PreciseType phi);
  /** get value of Real SH basis  */
  static PreciseType RealSH(const int l, const int m, const PreciseType theta, const PreciseType phi);

  /** get value of the triple integral of Complex SH basis  */
  static PreciseType ComplexTripleIntegration(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3);
  /** get value of the triple integral of Real SH basis  */
  static PreciseType RealTripleIntegration(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3, const bool is_precalculated=true);

  /** get value of the derivative of theta for Complex SH basis  */
  static std::complex<double> ComplexDerivativeOfTheta ( const int l, const int m, const double theta, const double phi);
  /** get value of the derivative of phi for Complex SH basis  */
  static std::complex<double> ComplexDerivativeOfPhi ( const int l, const int m, const double theta, const double phi);
  
  /** get value of the derivative of theta for Real SH basis  */
  static double RealDerivativeOfTheta ( const int l, const int m, const double theta, const double phi);
  /** get value of the derivative of phi for Real SH basis  */
  static double RealDerivativeOfPhi ( const int l, const int m, const double theta, const double phi);
  
protected:
  SphericalHarmonicsGenerator();
  virtual ~SphericalHarmonicsGenerator();

  void PrintSelf(std::ostream & os, Indent indent) const;

private :
  SphericalHarmonicsGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


}

  
// Define instantiation macro for this template.
#define ITK_TEMPLATE_SphericalHarmonicsGenerator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT SphericalHarmonicsGenerator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef SphericalHarmonicsGenerator< ITK_TEMPLATE_2 x > SphericalHarmonicsGenerator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphericalHarmonicsGenerator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSphericalHarmonicsGenerator_hxx)
# include "itkSphericalHarmonicsGenerator.hxx"
#endif

#endif
