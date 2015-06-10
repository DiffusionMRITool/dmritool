/**
 *       @file  itkCylinderModelGenerator.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "06-02-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkCylinderModelGenerator_h
#define __itkCylinderModelGenerator_h

#include "itkObject.h"
#include "itkDiffusionModelGenerator.h"
#include "itkFunctors.h"
#include "itkUnaryFunctorLookUpTable.h"

namespace itk
{

/**
 *   \class   CylinderModelGenerator
 *   \brief   Cylinder Model
 *
 *   Reference: 
 *   Resolution of complex tissue microarchitecture using the diffusion orientation transform (DOT), Evren Ozarslan, NeuroImage 2006
 *    
 * \ingroup DiffusionModels
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template<class PreciseType = double> 
class ITK_EXPORT CylinderModelGenerator : public DiffusionModelGenerator<PreciseType>
{
public:
  /** Standard class typedefs. */
  typedef CylinderModelGenerator Self;
  typedef DiffusionModelGenerator<PreciseType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(CylinderModelGenerator, DiffusionModelGenerator);
  
  typedef typename Superclass::MatrixType       MatrixType;
  typedef typename Superclass::VectorType       VectorType;
  typedef typename Superclass::MatrixPointer      MatrixPointer;
  typedef typename Superclass::VectorPointer      VectorPointer;
  typedef typename Superclass::STDVectorType    STDVectorType;
  typedef typename Superclass::STDVectorPointer   STDVectorPointer;
  
  typedef typename Superclass::SamplingSchemeQSpaceType        SamplingSchemeQSpaceType;
  typedef typename Superclass::SamplingSchemeQSpacePointer     SamplingSchemeQSpacePointer;
  typedef typename Superclass::SamplingSchemeRSpaceType        SamplingSchemeRSpaceType;
  typedef typename Superclass::SamplingSchemeRSpacePointer     SamplingSchemeRSpacePointer;
  typedef typename Superclass::PointType                       PointType;

  typedef Image<double>  Image3DType;
  typedef typename Image3DType::Pointer  Image3DPointer;
  
  typedef UnaryFunctorLookUpTable<Functor::EXP<double> >                   LUTExpType;
  typedef typename LUTExpType::Pointer                         LUTExpPointer;
  
  itkSetMacro(CylinderAxis, PointType);
  itkGetMacro(CylinderAxis, PointType);

  itkSetMacro(Length, double);
  itkGetMacro(Length, double);
  
  itkSetMacro(Radius, double);
  itkGetMacro(Radius, double);

  itkSetMacro(D0, double);
  itkGetMacro(D0, double);
  
  void Rotate (const MatrixType& mat);

  void ComputeDWISamples ();

  void ComputeEAPSamples ();
  
  void ComputeODFSamples ();
  
  void BuildTable ();
  
  // void VerifyInputParameters() const;

protected:
  /** CylinderModelGenerator constructor  */
  CylinderModelGenerator (); 

  virtual ~CylinderModelGenerator ()
    {
    }

  void PrintSelf(std::ostream& os, Indent indent) const;

  typename LightObject::Pointer InternalClone() const;

  double m_Length;
  double m_Radius;
  double m_D0;

  PointType m_CylinderAxis;

  LUTExpPointer m_LUTExp;

  Image3DPointer m_EAPVolumeForZAxis;
  Image3DPointer m_ODFVolumeForZAxis;

private :
  CylinderModelGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkCylinderModelGenerator_hxx)
#include "itkCylinderModelGenerator.hxx"
#endif

#endif 
