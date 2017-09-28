/**
 *       @file  itkTensorBasisMatrixGenerator.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-18-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkTensorBasisMatrixGenerator_hxx
#define __itkTensorBasisMatrixGenerator_hxx

#include "itkTensorBasisMatrixGenerator.h"
#include "utlRotationMatrixFromVectors.h"

namespace itk
{

template<typename TElement>
TensorBasisMatrixGenerator<TElement>
::TensorBasisMatrixGenerator() 
{
  m_EigenValue1=-1;
  m_EigenValue2=-1;
  m_EigenValue3=-1;
  m_EigenValueISO=-1;
}

template<typename TElement>
typename LightObject::Pointer
TensorBasisMatrixGenerator<TElement>
::InternalClone() const
{
  utlShowPosition(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_EigenValue1 = m_EigenValue1;
  rval->m_EigenValue2 = m_EigenValue2;
  rval->m_EigenValue3 = m_EigenValue3;
  rval->m_EigenValueISO = m_EigenValueISO;
  return loPtr;
}

template<typename TElement>
void
TensorBasisMatrixGenerator<TElement>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar3(true, m_EigenValue1, m_EigenValue2, m_EigenValue3, os<<indent);
  PrintVar1(true, m_EigenValueISO, os<<indent);
}

template<typename TElement>
void
TensorBasisMatrixGenerator<TElement>
::VerifyInputParameters() const
{
  Superclass::VerifyInputParameters();

  utlSAGlobalException(m_EigenValue1<0 || m_EigenValue2<0 || m_EigenValue3<0)(m_EigenValue1)(m_EigenValue2)(m_EigenValue3).msg("need to set m_EigenValue1, m_EigenValue2, m_EigenValue3");
  utlSAGlobalException(this->m_UseIsotropicTerm && m_EigenValueISO<0)(m_EigenValueISO).msg("need to set m_EigenValueISO");

  utlGlobalException(this->m_OutputType==Superclass::DWI && this->m_SamplingSchemeQSpace->GetBVector()->size()==0, "no bVector");
}


template<typename TElement>
void
TensorBasisMatrixGenerator<TElement>
::ComputeBasisMatrix()
{
  utlShowPosition(this->GetDebug());
  VerifyInputParameters();
  
  MatrixPointer qOrientations = this->m_SamplingSchemeQSpace->GetOrientationsCartesian();
  STDVectorPointer bVector = this->m_SamplingSchemeQSpace->GetBVector();
  MatrixPointer rOrientations = this->m_SamplingSchemeRSpace->GetOrientationsCartesian();
  STDVectorPointer rVector = this->m_SamplingSchemeRSpace->GetRadiusVector();
  double tau = this->m_SamplingSchemeQSpace->GetTau();

  int numberOfBasis = this->GetNumberOfBasis();
  int numberOfSamples = -1;
  // if (this->m_OutputType==Superclass::DWI)
  //   numberOfSamples =  this->GetNumberOfQSamples();
  // else if (this->m_OutputType==Superclass::EAP || this->m_OutputType==Superclass::ODF)
  //   numberOfSamples = this->GetNumberOfRSamples();
  numberOfSamples =  this->GetNumberOfSamples();
  MatrixPointer mat  (new MatrixType(numberOfSamples, numberOfBasis) );
  
  TensorType tensor;
  std::vector<double> vec(3,1);
  vec[0]=m_EigenValue1, vec[1]=m_EigenValue2, vec[2]=m_EigenValue3;

  Vector<double, 3> principalDirection(0.0);
  Vector<double, 3> e1;
  e1[0]=1.0, e1[1]=0.0, e1[2]=0.0;
  Matrix<double,3,3> rotation;

  VectorType out(numberOfSamples);
  for ( int i = 0; i < this->m_BasisOrientations->Rows(); i += 1 ) 
    {
    tensor.Fill(0.0);
    tensor.SetEigenValues(vec);

    principalDirection[0] = (*this->m_BasisOrientations)(i,0);
    principalDirection[1] = (*this->m_BasisOrientations)(i,1);
    principalDirection[2] = (*this->m_BasisOrientations)(i,2);
    utl::RotationMatrixFromVectors<Vector<double,3>, Matrix<double,3> >(e1, principalDirection, rotation); 
    tensor.Rotate(rotation);

    // std::cout << "tensor = " << tensor << std::endl << std::flush;
    if (this->m_OutputType==Superclass::DWI)
      tensor.GetDWISamples(out, *qOrientations, *bVector);
    else if (this->m_OutputType==Superclass::EAP)
      tensor.GetEAPSamples(out, *rOrientations, *rVector, tau);
    else if (this->m_OutputType==Superclass::ODF)
      tensor.GetODFSamples(out, *rOrientations, this->m_ODFOrder, false);

    // utl::PrintVnlVector(out, "out");
    if (this->m_UseIsotropicTerm)
      mat->SetColumn(i+1, out);
    else
      mat->SetColumn(i, out);
    }

  if (this->m_UseIsotropicTerm)
    {
    vec[0]=m_EigenValueISO, vec[1]=m_EigenValueISO, vec[2]=m_EigenValueISO;
    tensor.Fill(0.0);
    tensor.SetEigenValues(vec);
    
    if (this->m_OutputType==Superclass::DWI)
      tensor.GetDWISamples(out, *qOrientations, *bVector);
    else if (this->m_OutputType==Superclass::EAP)
      tensor.GetEAPSamples(out, *rOrientations, *rVector, tau);
    else if (this->m_OutputType==Superclass::ODF)
      tensor.GetODFSamples(out, *rOrientations, this->m_ODFOrder, false);
    
    mat->SetColumn(0, out);
    }

  if (this->m_OutputType==Superclass::DWI)
    this->m_QBasisMatrixForDWI = mat;
  else if (this->m_OutputType==Superclass::EAP)
    this->m_RBasisMatrixForEAP = mat;
  else if (this->m_OutputType==Superclass::ODF)
    this->m_RBasisMatrixForODF = mat;

}




}


#endif 




