/**
 *       @file  itkSphericalPolarFourierEstimationImageFilter.hxx
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

#ifndef __itkSphericalPolarFourierEstimationImageFilter_hxx
#define __itkSphericalPolarFourierEstimationImageFilter_hxx

#include "itkSphericalPolarFourierEstimationImageFilter.h"
#include "utl.h"

namespace itk
{

template< class TInputImage, class TOutputImage >
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::SphericalPolarFourierEstimationImageFilter() : Superclass(), 
  m_BasisCombinationMatrix(new MatrixType()), 
  m_BasisEnergyDL(new VectorType()), 
  m_BasisMatrixForB0(new MatrixType())
{
  m_BasisScale = -1.0;
  m_EstimationType = LS;
  m_LambdaSpherical = 0;
  m_LambdaRadial = 0;
  m_LambdaL1 = 0;
  m_LambdaL2 = 0;
  m_B0Weight = 1.0;
  m_IsAnalyticalB0 = true;
  m_BasisEnergyPowerDL = 1.0;

  m_MDImage=ScalarImageType::New();
  m_ScaleImage=ScalarImageType::New();

  // m_L1SolverType = FISTA_LS;
  m_L1SolverType = SPAMS;

  m_L2Solver = NULL;
  m_L1FISTASolver = NULL;
  m_L1SpamsSolver = NULL;

  m_IsOriginalBasis = true;

}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  Superclass::VerifyInputParameters();

  utlGlobalException(m_BasisScale<0, "negative scale");
  utlGlobalException(m_LambdaSpherical<0 || m_LambdaRadial<0, "negatvie regularization parameters");
  InputImagePointer input = const_cast<InputImageType *>(this->GetInput());
  utlGlobalException(!IsImageEmpty(m_MDImage) && !itk::VerifyImageInformation(input,m_MDImage,true), "wrong information in m_MDImage");
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  Superclass::GenerateOutputInformation();
  OutputImagePointer outputPtr = this->GetOutput();
  unsigned int numberOfComponentsPerPixel = this->RankToDim();  
  outputPtr->SetNumberOfComponentsPerPixel(numberOfComponentsPerPixel);
}

template< class TInputImage, class TOutputImage >
double
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::ComputeScale(const bool setScale)
{
  m_BasisScale = -1.0;
  return m_BasisScale;
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::SetBasisScale(const double scale)
{
  itkShowPositionThreadedLogger(this->GetDebug());
  double scale_old = m_BasisScale;  
  if (scale>0)
    m_BasisScale = scale;
  else
    this->ComputeScale(true);
  itkDebugMacro("setting m_BasisScale to " << m_BasisScale);

  if (scale>0 && std::fabs((scale_old-m_BasisScale)/m_BasisScale)>1e-8)
    {
    this->Modified();
    this->m_BasisRadialMatrix=MatrixPointer(new MatrixType());
    this->m_BasisMatrix=MatrixPointer(new MatrixType());
    this->m_BasisMatrixForB0=MatrixPointer(new MatrixType());
    }
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_LambdaSpherical = m_LambdaSpherical;
  rval->m_LambdaRadial = m_LambdaRadial;
  rval->m_LambdaL1 = m_LambdaL1;
  rval->m_LambdaL2 = m_LambdaL2;
  rval->m_BasisScale = m_BasisScale;
  rval->m_IsAnalyticalB0 = m_IsAnalyticalB0;
  rval->m_B0Weight = m_B0Weight;
  rval->m_EstimationType = m_EstimationType;
  rval->m_L1SolverType = m_L1SolverType;

  // NOTE: shared_ptr is thread safe, if the data is read in threads (not modified), thus do not need to copy the data block
  rval->m_BasisMatrixForB0 = m_BasisMatrixForB0;
  rval->m_IsOriginalBasis = m_IsOriginalBasis;
  
  rval->m_BasisCombinationMatrix = m_BasisCombinationMatrix;
  rval->m_BasisEnergyDL = m_BasisEnergyDL;
  rval->m_BasisEnergyPowerDL = m_BasisEnergyPowerDL;
  
  rval->m_MDImage = m_MDImage;
  rval->m_ScaleImage = m_ScaleImage;

  if (m_L2Solver)
    rval->m_L2Solver = m_L2Solver->Clone();
  if (m_L1FISTASolver)
    rval->m_L1FISTASolver = m_L1FISTASolver->Clone();
  if (m_L1SpamsSolver)
    rval->m_L1SpamsSolver = m_L1SpamsSolver->Clone();

  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  if (this->m_SamplingSchemeQSpace->GetBVector()->size()>0 && this->m_SamplingSchemeQSpace->GetRadiusVector()->size()==0)
    this->m_SamplingSchemeQSpace->ConvertBVectorToQVector();

  this->VerifyInputParameters();

  
  // create m_LoggerVector for multiple threads
  if (this->GetDebug() && this->GetNumberOfThreads()>1)
    this->CreateLoggerVector();
  
  if (this->GetDebug())
    this->Print(std::cout<<"this BeforeThreadedGenerateData = ");

  std::cout << "Use " << this->GetNumberOfThreads() << " threads!" << std::endl << std::flush;
  if (!IsImageEmpty(this->m_MaskImage))
    std::cout << "Use a mask" << std::endl << std::flush;
  
  if (this->m_EstimationType==Self::LS)
    {
    std::cout << "Use Least Square Estimation" << std::endl << std::flush;
    if (!this->m_L2Solver)
      this->m_L2Solver = L2SolverType::New();
    }
  else if (this->m_EstimationType==Self::L1_2 || this->m_EstimationType==Self::L1_DL )
    {
    if (this->m_EstimationType==Self::L1_2)
      std::cout << "Use L1 Estimation with two lambdas m_LambdaSpherical, m_LambdaRadial (L1_2)" << std::endl << std::flush;
    if (this->m_EstimationType==Self::L1_DL)
      std::cout << "Use L1 Estimation with learned dictionary, one m_LambdaL1" << std::endl << std::flush;
    utlGlobalException(!this->m_L1FISTASolver && !this->m_L1SpamsSolver, "must set a L1 solver first");
    utlGlobalException( (this->m_L1SolverType==Self::FISTA_LS) && !this->m_L1FISTASolver, "m_L1FISTASolver is needed.");
    utlGlobalException( (this->m_L1SolverType==Self::SPAMS) && !this->m_L1SpamsSolver, "m_L1SpamsSolver is needed.");

    if (this->m_L1SolverType==Self::FISTA_LS)
      {
      std::cout << "Use FISTA with least square initialization (FISTA_LS)" << std::endl << std::flush;
      this->m_L1FISTASolver->SetUseL2SolverForInitialization(true);
      }
    else if (this->m_L1SolverType==Self::SPAMS)
      std::cout << "Use SPAMS for weighted lasso (SPAMS)" << std::endl << std::flush;
    }
  else
    utlGlobalException(true, "wrong type");
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::InitializeThreadedLibraries()
{
  Superclass::InitializeThreadedLibraries();
#ifdef DMRITOOL_USE_OPENMP
  // if (this->m_L1SolverType==Self::SPAMS)
    // this->m_L1SpamsSolver->SetNumberOfThreads(this->GetNumberOfThreads()>1?1:-1);
  
  // NOTE: it seems that OMP_NUM_THREAD=1 in SpamsSolver works better than multiple thread, even if this->GetNumberOfThreads()==1, 
  // it is maybe becuase when this->GetNumberOfThreads()==1, we already use multiple-threads for blas.
  if (this->m_L1SolverType==Self::SPAMS)
    this->m_L1SpamsSolver->SetNumberOfThreads(1);
#endif

}


template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierEstimationImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar2(true, m_BasisScale, m_IsOriginalBasis, os<<indent);
  PrintVar2(true, m_IsAnalyticalB0, m_B0Weight, os<<indent);
  if (m_BasisCombinationMatrix->Rows()!=0)
    utl::PrintUtlMatrix(*m_BasisCombinationMatrix, "m_BasisCombinationMatrix", " ", os<<indent);
  if (m_EstimationType==LS)
    os << indent << "Least Square Estimation is used " << std::endl;
  else if (m_EstimationType==L1_2)
    os << indent << "L1 Estimation is used (two lambdas, Laplace-Beltrami like regularization)" << std::endl;
  else if (m_EstimationType==L1_DL)
    os << indent << "L1 Estimation with learned dictionary is used " << std::endl;
  PrintVar4(true, m_LambdaSpherical, m_LambdaRadial, m_LambdaL1, m_LambdaL2, os<<indent);

}

}

#endif 



