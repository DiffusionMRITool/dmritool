/**
 *       @file  itkSpecialFunctionGeneratorGTest.cxx
 *      @brief  
 *     Created  "11-08-2015
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "utlGTest.h"
#include "itkSpecialFunctionGenerator.h"
#include "utl.h"
#include "itkDiffusionTensor.h"

#include "itkSHCoefficientsFit.h"
#include "itkSHCoefficientsRotation.h"
  
typedef utl::NDArray<double,1>  VectorType;
typedef utl_shared_ptr<VectorType> VectorPointer;
typedef utl::NDArray<double,2>  MatrixType;
typedef utl_shared_ptr<utl::NDArray<double,2> >  MatrixPointer;


TEST(itkSpecialFunctionGenerator, BesselJInteger)
{
  EXPECT_NEAR(0.9275303677,utl::BesselJInteger(0,0.5434), 1e-8);
  EXPECT_NEAR(0.2617940621,utl::BesselJInteger(1,0.5434), 1e-8);
  EXPECT_NEAR(-0.2617940621,utl::BesselJInteger(-1,0.5434), 1e-8);
}

TEST(itkSpecialFunctionGenerator, BesselJa)
{
  EXPECT_NEAR(0.9275303677,utl::BesselJa(0,0.5434), 1e-8);
  EXPECT_NEAR(0.2617940621,utl::BesselJa(1,0.5434), 1e-8);
  EXPECT_NEAR(-0.2617940621,utl::BesselJa(-1,0.5434), 1e-8);

  EXPECT_NEAR(0.5596443646,utl::BesselJa(0.5,0.5434), 1e-8);
  EXPECT_NEAR(0.1034236073,utl::BesselJa(1.5,0.5434), 1e-8);
}

TEST(itkSpecialFunctionGenerator, BesselJIntegerPrime)
{
  EXPECT_NEAR(-0.2617940621,utl::BesselJIntegerPrime(0,0.5434), 1e-8);
  EXPECT_NEAR(0.4457599184,utl::BesselJIntegerPrime(1,0.5434), 1e-8);
  EXPECT_NEAR(-0.4457599184,utl::BesselJIntegerPrime(-1,0.5434), 1e-8);
}

TEST(utlMath, GetExpProductLegendreCoef)
{
  EXPECT_NEAR(0.3678794412, utl::GetExpProductLegendreCoef(1.0, 1.0, 0), 1e-8);
  EXPECT_NEAR(27.11263892, utl::GetExpProductLegendreCoef(-3.3, -3.3, 0), 1e-8);

  EXPECT_NEAR(0.1630132586, utl::GetExpProductLegendreCoef(3.3, 1.3, 0), 1e-8);
  EXPECT_NEAR(-0.1710392255, utl::GetExpProductLegendreCoef(3.3, 1.3, 2), 1e-8);
  EXPECT_NEAR(0.05407323817, utl::GetExpProductLegendreCoef(3.3, 1.3, 4), 1e-8);
  EXPECT_NEAR(-0.01048537144, utl::GetExpProductLegendreCoef(3.3, 1.3, 6), 1e-8);
  EXPECT_NEAR(0.001468833732, utl::GetExpProductLegendreCoef(3.3, 1.3, 8), 1e-8);
  EXPECT_NEAR(-0.0001610249754, utl::GetExpProductLegendreCoef(3.3, 1.3, 10), 1e-8);
  EXPECT_NEAR(0.00001449640342, utl::GetExpProductLegendreCoef(3.3, 1.3, 12), 1e-8);

  EXPECT_NEAR(0.08720854874, utl::GetExpProductLegendreCoef(1.3, 3.3, 0), 1e-8);
  EXPECT_NEAR(0.1294597112, utl::GetExpProductLegendreCoef(1.3, 3.3, 2), 1e-8);
  EXPECT_NEAR(0.04519205145, utl::GetExpProductLegendreCoef(1.3, 3.3, 4), 1e-8);
  EXPECT_NEAR(0.009186405117, utl::GetExpProductLegendreCoef(1.3, 3.3, 6), 1e-8);
  EXPECT_NEAR(0.001322836219, utl::GetExpProductLegendreCoef(1.3, 3.3, 8), 1e-8);
  EXPECT_NEAR(0.0001476644949, utl::GetExpProductLegendreCoef(1.3, 3.3, 10), 1e-8);
  EXPECT_NEAR(0.00001346433608, utl::GetExpProductLegendreCoef(1.3, 3.3, 12), 1e-8);
}

TEST(itkSpecialFunctionGenerator, GetSymmetricTensorSHCoef)
{
  itk::InitializeLUTExp();

  MatrixPointer gradT4 = utl::ReadGrad<double>(4,DIRECTION_NODUPLICATE,CARTESIAN_TO_CARTESIAN);
  int N = gradT4->Rows();
  int lMax = 12;
  double b = 4000, e1=1.7e-3, e2=0.2e-3;
  std::vector<double> bVector(N, b);
  MatrixPointer shMatrix = utl::ComputeSHMatrix(lMax, *gradT4, CARTESIAN_TO_SPHERICAL);

  typedef itk::DiffusionTensor<double> TensorType;
  std::vector<TensorType>  tensors(2);
  std::vector<double> vec(3,1);
  MatrixType matEigen(3,3);

  tensors[0].Fill(0.0); 
  vec[0]=e1, vec[1]=e2, vec[2]=e2;
  tensors[0].SetEigenValues(vec);
  tensors[1].Fill(0.0); 
  vec[0]=e2, vec[1]=e2, vec[2]=e1;
  tensors[1].SetEigenValues(vec);

  VectorType dwiVec(N), dwiVecRe(N); 
  VectorPointer dwiSHCoef;
  double rr=1, theta, phi;
  for ( int i = 0; i < tensors.size(); ++i ) 
    {
    TensorType tensor = tensors[i];
    tensor.GetDWISamples(dwiVec, *gradT4, bVector);
    tensor.GetEigenValuesVectorsAnalytic(vec, matEigen);

    utl::cartesian2Spherical(matEigen(0,2), matEigen(1,2), matEigen(2,2), rr, theta, phi);

    dwiSHCoef = utl::GetSymmetricTensorSHCoef ( b, e1, e2, lMax, theta, phi );
    dwiVecRe = (*shMatrix)* (*dwiSHCoef);
    VectorType diff = dwiVec-dwiVecRe;
    // utl::PrintUtlVector(diff, "diff");
    // std::cout << "diff = " << diff.GetTwoNorm()/dwiVec.GetTwoNorm() << std::endl << std::flush;
    EXPECT_NEAR_VECTOR(dwiVec, dwiVecRe, N, 0.005);
    EXPECT_NEAR(0.0, diff.GetTwoNorm()/dwiVec.GetTwoNorm(), 0.005);
    }
}

TEST(itkSpecialFunctionGenerator, Hyperg1F1)
{
  EXPECT_NEAR(8.22631388275362, utl::Hyperg1F1(0.5,1.5,4.0), 1e-8);
  EXPECT_NEAR(54.5981500331442, utl::Hyperg1F1(1.5,1.5,4.0), 1e-8);
  EXPECT_NEAR(200.193216788196, utl::Hyperg1F1(2.5,1.5,4.0), 1e-8);
}


