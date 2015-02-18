/**
 *       @file  itkDWISingleVoxelGenerator.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkDWISingleVoxelGenerator_hxx
#define __itkDWISingleVoxelGenerator_hxx

#include "itkDWISingleVoxelGenerator.h"
#include "itkDiffusionTensor.h"
#include "utlRotationMatrixFromVectors.h"

namespace itk
{

template <class TOutputImage, class TScalarImage>
DWISingleVoxelGenerator<TOutputImage, TScalarImage>
::DWISingleVoxelGenerator() 
: Superclass(), 
  m_StoredOrientationMatrix(new MatrixType())
{
  m_RandomType = FIXED;
}

template <class TOutputImage, class TScalarImage>
void DWISingleVoxelGenerator<TOutputImage, TScalarImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar1(true, m_RandomType, os<<indent);
  if (m_StoredOrientationMatrix->Rows()>0)
    utl::PrintUtlMatrix(*m_StoredOrientationMatrix, "m_StoredOrientationMatrix");
}

template <class TOutputImage, class TScalarImage>
typename LightObject::Pointer
DWISingleVoxelGenerator<TOutputImage, TScalarImage>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_DiffusionParameterValues = m_DiffusionParameterValues;
  rval->m_RandomType = m_RandomType;
  rval->m_StoredOrientationMatrix = m_StoredOrientationMatrix;
  return loPtr;
}

template <class TOutputImage, class TScalarImage>
void DWISingleVoxelGenerator<TOutputImage, TScalarImage>
::Initialization()
{
  Superclass::Initialization(); 
  if (m_RandomType==UNIFORM && m_StoredOrientationMatrix->Rows())
    {
    for ( int i = 0; i < TOutputImage::ImageDimension; i += 1 ) 
      this->m_OutputSize[i] = 1.0;
    this->m_OutputSize[0] = m_StoredOrientationMatrix->Rows();
    }
}

template <class TOutputImage, class TScalarImage>
void DWISingleVoxelGenerator<TOutputImage, TScalarImage>
::GenerateData()
{
  utlShowPosition(this->GetDebug());
  this->Initialization();

  bool isOutputPeak = this->m_MaxNumberOfPeaks>0;
  
  MatrixPointer qSpaceOrientationMatrix = this->m_SamplingSchemeQSpace->GetOrientationsCartesian();
  STDVectorPointer bVector = this->m_SamplingSchemeQSpace->GetBVector();
  MatrixPointer rSpaceOrientationMatrix = this->m_SamplingSchemeRSpace->GetOrientationsCartesian();
  STDVectorPointer rVector = this->m_SamplingSchemeRSpace->GetRadiusVector();
  double tau = this->m_SamplingSchemeQSpace->GetTau();

  // allocate all outputs
  this->AllocateOutputs();

  int numberOfDiffusivityComponentsPerParameterSet = -1;
  unsigned int numberOfOrientationComponentsPerParameterSet = 3;
  if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS)
    {
    numberOfDiffusivityComponentsPerParameterSet=2;
    numberOfOrientationComponentsPerParameterSet = 3;
    }
  else if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS)
    {
    numberOfDiffusivityComponentsPerParameterSet=2;
    numberOfOrientationComponentsPerParameterSet = 2;
    }
  else if (this->m_ModelType==Superclass::TENSOR_IN_EULER_ANGLES)
    {
    numberOfDiffusivityComponentsPerParameterSet=3;
    numberOfOrientationComponentsPerParameterSet = 3;
    utlGlobalException(true, "TODO");
    }
  else if (this->m_ModelType==Superclass::CYLINDER_SPHERICAL_MODEL)
    {
    numberOfOrientationComponentsPerParameterSet = 2;
    utlGlobalException(true, "TODO");
    }

  int numberOfComponentsPerParameterSet = numberOfDiffusivityComponentsPerParameterSet + numberOfOrientationComponentsPerParameterSet+1;
  // utlSAPrint(m_DiffusionParameterValues.size())(numberOfDiffusivityComponentsPerParameterSet)(numberOfOrientationComponentsPerParameterSet)(numberOfComponentsPerParameterSet);
  utlSAGlobalException(!utl::IsInt(m_DiffusionParameterValues.size()*1.0 / numberOfComponentsPerParameterSet))(m_DiffusionParameterValues.size())(numberOfComponentsPerParameterSet).msg("wrong diffusion parameter size");
  int numberOfParameterSets = m_DiffusionParameterValues.size() /numberOfComponentsPerParameterSet;
  int numberOfComponentsPerPixel4Peak = isOutputPeak?PeakContainerHelper::GetDimension(this->m_PeakType,this->m_MaxNumberOfPeaks):-1;
  utlSAGlobalException(isOutputPeak && numberOfParameterSets>PeakContainerHelper::GetNumberOfPeaks(this->m_PeakType,numberOfComponentsPerPixel4Peak))(numberOfParameterSets)(numberOfComponentsPerPixel4Peak).msg("wrong peak size");

  std::vector<DiffusionTensor<double> > tensorVec;
  std::vector<double> weightVec;

  DiffusionTensor<double> tensor;
  Vector<PrecisionType, 3> principalDirection(0.0);
  double norm, partialVolumeWeight;
  Vector<PrecisionType, 3> e1;
  e1[0]=1.0, e1[1]=0.0, e1[2]=0.0;
  PrecisionType sumPartialVolumeWeight=0.0;
  STDVectorType diffusivities(numberOfDiffusivityComponentsPerParameterSet);
  // utl::PrintVector(m_DiffusionParameterValues, "m_DiffusionParameterValues");
  // utlSAPrint(numberOfParameterSets)(numberOfOrientationComponentsPerParameterSet)(numberOfComponentsPerParameterSet);

  OutputImagePixelType pixelPeakInput, pixelPeak;
  int numberOfPeaks=0;
  if (isOutputPeak)
    {
    pixelPeakInput.SetSize(numberOfComponentsPerPixel4Peak);
    pixelPeakInput.Fill(0.0);
    if (this->m_PeakType==NXYZ || this->m_PeakType==NXYZV)
      pixelPeakInput[0] = numberOfParameterSets;
    pixelPeak.SetSize(numberOfComponentsPerPixel4Peak);
    }
  for ( unsigned int s=0; s<numberOfParameterSets; s++ )
    {
    for ( unsigned int d=0; d<numberOfOrientationComponentsPerParameterSet; d++ )
      principalDirection[d] = m_DiffusionParameterValues[s* numberOfComponentsPerParameterSet + d];

    norm = principalDirection.GetNorm();
    if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS || this->m_ModelType==Superclass::CYLINDER_SPHERICAL_MODEL)
      {
      principalDirection *= vnl_math::pi/180.0;
      double theta=principalDirection[0], phi=principalDirection[1];
      utl::spherical2Cartesian(1.0,theta,phi,principalDirection[0],principalDirection[1],principalDirection[2]);
      }
    else
      {
      for ( unsigned int d=0; d<numberOfOrientationComponentsPerParameterSet; d++ )
        principalDirection[d] /= norm;
      }

    for ( unsigned int d=0; d<numberOfDiffusivityComponentsPerParameterSet; d++ )
      {
      diffusivities[d] = m_DiffusionParameterValues[s*numberOfComponentsPerParameterSet + numberOfOrientationComponentsPerParameterSet + d ];
      }
    std::sort(diffusivities.begin(), diffusivities.end(), std::greater<double>());

    partialVolumeWeight = m_DiffusionParameterValues[s*numberOfComponentsPerParameterSet + numberOfOrientationComponentsPerParameterSet + numberOfDiffusivityComponentsPerParameterSet];
    sumPartialVolumeWeight += partialVolumeWeight;
    weightVec.push_back(partialVolumeWeight);

    // utlPrintVar2(true, s, PeakContainerHelper::GetString(this->m_PeakType));
    // std::cout << "principalDirection = " << principalDirection << std::endl << std::flush;
    // utl::PrintVector(diffusivities, "diffusivities");

    std::vector<double> vec(3,1);
    if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS || this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS)
      vec[0]=diffusivities[0], vec[1]=diffusivities[1], vec[2]=diffusivities[1];
    else if (this->m_ModelType==Superclass::TENSOR_IN_EULER_ANGLES)
      vec[0]=diffusivities[0], vec[1]=diffusivities[1], vec[2]=diffusivities[2];
    tensor.Fill(0.0);
    tensor.SetEigenValues(vec);
    // utlPrintVar2(true, "before", tensor);
    Matrix<double,3,3> rotation;
    utl::RotationMatrixFromVectors<Vector<double,3>, Matrix<double,3> >(e1, principalDirection, rotation); 
    tensor.Rotate(rotation);

    // utlPrintVar2(true, "after", tensor);
    tensorVec.push_back(tensor);
    if (isOutputPeak)
      {
      PeakContainerHelper::SetPeak<Vector<double,3>, OutputImagePixelType >(principalDirection, pixelPeakInput, s, this->m_PeakType);
      if (this->m_PeakType==NXYZV || this->m_PeakType==XYZV)
        PeakContainerHelper::SetPeakValue<OutputImagePixelType >(partialVolumeWeight, pixelPeakInput, s, this->m_PeakType);
      numberOfPeaks++;
      }
    }

  if (isOutputPeak && (this->m_PeakType==NXYZ || this->m_PeakType==NXYZV))
    pixelPeakInput[0] = numberOfPeaks;

  utlSAException(sumPartialVolumeWeight<1e-9)(sumPartialVolumeWeight).msg("wrong partialVolumeWeight");
  for ( int s = 0; s < numberOfParameterSets; s += 1 ) 
    weightVec[s] /= sumPartialVolumeWeight;


  // typename OutputImageType::Pointer outputPtr = this->GetOutput();
  OutputImagePointer outputDWI=this->GetDWIImage(), outputODF=this->GetODFImage(), outputEAP=this->GetEAPImage(), outputPeak=this->GetPeakImage();
  ScalarImagePointer b0Image = this->GetB0Image(), outputRTO = this->GetRTOImage(), outputMSD = this->GetMSDImage();
  OutputImagePixelType pixelDWI, tempPixelDWI, pixelEAP, tempPixelEAP, pixelODF, tempPixelODF;
  double pixelRTO=0, tempPixelRTO=0, pixelMSD=0, tempPixelMSD=0;  
  int numberOfComponentsPerPixel4DWI = this->GetNumberOfQSpaceSamples();
  int numberOfComponentsPerPixel4EAP = this->GetNumberOfRSpaceSamples();
  int numberOfComponentsPerPixel4ODF = this->GetNumberOfRSpaceSamples();
  OutputImagePixelType zeroPixel;
  if (this->m_IsOutputDWI)
    {
    pixelDWI.SetSize( numberOfComponentsPerPixel4DWI );
    pixelDWI.Fill( 0 );
    tempPixelDWI.SetSize( numberOfComponentsPerPixel4DWI );
    tempPixelDWI.Fill( 0 );
    }
  if (this->m_IsOutputEAP)
    {
    pixelEAP.SetSize( numberOfComponentsPerPixel4EAP );
    pixelEAP.Fill( 0 );
    tempPixelEAP.SetSize( numberOfComponentsPerPixel4EAP );
    tempPixelEAP.Fill( 0 );
    }
  if (this->m_IsOutputODF)
    {
    pixelODF.SetSize( numberOfComponentsPerPixel4ODF );
    pixelODF.Fill( 0 );
    tempPixelODF.SetSize( numberOfComponentsPerPixel4ODF );
    tempPixelODF.Fill( 0 );
    }

  ImageRegionIteratorWithIndex<OutputImageType> outputDWIIt, outputODFIt, outputEAPIt, outputPeakIt;
  ImageRegionIteratorWithIndex<ScalarImageType> b0It, outputRTOIt, outputMSDIt;
  if ( this->m_IsOutputDWI )
    {
  // outputDWI->Print(std::cout<<"outputDWI");
    outputDWIIt = ImageRegionIteratorWithIndex<OutputImageType>( outputDWI, outputDWI->GetRequestedRegion() );
    b0It = ImageRegionIteratorWithIndex<ScalarImageType>( b0Image, b0Image->GetRequestedRegion() );
    }
  if ( this->m_IsOutputODF )
    outputODFIt = ImageRegionIteratorWithIndex<OutputImageType>( outputODF, outputODF->GetRequestedRegion() );
  if ( this->m_IsOutputEAP )
    outputEAPIt = ImageRegionIteratorWithIndex<OutputImageType>( outputEAP, outputEAP->GetRequestedRegion() );
  if ( isOutputPeak )
    outputPeakIt = ImageRegionIteratorWithIndex<OutputImageType>( outputPeak, outputPeak->GetRequestedRegion() );
  if ( this->m_IsOutputRTO )
    outputRTOIt = ImageRegionIteratorWithIndex<ScalarImageType>( outputRTO, outputRTO->GetRequestedRegion() );
  if ( this->m_IsOutputMSD )
    outputMSDIt = ImageRegionIteratorWithIndex<ScalarImageType>( outputMSD, outputMSD->GetRequestedRegion() );

  
  b0It.GoToBegin();
  outputDWIIt.GoToBegin();
  outputODFIt.GoToBegin();
  outputEAPIt.GoToBegin();
  outputPeakIt.GoToBegin();
  OutputImagePixelType pixel_noise=pixelDWI;
  Matrix<double,3,3> rotation;
  int indexGrad = 0;
  std::vector<double> peak;
  while ( !outputDWIIt.IsAtEnd() )
    {
    if (this->m_IsOutputDWI)
      pixelDWI.Fill( 0 );
    if (this->m_IsOutputODF)
      pixelODF.Fill( 0 );
    if (this->m_IsOutputEAP)
      pixelEAP.Fill( 0 );
    pixelMSD=0, pixelRTO=0;
    if (isOutputPeak)
      pixelPeak = pixelPeakInput;
      
    if (m_RandomType==UNIFORM)
      {
      std::vector<double> vv(3,0);
      if (m_StoredOrientationMatrix->Rows()>0)
        {
        vv[0] = (*m_StoredOrientationMatrix)(indexGrad,0);
        vv[1] = (*m_StoredOrientationMatrix)(indexGrad,1);
        vv[2] = (*m_StoredOrientationMatrix)(indexGrad,2);
        indexGrad++;
        }
      else
        vv = utl::RandomPointInSphere(false);
      for ( int i = 0; i < 3; i += 1 ) 
        principalDirection[i] = vv[i];
      utl::RotationMatrixFromVectors<Vector<double,3>, Matrix<double,3> >(e1, principalDirection, rotation); 
      // utl::PrintContainer(&e1[0], &e1[3], "e1");
      // utl::PrintContainer(&principalDirection[0], &principalDirection[3], "principalDirection");
      // utl::PrintMatrix(rotation, 3,3, "rotation");
      }

    for ( int s = 0; s < numberOfParameterSets; s += 1 ) 
      {
      tensor = tensorVec[s];
      partialVolumeWeight = weightVec[s];
      if (isOutputPeak)
        peak = PeakContainerHelper::GetPeak(pixelPeakInput, s, this->m_PeakType);
      // tensor.Print(std::cout<<"tensor 0=");

      if (m_RandomType==UNIFORM)
        {
        tensor.Rotate(rotation);
        if (isOutputPeak)
          {
          peak = PeakContainerHelper::GetPeak(pixelPeakInput, s, this->m_PeakType);
          vnl_vector<double> peakVnl = utl::StdVectorToVnlVector(peak);
          vnl_vector<double> peakRotated = rotation.GetVnlMatrix()*peakVnl;
          PeakContainerHelper::SetPeak<vnl_vector<double>, OutputImagePixelType >(peakRotated, pixelPeak, s, this->m_PeakType);
          }
        }
      // tensor.Print(std::cout<<"tensor 1=");

      if ( this->m_IsOutputDWI )
        {
        tensor.GetDWISamples(tempPixelDWI, *qSpaceOrientationMatrix, *bVector);
        pixelDWI += tempPixelDWI*this->m_B0Scale*partialVolumeWeight;
        }
      if ( this->m_IsOutputODF )
        {
        tensor.GetODFSamples(tempPixelODF, *rSpaceOrientationMatrix,this->m_ODFOrder, false);
        pixelODF += tempPixelODF*partialVolumeWeight;
        }
      if ( this->m_IsOutputEAP)
        {
        tensor.GetEAPSamples(tempPixelEAP, *rSpaceOrientationMatrix, *rVector, tau);
        pixelEAP += tempPixelEAP*this->m_B0Scale*partialVolumeWeight;
        }
      if ( this->m_IsOutputRTO)
        {
        tempPixelRTO = tensor.GetReturnToOrigin(tau);
        pixelRTO += tempPixelRTO*partialVolumeWeight;
        }
      if ( this->m_IsOutputMSD)
        {
        tempPixelMSD = tensor.GetMeanSquaredDisplacement(tau);
        pixelMSD += tempPixelMSD*partialVolumeWeight;
        }
      }

    if (this->m_IsOutputODF && this->m_ODFOrder!=2)
      {
      PrecisionType normFactor = 4*vnl_math::pi*utl::GetSumOfVector<OutputImagePixelType>(pixelODF,pixelODF.Size()) / pixelODF.Size();
      if (normFactor!=0)
        pixelODF /= normFactor;
      }

    if ( this->m_IsOutputDWI )
      {
      if ( (this->m_NoiseSigma>0 || this->m_SNR>0)  && pixelDWI.GetSquaredNorm()>0 )
        {
        PrecisionType sigma_real = -1;
        if ( this->m_SNR > 0 && this->m_NoiseSigma <= 0)
          {
          double mean = 0;
          for ( unsigned int k=0; k<numberOfComponentsPerPixel4DWI; k++ )
            mean += pixelDWI[k];
          mean /= numberOfComponentsPerPixel4DWI;
          sigma_real = mean/this->m_SNR;
          }
        else if ( this->m_SNR <= 0 && this->m_NoiseSigma > 0)
          sigma_real = this->m_NoiseSigma;
        pixelDWI = utl::AddNoise<OutputImagePixelType>(pixelDWI, pixelDWI.GetSize(), sigma_real, true);
        }
      outputDWIIt.Set(pixelDWI);
      b0It.Set(this->m_B0Scale);
      ++b0It;
      ++outputDWIIt;
      }

    if ( this->m_IsOutputODF )
      {
      outputODFIt.Set(pixelODF);
      ++outputODFIt;
      }

    if ( this->m_IsOutputEAP )
      {
      outputEAPIt.Set(pixelEAP);
      ++outputEAPIt;
      }
    if ( isOutputPeak )
      {
      outputPeakIt.Set(pixelPeak);
      ++outputPeakIt;
      }
    if ( this->m_IsOutputRTO )
      {
      outputRTOIt.Set(pixelRTO);
      ++outputRTOIt;
      }
    if ( this->m_IsOutputMSD )
      {
      outputMSDIt.Set(pixelMSD);
      ++outputMSDIt;
      }
    }
}

}


#endif 



