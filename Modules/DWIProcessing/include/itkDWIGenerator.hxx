/*=========================================================================

 Program:   DWI Generator

 Copyright (c) Pew-Thian Yap, Jian Cheng. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkDWIGenerator_hxx
#define __itkDWIGenerator_hxx

#include "itkDWIGenerator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkDiffusionTensor.h"
#include "utlRotationMatrixFromVectors.h"
#include "itksys/SystemTools.hxx"

#include "utl.h"

namespace itk
{

template <class TOutputImage, class TScalarImage>
DWIGenerator<TOutputImage, TScalarImage>
::DWIGenerator() : Superclass()
{
}

template <class TOutputImage, class TScalarImage>
DWIGenerator<TOutputImage, TScalarImage>
::~DWIGenerator()
{
}

template <class TOutputImage, class TScalarImage>
void DWIGenerator<TOutputImage, TScalarImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " << m_FileName << std::endl;
}

template <class TOutputImage, class TScalarImage>
typename LightObject::Pointer
DWIGenerator<TOutputImage, TScalarImage>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  
  rval->m_BackgroundDiffusionParameterValues = m_BackgroundDiffusionParameterValues;
  rval->m_FileName = m_FileName;
  return loPtr;
}

template <class TOutputImage, class TScalarImage>
void DWIGenerator<TOutputImage, TScalarImage>
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
  if (this->GetDebug())
    std::cout << "tau = " << tau << std::endl << std::flush;
  
  int numberOfComponentsPerPixel4Peak = isOutputPeak?PeakContainerHelper::GetDimension(this->m_PeakType,this->m_MaxNumberOfPeaks):-1;
  // utlSAGlobalException(isOutputPeak && numberOfParameterSets>PeakContainerHelper::GetNumberOfPeaks(this->m_PeakType,numberOfComponentsPerPixel4Peak))(numberOfParameterSets)(numberOfComponentsPerPixel4Peak).msg("wrong peak size");
  OutputImagePixelType pixelPeak;
  if (isOutputPeak)
    {
    pixelPeak.SetSize(numberOfComponentsPerPixel4Peak);
    }
      
  // Parse parameter file
  // This section should be rewritten to support generic image dimension
//  char * pch;
//  unsigned int idx;
  unsigned int ndims = 0;
//  char tempLine[256];
  std::string line, extractedLine;
  double scale = 1;
  double sigma = -1;
  double snr = -1;
  double backgroundScale = -1;

  DiffusionParameterValuesType diffusionParameterValues;
  DiffusionParameterContainerType diffusionParameterContainer;
  diffusionParameterContainer.clear();
  DiffusionParameterValuesType backgroundDiffusionParameterValues;
  backgroundDiffusionParameterValues.clear();

  // Read Header
  utlGlobalException(!utl::IsFileExist(m_FileName), "The file does not exist! m_FileName="<<m_FileName);

  std::vector<std::vector<std::string> > stringMatrix;
  utl::ReadLines(m_FileName, stringMatrix, " ");

  for ( int i = 0; i < stringMatrix.size(); i += 1 ) 
    {
    if (stringMatrix[i][0]=="NDims" && stringMatrix[i][1]=="=")
      std::istringstream ( stringMatrix[i][2] ) >> ndims;

    if (stringMatrix[i][0]=="DimSize" && stringMatrix[i][1]=="=")
      {
      for ( int kk = 0; kk < stringMatrix[i].size()-2; kk += 1 ) 
        std::istringstream ( stringMatrix[i][kk+2] ) >> this->m_OutputSize[kk];
      }

    if (stringMatrix[i][0]=="ElementSpacing" && stringMatrix[i][1]=="=")
      {
      for ( int kk = 0; kk < stringMatrix[i].size()-2; kk += 1 ) 
        std::istringstream ( stringMatrix[i][kk+2] ) >> this->m_OutputSpacing[kk];
      }

    // if (stringMatrix[i][0]=="BValues" && stringMatrix[i][1]=="=")
    //   {
    //   if (utl::IsNumber(stringMatrix[i][2])) // float numbers
    //     {
    //     bVector->SetSize(stringMatrix[i].size()-2);
    //     for ( int kk = 0; kk < stringMatrix[i].size()-2; kk += 1 ) 
    //       std::istringstream ( stringMatrix[i][kk+2] ) >> (*bVector)[kk];
    //     }
    //   else // bVector in txt file
    //     {
    //     utlException(stringMatrix[i].size()!=3, "need to set only one txt file");
    //     std::string ext, file;
    //     utl::getFileExtension(stringMatrix[i][2], ext, file);
    //     utlException(ext!="txt", "wrong bVector file");
    //     }
    //   }

    if (stringMatrix[i][0]=="Scale" && stringMatrix[i][1]=="=")
      std::istringstream ( stringMatrix[i][2] ) >> this->m_B0Scale;

    if (stringMatrix[i][0]=="ModelType" && stringMatrix[i][1]=="=")
      {
      if (stringMatrix[i][2]=="SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS") 
        this->m_ModelType = Superclass::SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS;
      else if (stringMatrix[i][2]=="SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS")
        this->m_ModelType = Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS;
      else if (stringMatrix[i][2]=="TENSOR_IN_EULER_ANGLES") 
        {
        this->m_ModelType = Superclass::TENSOR_IN_EULER_ANGLES;
        utlGlobalException(true, "TODO");
        }
      else if (stringMatrix[i][2]=="CYLINDER_SPHERICAL_MODEL") 
        this->m_ModelType = Superclass::CYLINDER_SPHERICAL_MODEL;
      else
        utlGlobalException(true, "wrong model type");
      }

    if (stringMatrix[i][0]=="BackgroundScale" && stringMatrix[i][1]=="=")
      std::istringstream ( stringMatrix[i][2] ) >> backgroundScale;

    if (stringMatrix[i][0]=="DiffusionParameters" && stringMatrix[i][1]=="=")
      {
      diffusionParameterValues.clear();
      for ( int kk = 0; kk < stringMatrix[i].size()-2; kk += 1 ) 
        diffusionParameterValues.push_back( atof(stringMatrix[i][kk+2].c_str()) );

      diffusionParameterContainer.push_back( diffusionParameterValues );
      }

    if (stringMatrix[i][0]=="BackgroundDiffusionParameters" && stringMatrix[i][1]=="=")
      {
      backgroundDiffusionParameterValues.clear();
      for ( int kk = 0; kk < stringMatrix[i].size()-2; kk += 1 ) 
        backgroundDiffusionParameterValues.push_back( atof(stringMatrix[i][kk+2].c_str()) );
      }

    if (stringMatrix[i][0]=="RicianNoiseSigma" && stringMatrix[i][1]=="=")
      std::istringstream ( stringMatrix[i][2] ) >> sigma;

    if (stringMatrix[i][0]=="SNR" && stringMatrix[i][1]=="=")
      std::istringstream ( stringMatrix[i][2] ) >> snr;
    }

  if ( this->m_NoiseSigma >= 0 )
    {
    sigma = this->m_NoiseSigma;
    }
 
  if ( this->m_SNR >= 0 )
    {
    snr = this->m_SNR;
    }

  utlException(snr>0 && sigma>0, "snr and sigma can not be set simultaneously. snr="<<snr<<", sigma="<<sigma);
  
  scale = this->m_B0Scale;
  if (backgroundScale<=0)
    backgroundScale = this->m_B0Scale;
       
  // NOTE: for VectorImage, the spacing is determined by the size
  for ( int i = this->m_OutputSize.GetSizeDimension()-1; i >= 0; i -= 1 ) 
    {
    if (this->m_OutputSize[i]==1)
      this->m_OutputSpacing[i] = 1;
    else
      break;
    }

  // allocate all outputs
  this->AllocateOutputs();

  // typename OutputImageType::Pointer outputPtr = this->GetOutput();
  OutputImagePointer outputDWI=this->GetDWIImage(), outputODF=this->GetODFImage(), outputEAP=this->GetEAPImage(), outputPeak=this->GetPeakImage();
  ScalarImagePointer b0Image = this->GetB0Image(), outputRTO = this->GetRTOImage(), outputMSD = this->GetMSDImage();
  OutputImagePixelType pixelDWI, tempPixelDWI, pixelEAP, tempPixelEAP, pixelODF, tempPixelODF;
  double pixelRTO=0, tempPixelRTO=0, pixelMSD=0, tempPixelMSD=0;  
  unsigned int numberOfComponentsPerPixel4DWI = this->GetNumberOfQSpaceSamples();
  unsigned int numberOfComponentsPerPixel4EAP = this->GetNumberOfRSpaceSamples();
  unsigned int numberOfComponentsPerPixel4ODF = this->GetNumberOfRSpaceSamples();
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

  OutputImageIndexType index;

    
  // Convenient variables
  int numberOfParameterSets = -1;
  unsigned int numberOfOrientationComponentsPerParameterSet = 3;
  int numberOfComponentsPerParameterSet = -1; // Added partial volume fractions
  int numberOfDiffusivityComponentsPerParameterSet = -1;
  
  bool isTensorModel=false;
  DiffusionTensor<double> tensor;

  bool isCylinderModel=false;
  typename CylinderModelType::PointType cylinderAxis;
  typename CylinderModelType::VectorPointer tempPixelDWIVec, tempPixelODFVec, tempPixelEAPVec;

  if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS)
    {
    numberOfDiffusivityComponentsPerParameterSet=2;
    numberOfOrientationComponentsPerParameterSet = 3;
    numberOfComponentsPerParameterSet = 6;
    isTensorModel = true;
    }
  else if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS)
    {
    numberOfDiffusivityComponentsPerParameterSet=2;
    numberOfOrientationComponentsPerParameterSet = 2;
    numberOfComponentsPerParameterSet = 5;
    isTensorModel = true;
    }
  else if (this->m_ModelType==Superclass::TENSOR_IN_EULER_ANGLES)
    {
    numberOfDiffusivityComponentsPerParameterSet=3;
    numberOfOrientationComponentsPerParameterSet = 3;
    numberOfComponentsPerParameterSet = 7;
    isTensorModel = true;
    utlGlobalException(true, "TODO");
    }
  else if (this->m_ModelType==Superclass::CYLINDER_SPHERICAL_MODEL)
    {
    numberOfDiffusivityComponentsPerParameterSet=0;
    numberOfOrientationComponentsPerParameterSet = 2;
    numberOfComponentsPerParameterSet = 3;
    isCylinderModel = true;
    }

  utlGlobalException(!isTensorModel && !isCylinderModel, "wrong model type");
  
  if (isCylinderModel)
    {
    this->m_CylinderModel->SetDebug(this->GetDebug());
    this->m_CylinderModel->BuildTable();
    this->m_CylinderModel->SetSamplingSchemeQSpace(this->m_SamplingSchemeQSpace);
    this->m_CylinderModel->SetSamplingSchemeRSpace(this->m_SamplingSchemeRSpace);
    }
  
  Vector<PrecisionType, 3> principalDirection(0.0);
  STDVectorType diffusivities;
  if (numberOfDiffusivityComponentsPerParameterSet>0)
    diffusivities.resize(numberOfDiffusivityComponentsPerParameterSet);
  Vector<PrecisionType, 3> e1;
  e1[0]=1.0, e1[1]=0.0, e1[2]=0.0;
  double partialVolumeWeight = 0;

  double norm;
//  double temp;
  
  // std::cout << "scale = " << scale << std::endl;
  // std::cout << "backgroundScale = " << backgroundScale << std::endl;
  // utlPrintVar1(true, numberOfDiffusivityComponentsPerParameterSet);

  // Simulated DWI data
  for ( unsigned int k=0; k<diffusionParameterContainer.size(); k++ )
    {
    diffusionParameterValues = diffusionParameterContainer[k];

    for ( unsigned int dim=0; dim<ndims; dim++ )
      {
      index[dim] = diffusionParameterValues[dim];
      }
    if (this->GetDebug())
      {
      std::cout << "index = " << index << std::endl << std::flush;
      utl::PrintVector(diffusionParameterValues, "diffusionParameterValues");
      }
    
    // 3 elements for directions, 2 elements for diffusivities
    utlSAGlobalException(!utl::IsInt( (1.0*diffusionParameterValues.size()-ndims)/numberOfComponentsPerParameterSet) )(diffusionParameterValues.size())(ndims)(numberOfComponentsPerParameterSet).msg("wrong diffusion parameter size");
    numberOfParameterSets = ( diffusionParameterValues.size() - ndims ) /numberOfComponentsPerParameterSet;
    if (this->m_IsOutputDWI)
      pixelDWI.Fill( 0 );
    if (this->m_IsOutputEAP)
      pixelEAP.Fill( 0 );
    if (this->m_IsOutputODF)
      pixelODF.Fill( 0 );
    if (isOutputPeak)
      pixelPeak.Fill( 0 );
    pixelMSD=0, pixelRTO=0;

    // std::cout << "k = "<< k << ", pixel = " << pixel << std::endl << std::flush;
    // utlPrintVar4(true, numberOfParameterSets, diffusionParameterValues.size(), ndims, numberOfComponentsPerParameterSet);
    PrecisionType sumPartialVolumeWeight=0.0;
    int numberOfPeaks=0;
    for ( unsigned int s=0; s<numberOfParameterSets; s++ )
      {
      for ( unsigned int d=0; d<numberOfOrientationComponentsPerParameterSet; d++ )
        principalDirection[d] = diffusionParameterValues[ndims + s* numberOfComponentsPerParameterSet + d];

      norm = principalDirection.GetNorm();
      if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS || this->m_ModelType==Superclass::CYLINDER_SPHERICAL_MODEL)
        {
        principalDirection *= vnl_math::pi/180.0;
        double theta=principalDirection[0], phi=principalDirection[1];
        utl::spherical2Cartesian(1.0,theta,phi,principalDirection[0],principalDirection[1],principalDirection[2]);
        }
      else
        {
        for ( unsigned int d=0; d<3; d++ )
          principalDirection[d] /= norm;
        }

      for ( unsigned int d=0; d<numberOfDiffusivityComponentsPerParameterSet; d++ )
        {
        diffusivities[d] = diffusionParameterValues[ndims + s 
            * numberOfComponentsPerParameterSet 
            + numberOfOrientationComponentsPerParameterSet + d ];
        }
      if (diffusivities.size()>0)
        std::sort(diffusivities.begin(), diffusivities.end(), std::greater<double>());
      
      partialVolumeWeight = diffusionParameterValues[ndims + s 
          * numberOfComponentsPerParameterSet 
          + numberOfOrientationComponentsPerParameterSet 
          + numberOfDiffusivityComponentsPerParameterSet];

      sumPartialVolumeWeight += partialVolumeWeight;
  
      // utlPrintVar4(true, k, index[0], index[1], index[2]);
      // utl::PrintVector(diffusivities, "diffusivities");

      if (isOutputPeak)
        {
        PeakContainerHelper::SetPeak<Vector<double,3>, OutputImagePixelType >(principalDirection, pixelPeak, s, this->m_PeakType);
        if (this->m_PeakType==NXYZV || this->m_PeakType==XYZV)
          PeakContainerHelper::SetPeakValue<OutputImagePixelType >(partialVolumeWeight, pixelPeak, s, this->m_PeakType);
        numberOfPeaks++;
        }

      if (isTensorModel)
        {
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

        if ( this->m_IsOutputDWI )
          {
          tensor.GetDWISamples(tempPixelDWI, *qSpaceOrientationMatrix, *bVector);
          // utlPrintVar1(true, tempPixel);
          // utlPrintVar1(true, tempPixel*scale*partialVolumeWeight);
          pixelDWI += tempPixelDWI*scale*partialVolumeWeight;
          }
        if ( this->m_IsOutputODF )
          {
          tensor.GetODFSamples(tempPixelODF, *rSpaceOrientationMatrix,this->m_ODFOrder, false);
          pixelODF += tempPixelODF*partialVolumeWeight;
          }
        if ( this->m_IsOutputEAP)
          {
          tensor.GetEAPSamples(tempPixelEAP, *rSpaceOrientationMatrix, *rVector,tau);
          pixelEAP += tempPixelEAP*partialVolumeWeight;
          // std::cout << "\n\nindex = " << index << std::endl << std::flush;
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
      else if (isCylinderModel)
        {
        cylinderAxis[0] = principalDirection[0];
        cylinderAxis[1] = principalDirection[1];
        cylinderAxis[2] = principalDirection[2];
        this->m_CylinderModel->SetCylinderAxis(cylinderAxis);
        if (this->m_IsOutputDWI)
          {
          if (this->GetDebug())
            this->m_CylinderModel->Print(std::cout<<"this->m_CylinderModel=");
          this->m_CylinderModel->ComputeDWISamples();
          tempPixelDWIVec = this->m_CylinderModel->GetDWISamples();
          for ( int i = 0; i < tempPixelDWIVec->size(); i += 1 ) 
            pixelDWI[i] += (*tempPixelDWIVec)[i]*scale*partialVolumeWeight;
          }
        if ( this->m_IsOutputODF )
          {
          utlGlobalException(true, "TODO: to be implemented");
          }
        if ( this->m_IsOutputEAP )
          {
          utlGlobalException(true, "TODO: to be implemented");
          }
        if ( this->m_IsOutputRTO )
          {
          utlGlobalException(true, "TODO: to be implemented");
          }
        if ( this->m_IsOutputMSD )
          {
          utlGlobalException(true, "TODO: to be implemented");
          }
        }

      }
    if (isOutputPeak)
      {
      if (this->m_PeakType==NXYZ || this->m_PeakType==NXYZV)
        pixelPeak[0] = numberOfPeaks;
      }

    utlException(sumPartialVolumeWeight<1e-9, "wrong partialVolumeWeight");
    if ( this->m_IsOutputDWI )
      pixelDWI /= sumPartialVolumeWeight;
    if ( this->m_IsOutputODF )
      pixelODF /= sumPartialVolumeWeight;
    if ( this->m_IsOutputEAP )
      pixelEAP /= sumPartialVolumeWeight;
    if ( this->m_IsOutputRTO )
      pixelRTO /= sumPartialVolumeWeight;
    if ( this->m_IsOutputMSD )
      pixelMSD /= sumPartialVolumeWeight;
    if (this->m_IsOutputODF && this->m_ODFOrder!=2)
      {
      PrecisionType normFactor = 4*vnl_math::pi*utl::GetSumOfVector<OutputImagePixelType>(pixelODF,pixelODF.Size()) / pixelODF.Size();
      if (normFactor!=0)
        pixelODF /= normFactor;
      }

    // Add Rician noise
    if ( this->m_IsOutputDWI && (sigma>0 || snr>0)  && pixelDWI.GetSquaredNorm()>0 )
      {
      PrecisionType sigma_real = -1;
      if ( snr > 0 && sigma <= 0)
        {
        double mean = 0;
        for ( unsigned int k=0; k<numberOfComponentsPerPixel4DWI; k++ )
          mean += pixelDWI[k];
        mean /= numberOfComponentsPerPixel4DWI;
        sigma_real = mean/snr;
        }
      else if ( snr <= 0 && sigma > 0)
//        sigma_real = scale * sigma;
        sigma_real = sigma;
      pixelDWI = utl::AddNoise<OutputImagePixelType>(pixelDWI, pixelDWI.GetSize(), sigma_real, true);
      }
    
    if ( isOutputPeak )
      {
      outputPeak->SetPixel(index, pixelPeak);
      }
    if ( this->m_IsOutputDWI )
      {
      outputDWI->SetPixel(index, pixelDWI); 
      b0Image->SetPixel(index, scale);
      }
    if ( this->m_IsOutputODF )
      outputODF->SetPixel(index, pixelODF); 
    if ( this->m_IsOutputEAP )
      outputEAP->SetPixel(index, pixelEAP); 
    if ( this->m_IsOutputRTO )
      outputRTO->SetPixel(index, pixelRTO); 
    if ( this->m_IsOutputMSD )
      outputMSD->SetPixel(index, pixelMSD); 
    }
    
  ImageRegionIteratorWithIndex<OutputImageType> outputDWIIt, outputODFIt, outputEAPIt, outputPeakIt;
  ImageRegionIteratorWithIndex<ScalarImageType> b0It, outputRTOIt, outputMSDIt;
  if ( this->m_IsOutputDWI )
    {
    outputDWIIt = ImageRegionIteratorWithIndex<OutputImageType>( outputDWI, outputDWI->GetRequestedRegion() );
    b0It = ImageRegionIteratorWithIndex<ScalarImageType>( b0Image, outputDWI->GetRequestedRegion() );
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
        
      
  if (m_BackgroundDiffusionParameterValues.size()>0)
    backgroundDiffusionParameterValues = m_BackgroundDiffusionParameterValues;

  // Background      
  if ( backgroundDiffusionParameterValues.size() > 0 )
  {
    numberOfParameterSets = backgroundDiffusionParameterValues.size() /numberOfComponentsPerParameterSet;
    // utl::PrintVector(backgroundDiffusionParameterValues, "backgroundDiffusionParameterValues");
    // utlPrintVar3(true, numberOfComponentsPerParameterSet, numberOfParameterSets, numberOfOrientationComponentsPerParameterSet);
    if (this->m_IsOutputDWI)
      pixelDWI.Fill( 0 );
    if (this->m_IsOutputEAP)
      pixelEAP.Fill( 0 );
    if (this->m_IsOutputODF)
      pixelODF.Fill( 0 );
    pixelRTO=0, pixelMSD=0;
    PrecisionType sumPartialVolumeWeight=0.0;

    for ( unsigned int s=0; s<numberOfParameterSets; s++ )
      {
      for ( unsigned int d=0; d<numberOfOrientationComponentsPerParameterSet; d++ )
        principalDirection[d] = backgroundDiffusionParameterValues[s* numberOfComponentsPerParameterSet + d];

      // std::cout << "principalDirection = " << principalDirection << std::endl << std::flush;
      // utlException (std::abs(norm)<1e-8, "directions have zero norm");
      norm = principalDirection.GetNorm();
      if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS || this->m_ModelType==Superclass::CYLINDER_SPHERICAL_MODEL)
        {
        principalDirection *= vnl_math::pi/180.0;
        double theta=principalDirection[0], phi=principalDirection[1];
        utl::spherical2Cartesian(1.0,theta,phi,principalDirection[0],principalDirection[1],principalDirection[2]);
        }
      else
        {
        for ( unsigned int d=0; d<3; d++ )
          principalDirection[d] /= norm;
        }

      for ( unsigned int d=0; d<numberOfDiffusivityComponentsPerParameterSet; d++ )
        {
        diffusivities[d] = backgroundDiffusionParameterValues[s* numberOfComponentsPerParameterSet 
            + numberOfOrientationComponentsPerParameterSet + d ];
        }
      if (diffusivities.size()>0)
        std::sort(diffusivities.begin(), diffusivities.end(), std::greater<double>());
      
      partialVolumeWeight = backgroundDiffusionParameterValues[s* numberOfComponentsPerParameterSet 
          + numberOfOrientationComponentsPerParameterSet 
          + numberOfDiffusivityComponentsPerParameterSet];

      sumPartialVolumeWeight += partialVolumeWeight;

      if (isTensorModel)
        {
        std::vector<double> vec(3,1);
        if (this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS || this->m_ModelType==Superclass::SYMMETRICAL_TENSOR_IN_SPHERICAL_COORDS)
          vec[0]=diffusivities[0], vec[1]=diffusivities[1], vec[2]=diffusivities[1];
        else if (this->m_ModelType==Superclass::TENSOR_IN_EULER_ANGLES)
          vec[0]=diffusivities[0], vec[1]=diffusivities[1], vec[2]=diffusivities[2];
        tensor.Fill(0.0);
        tensor.SetEigenValues(vec);
        // tensor.Print(std::cout<<"tensor = ");
        Matrix<double,3,3> rotation;
        utl::RotationMatrixFromVectors<Vector<double,3>, Matrix<double,3> >(e1, principalDirection, rotation); 
        tensor.Rotate(rotation);
        // tensor.Print(std::cout<<"tensor = ");

        if ( this->m_IsOutputDWI )
          {
          tensor.GetDWISamples(tempPixelDWI, *qSpaceOrientationMatrix, *bVector);
          pixelDWI += tempPixelDWI*backgroundScale*partialVolumeWeight;
          }
        if ( this->m_IsOutputODF )
          {
          tensor.GetODFSamples(tempPixelODF, *rSpaceOrientationMatrix,this->m_ODFOrder, false);
          pixelODF += tempPixelODF*partialVolumeWeight;
          }
        if ( this->m_IsOutputEAP)
          {
          tensor.GetEAPSamples(tempPixelEAP, *rSpaceOrientationMatrix, *rVector,tau);
          pixelEAP += tempPixelEAP*partialVolumeWeight;
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
      else if (isCylinderModel)
        {
        cylinderAxis[0] = principalDirection[0];
        cylinderAxis[1] = principalDirection[1];
        cylinderAxis[2] = principalDirection[2];
        this->m_CylinderModel->SetCylinderAxis(cylinderAxis);
        if (this->m_IsOutputDWI)
          {
          // cyl->Print(std::cout<<"cy1=");
          this->m_CylinderModel->ComputeDWISamples();
          tempPixelDWIVec = this->m_CylinderModel->GetDWISamples();
          for ( int i = 0; i < tempPixelDWIVec->size(); i += 1 ) 
            pixelDWI[i] += (*tempPixelDWIVec)[i]*backgroundScale*partialVolumeWeight;
          }

        }

      }
    utlException(sumPartialVolumeWeight<1e-9, "wrong partialVolumeWeight");
    if ( this->m_IsOutputDWI )
      pixelDWI /= sumPartialVolumeWeight;
    if ( this->m_IsOutputODF )
      pixelODF /= sumPartialVolumeWeight;
    if ( this->m_IsOutputEAP )
      pixelEAP /= sumPartialVolumeWeight;
    if ( this->m_IsOutputRTO )
      pixelRTO /= sumPartialVolumeWeight;
    if ( this->m_IsOutputMSD )
      pixelMSD /= sumPartialVolumeWeight;
    if (this->m_IsOutputODF && this->m_ODFOrder!=2)
      {
      PrecisionType normFactor = 4*vnl_math::pi*utl::GetSumOfVector<OutputImagePixelType>(pixelODF,pixelODF.Size()) / pixelODF.Size();
      if (normFactor!=0)
        pixelODF /= normFactor;
      }
  
    // Add Rician noise
    PrecisionType sigma_real = -1;
    if ( this->m_IsOutputDWI && (sigma>0 || snr>0)  && pixelDWI.GetSquaredNorm()>0 )
      {
      if ( snr > 0 && sigma <= 0)
        {
        double mean = 0;
        for ( unsigned int k=0; k<numberOfComponentsPerPixel4DWI; k++ )
          mean += pixelDWI[k];
        mean /= numberOfComponentsPerPixel4DWI;
        sigma_real = mean/snr;
        }
      else if ( snr <= 0 && sigma > 0)
//        sigma_real = backgroundScale * sigma;
        sigma_real = sigma;
      }

    // Set pixels
    if (this->m_IsOutputDWI)
      {
      b0It.GoToBegin();
      outputDWIIt.GoToBegin();
      OutputImagePixelType pixel_noise=pixelDWI;
      while ( !b0It.IsAtEnd() )
        {
        if ( outputDWIIt.Get().GetNorm() == 0 )
          {
          if (sigma_real>0)
            pixel_noise = utl::AddNoise<OutputImagePixelType>(pixelDWI, pixelDWI.GetSize(), sigma_real, true);
          outputDWIIt.Set( pixel_noise );
          b0It.Set( backgroundScale );
          }
        ++outputDWIIt;
        ++b0It;
        }
      }

    if (this->m_IsOutputODF)
      {
      outputODFIt.GoToBegin();
      while ( !outputODFIt.IsAtEnd() )
        {
        if ( outputODFIt.Get().GetNorm() == 0 )
          outputODFIt.Set( pixelODF );
        ++outputODFIt;
        }
      }

    if (this->m_IsOutputEAP)
      {
      outputEAPIt.GoToBegin();
      while ( !outputEAPIt.IsAtEnd() )
        {
        if ( outputEAPIt.Get().GetNorm() == 0 )
          outputEAPIt.Set( pixelEAP );
        ++outputEAPIt;
        }
      }
    
    if (this->m_IsOutputRTO)
      {
      outputRTOIt.GoToBegin();
      while ( !outputRTOIt.IsAtEnd() )
        {
        if ( outputRTOIt.Get() == 0 )
          outputRTOIt.Set( pixelRTO );
        ++outputRTOIt;
        }
      }
    if (this->m_IsOutputMSD)
      {
      outputMSDIt.GoToBegin();
      while ( !outputMSDIt.IsAtEnd() )
        {
        if ( outputMSDIt.Get() == 0 )
          outputMSDIt.Set( pixelMSD );
        ++outputMSDIt;
        }
      }
  }

}


} //namespace ITK

#endif
