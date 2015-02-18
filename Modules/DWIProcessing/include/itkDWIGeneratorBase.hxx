/**
 *       @file  itkDWIGeneratorBase.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-22-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkDWIGeneratorBase_hxx
#define __itkDWIGeneratorBase_hxx

#include "itkDWIGeneratorBase.h"
#include "utl.h"

namespace itk
{

template <class TOutputImage, class TScalarImage>
DWIGeneratorBase<TOutputImage, TScalarImage>
::DWIGeneratorBase() 
{
  this->SetNumberOfRequiredOutputs(1); // at least 1 output
  this->SetNthOutput( 0, ( TOutputImage::New() ).GetPointer() ); // dwi
  this->SetNthOutput( 1, ( TScalarImage::New() ).GetPointer() ); // b0
  this->SetNthOutput( 2, ( TOutputImage::New() ).GetPointer() ); // ODF
  this->SetNthOutput( 3, ( TOutputImage::New() ).GetPointer() ); // EAP
  this->SetNthOutput( 4, ( TOutputImage::New() ).GetPointer() ); // peak
  this->SetNthOutput( 5, ( TScalarImage::New() ).GetPointer() ); // RTO
  this->SetNthOutput( 6, ( TScalarImage::New() ).GetPointer() ); // MSD
  
  m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
  m_SamplingSchemeRSpace = SamplingSchemeQSpaceType::New();

  m_CylinderModel = CylinderModelType::New();

  m_NoiseSigma = -1; // Unset
  m_SNR = -1; // Unset
  m_ODFOrder = 2;
  m_ModelType = SYMMETRICAL_TENSOR_IN_CARTESIAN_COORDS;
  m_B0Scale = 1.0;  // if it is not set, 1 is used
  m_MaxNumberOfPeaks = -1;
  m_PeakType = NXYZ;
  
  m_OutputOrigin.Fill(0);
  m_OutputDirection.SetIdentity();
  for ( int i = 0; i < OutputImageType::ImageDimension; i += 1 ) 
    {
    m_OutputSize[i]=1.0;
    m_OutputSpacing[i]=1.0;
    }

  m_IsOutputDWI=false;
  m_IsOutputEAP=false;
  m_IsOutputODF=false;
  m_IsOutputRTO=false;
  m_IsOutputMSD=false;
}

template <class TOutputImage, class TScalarImage>
DWIGeneratorBase<TOutputImage, TScalarImage>
::~DWIGeneratorBase()
{
}

template <class TOutputImage, class TScalarImage>
typename LightObject::Pointer
DWIGeneratorBase<TOutputImage, TScalarImage>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_NoiseSigma = m_NoiseSigma;
  rval->m_SNR = m_SNR;

  rval->m_ODFOrder = m_ODFOrder;

  rval->m_ModelType = m_ModelType;
  rval->m_B0Scale = m_B0Scale;

  rval->m_MaxNumberOfPeaks = m_MaxNumberOfPeaks;
  rval->m_PeakType = m_PeakType;
  
  rval->m_CylinderModel = m_CylinderModel;

  rval->m_OutputSize = m_OutputSize;
  rval->m_OutputSpacing = m_OutputSpacing;
  rval->m_OutputOrigin = m_OutputOrigin;
  rval->m_OutputDirection = m_OutputDirection;

  rval->m_SamplingSchemeQSpace = m_SamplingSchemeQSpace;
  rval->m_SamplingSchemeRSpace = m_SamplingSchemeRSpace;
  
  rval->m_IsOutputDWI=m_IsOutputDWI;
  rval->m_IsOutputEAP=m_IsOutputEAP;
  rval->m_IsOutputODF=m_IsOutputODF;
  rval->m_IsOutputRTO=m_IsOutputRTO;
  rval->m_IsOutputMSD=m_IsOutputMSD;

  return loPtr;
}

template <class TOutputImage, class TScalarImage>
void 
DWIGeneratorBase<TOutputImage, TScalarImage>
::Initialization()
{
  utlShowPosition(this->GetDebug());
  utlGlobalException(!m_IsOutputDWI && !m_IsOutputODF && !m_IsOutputEAP && !m_IsOutputRTO && !m_IsOutputMSD && m_MaxNumberOfPeaks<=0, "no output!");
  utlGlobalException(m_IsOutputDWI && this->GetNumberOfQSpaceSamples()==0, "need to set sampling schemes in q-space for DWI");
  utlGlobalException((m_IsOutputODF || m_IsOutputEAP) && this->GetNumberOfRSpaceSamples()==0, "need to set sampling schemes in r-space for EAP or ODF");

  MatrixPointer qSpaceOrientationMatrix = m_SamplingSchemeQSpace->GetOrientationsCartesian();
  STDVectorPointer bVector = m_SamplingSchemeQSpace->GetBVector();
  MatrixPointer rSpaceOrientationMatrix = m_SamplingSchemeRSpace->GetOrientationsCartesian();
  STDVectorPointer rVector = m_SamplingSchemeRSpace->GetRadiusVector();
  utlGlobalException(m_IsOutputDWI && qSpaceOrientationMatrix->Rows()!=bVector->size(), "Need to set OrientationsCartesian and BVector. qSpaceOrientationMatrix->Rows()=" << qSpaceOrientationMatrix->Rows() << ", bVector->size()=" << bVector->size());
  utlGlobalException( (m_IsOutputODF || m_IsOutputEAP) && rSpaceOrientationMatrix->Rows()!=rVector->size(), "Need to set rSpaceOrientationMatrix and rVector");
  utlGlobalException(m_SNR>0 && m_NoiseSigma>0, "Only one of m_SNR and m_NoiseSigma is allowed.");
  
  double tau = this->m_SamplingSchemeQSpace->GetTau();
  utlGlobalException(tau!=this->m_SamplingSchemeRSpace->GetTau(), "need to use the same tau in m_SamplingSchemeQSpace and m_SamplingSchemeRSpace");
}



template <class TOutputImage, class TScalarImage>
void 
DWIGeneratorBase<TOutputImage, TScalarImage>
::AllocateOutputs()
{
  utlShowPosition(this->GetDebug());
  
  bool isOutputPeak = this->m_MaxNumberOfPeaks>0;
  
  OutputImagePointer outputDWI=this->GetDWIImage(), outputODF=this->GetODFImage(), outputEAP=this->GetEAPImage(), outputPeak=this->GetPeakImage();
  ScalarImagePointer b0Image = this->GetB0Image();
  ScalarImagePointer rtoImage = this->GetRTOImage();
  ScalarImagePointer msdImage = this->GetMSDImage();
  OutputImageIndexType outputStartIndex;
  OutputImageRegionType outputRegion;
  outputStartIndex.Fill(0);
  outputRegion.SetIndex( outputStartIndex );
  outputRegion.SetSize( m_OutputSize );  
  
  int numberOfComponentsPerPixel4DWI = this->GetNumberOfQSpaceSamples();
  int numberOfComponentsPerPixel4EAP = this->GetNumberOfRSpaceSamples();
  int numberOfComponentsPerPixel4ODF = this->GetNumberOfRSpaceSamples();
  int numberOfComponentsPerPixel4Peak = isOutputPeak?PeakContainerHelper::GetDimension(m_PeakType,m_MaxNumberOfPeaks):-1;
  OutputImagePixelType zeroPixel;
  
  if (m_IsOutputDWI)
    {
    // outputDWI = OutputImageType::New();
    outputDWI->SetRegions( outputRegion );
    outputDWI->SetSpacing( m_OutputSpacing );
    outputDWI->SetOrigin( m_OutputOrigin );
    outputDWI->SetDirection( m_OutputDirection );
    outputDWI->SetNumberOfComponentsPerPixel( numberOfComponentsPerPixel4DWI );
    zeroPixel.SetSize( numberOfComponentsPerPixel4DWI );
    zeroPixel.Fill(0);
    outputDWI->Allocate();
    outputDWI->FillBuffer( zeroPixel );
  
    // b0Image = ScalarImageType::New();
    // b0Image = this->GetB0Image();
    b0Image->SetRegions( outputRegion );
    b0Image->SetSpacing( m_OutputSpacing );
    b0Image->SetOrigin( m_OutputOrigin );
    b0Image->SetDirection( m_OutputDirection );
    b0Image->Allocate();
    b0Image->FillBuffer( 0 );
    }
  if (m_IsOutputEAP)
    {
    // outputEAP = OutputImageType::New();
    outputEAP->SetRegions( outputRegion );
    outputEAP->SetSpacing( m_OutputSpacing );
    outputEAP->SetOrigin( m_OutputOrigin );
    outputEAP->SetDirection( m_OutputDirection );
    outputEAP->SetNumberOfComponentsPerPixel( numberOfComponentsPerPixel4EAP );
    zeroPixel.SetSize( numberOfComponentsPerPixel4EAP );
    zeroPixel.Fill(0);
    outputEAP->Allocate();
    outputEAP->FillBuffer( zeroPixel );
    }
  if (m_IsOutputODF)
    {
    // outputODF = OutputImageType::New();
    outputODF->SetRegions( outputRegion );
    outputODF->SetSpacing( m_OutputSpacing );
    outputODF->SetOrigin( m_OutputOrigin );
    outputODF->SetDirection( m_OutputDirection );
    outputODF->SetNumberOfComponentsPerPixel( numberOfComponentsPerPixel4ODF );
    zeroPixel.SetSize( numberOfComponentsPerPixel4ODF );
    zeroPixel.Fill(0);
    outputODF->Allocate();
    outputODF->FillBuffer( zeroPixel );
    }
  if (isOutputPeak)
    {
    // outputODF = OutputImageType::New();
    outputPeak->SetRegions( outputRegion );
    outputPeak->SetSpacing( m_OutputSpacing );
    outputPeak->SetOrigin( m_OutputOrigin );
    outputPeak->SetDirection( m_OutputDirection );
    outputPeak->SetNumberOfComponentsPerPixel( numberOfComponentsPerPixel4Peak );
    zeroPixel.SetSize( numberOfComponentsPerPixel4Peak );
    zeroPixel.Fill(0);
    outputPeak->Allocate();
    outputPeak->FillBuffer( zeroPixel );
    }
  if (m_IsOutputRTO)
    {
    rtoImage->SetRegions( outputRegion );
    rtoImage->SetSpacing( m_OutputSpacing );
    rtoImage->SetOrigin( m_OutputOrigin );
    rtoImage->SetDirection( m_OutputDirection );
    rtoImage->Allocate();
    rtoImage->FillBuffer( 0 );
    }
  if (m_IsOutputMSD)
    {
    msdImage->SetRegions( outputRegion );
    msdImage->SetSpacing( m_OutputSpacing );
    msdImage->SetOrigin( m_OutputOrigin );
    msdImage->SetDirection( m_OutputDirection );
    msdImage->Allocate();
    msdImage->FillBuffer( 0 );
    }
}

template <class TOutputImage, class TScalarImage>
void 
DWIGeneratorBase<TOutputImage, TScalarImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Gneratortype: " << m_ModelType  << std::endl;
  os << indent << "ODFOrder: " << m_ODFOrder  << std::endl;
  os << indent << "NoiseSigma: " << m_NoiseSigma  << std::endl;
  os << indent << "SNR: " << m_SNR  << std::endl;
  if (m_MaxNumberOfPeaks>0)
    PrintVar2(true, m_MaxNumberOfPeaks, PeakContainerHelper::GetString(m_PeakType), os<<indent);
  PrintVar1(true, m_OutputSize, os<<indent);
  os << indent << "m_SamplingSchemeQSpace: " << m_SamplingSchemeQSpace << std::endl;
  os << indent << "m_SamplingSchemeRSpace: " << m_SamplingSchemeRSpace << std::endl;
}

}


#endif 


