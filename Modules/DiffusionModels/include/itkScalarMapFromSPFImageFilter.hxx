/**
 *       @file  itkScalarMapFromSPFImageFilter.hxx
 *      @brief  
 *     Created  "03-18-2014
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkScalarMapFromSPFImageFilter_hxx
#define __itkScalarMapFromSPFImageFilter_hxx

#include <itkProgressReporter.h>
#include "itkScalarMapFromSPFImageFilter.h"
#include "utl.h"
#include "itkSphericalPolarFourierGenerator.h"
#include "itkSpecialFunctionGenerator.h"

namespace itk
{  

template< class TInputImage, class TOutputImage >
void
ScalarMapFromSPFImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();
  typename TInputImage::ConstPointer inputPtr = this->GetInput();
  typename TOutputImage::Pointer outputPtr = this->GetOutput();
  itk::CopyImageInformation(inputPtr, outputPtr);
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
ScalarMapFromSPFImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_SumWeight = m_SumWeight;
  rval->m_MapType = m_MapType;
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
ScalarMapFromSPFImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData ( )
{
  // utlShowPosition(true);
  m_SumWeight.ReSize(this->m_RadialRank+1);
  if (m_MapType==RTO)
    {
    for ( int i = 0; i <= this->m_RadialRank ; ++i ) 
      {
      double sign = utl::IsEven(i) ? 1.0 : -1.0;
      m_SumWeight[i] = sign * std::sqrt( utl::GammaHalfInteger(i+1.5) / utl::Factorial(i) );
      }
    m_SumWeight %= 4*std::sqrt(M_PI);
    }
  else if (m_MapType==MSD)
    {
    for ( int n = 0; n <= this->m_RadialRank; ++n  ) 
      {
      if (n==0)
        m_SumWeight[n] = utl::Lagurre(n,0.5,0.0);
      else
        m_SumWeight[n] = 2*utl::Lagurre(n-1,1.5,0.0) + utl::Lagurre(n,0.5,0.0);
      }
    m_SumWeight.Scale(3.0/(8.0*M_PI*M_PI*std::sqrt(M_PI)));
    }
  else if (m_MapType==PFA)
    {

    }

  Superclass::BeforeThreadedGenerateData();
}

template <class TInputImage, class TOutputImage>
void
ScalarMapFromSPFImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData( const typename TOutputImage::RegionType &outputRegionForThread, ThreadIdType threadId)
{
  // utlShowPosition(true);
  typename TInputImage::ConstPointer  inputPtr = this->GetInput();
  typename TOutputImage::Pointer outputPtr = this->GetOutput();

  // Define the iterators
  ImageRegionConstIteratorWithIndex<TInputImage>  inputIt(inputPtr, outputRegionForThread);
  ImageRegionIteratorWithIndex<TOutputImage> outputIt(outputPtr, outputRegionForThread);
  ImageRegionIteratorWithIndex<MaskImageType> maskIt;
  if (this->IsMaskUsed())
    maskIt = ImageRegionIteratorWithIndex<MaskImageType>(this->m_MaskImage, outputRegionForThread);

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  typename TInputImage::PixelType inputPixel;
  typename TInputImage::IndexType inputIndex;
  typename TOutputImage::PixelType outputPixel;
  ScalarImagePointer scaleImage = this->GetScaleImage();
  
  unsigned int inputDim = inputPtr->GetNumberOfComponentsPerPixel();
  inputPixel.SetSize(inputDim);

  inputIt.GoToBegin();
  outputIt.GoToBegin();

  double scale = this->m_BasisScale;
  int dimSH = utl::RankToDimSH(this->m_SHRank);

  while( !inputIt.IsAtEnd() ) 
    {
    if (!this->IsMaskUsed() || (this->IsMaskUsed() && maskIt.Get()>0))
      {
      inputPixel = inputIt.Get();
      if (inputPixel.GetSquaredNorm()>1e-8)
        {
        inputIndex = inputIt.GetIndex();
        if (this->GetDebug())
          std::cout << "index = " << inputIndex << std::endl << std::flush;

        if (!IsImageEmpty(scaleImage))
          scale = scaleImage->GetPixel(inputIndex);

        if (m_MapType==RTO)
          {
          outputPixel=0;
          for ( int i = 0; i <= this->m_RadialRank; ++i ) 
            outputPixel += inputPixel[i*dimSH]*m_SumWeight[i];
          outputPixel *= std::pow(scale, 0.75);
          }
        else if (m_MapType==MSD)
          {
          VectorType scaleWeight(this->m_RadialRank+1);
          for ( int n = 0; n <= this->m_RadialRank; n += 1 ) 
            {
            scaleWeight[n] = 2.0 / std::pow(scale,(double)1.5) * utl::Factorial(n) / utl::Gamma(n+1.5);
            scaleWeight[n] = std::sqrt(scaleWeight[n])/scale;
            }
          outputPixel=0;
          for ( int i = 0; i <= this->m_RadialRank; ++i ) 
            outputPixel += inputPixel[i*dimSH]*m_SumWeight[i]*scaleWeight[i];
          }
        else if (m_MapType==PFA)
          {
          double sum = 0;
          for ( int i = 0; i <= this->m_RadialRank; ++i ) 
            sum += inputPixel[i*dimSH]*inputPixel[i*dimSH];
          outputPixel = std::sqrt(1- sum/inputPixel.GetSquaredNorm() );
          }

        }
      else
        outputPixel=0;
      }
    else
      outputPixel=0;

    outputIt.Set(outputPixel);
    progress.CompletedPixel();    
    if (this->IsMaskUsed())
      ++maskIt;
    ++inputIt;
    ++outputIt;
    progress.CompletedPixel();  // potential exception thrown here
    }
}

template< class TInputImage, class TOutputImage >
void 
ScalarMapFromSPFImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  utl::PrintUtlVector(m_SumWeight, "m_SumWeight", " ", os<<indent);
}

}

#endif 



