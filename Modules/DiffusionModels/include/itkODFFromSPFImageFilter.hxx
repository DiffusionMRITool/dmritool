/**
 *       @file  itkODFFromSPFImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-08-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkODFFromSPFImageFilter_hxx
#define __itkODFFromSPFImageFilter_hxx

#include <itkProgressReporter.h>
#include "itkODFFromSPFImageFilter.h"
#include "utl.h"
#include "itkSphericalPolarFourierGenerator.h"
#include "itkSpecialFunctionGenerator.h"

namespace itk
{  

template< class TInputImage, class TOutputImage >
void
ODFFromSPFImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();
  typename TOutputImage::Pointer outputPtr = this->GetOutput();
  unsigned int numberOfComponentsPerPixel = utl::RankToDimSH(this->GetSHRank());  
  outputPtr->SetNumberOfComponentsPerPixel(numberOfComponentsPerPixel);
}

template< class TInputImage, class TOutputImage >
void
ODFFromSPFImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters()
{
  Superclass::VerifyInputParameters();
  utlGlobalException(this->m_ODFOrder<-1, "m_ODFOrder should be -1, 0, 1, ...");
  utlGlobalException(this->m_ODFOrder==-1 && this->m_BMax<=0,"for Funk-Radon transform, bvalue can not be negative!");
}

template< class TInputImage, class TOutputImage >
void
ODFFromSPFImageFilter< TInputImage, TOutputImage >
::ComputeSPFToFeatureTransform()
{
  utlException(this->m_SHRank<0 || this->m_RadialRank<0, "need to set this->m_SHRank and this->m_RadialRank");

  double qValue;
  if (this->m_BMax>0)
    qValue = std::sqrt(this->m_BMax/(4*M_PI*M_PI*this->m_Tau));
  else
    {
    utlException(this->m_ODFOrder==-1,"for Funk-Radon transform, this->m_BMax can not be negative");
    qValue = -1;
    }
  double C = qValue>0 ? qValue*qValue/this->m_BasisScale : -1;
  if ( C>0)
    {
    if (this->GetDebug())
      {
      std::cout << "integration in a disk" << std::endl;
      std::cout << "C = " << C << std::endl;
      std::cout << "qValue = " << qValue << ", m_BasisScale = " << this->m_BasisScale << std::endl;
      }
    }
  else
    {
    if (this->GetDebug())
      {
      std::cout << "integration in the whole plane" << std::endl;
      // C should be -1
      std::cout << "C = " << C << ", m_BasisScale = " << this->m_BasisScale << std::endl;
      }
    }

  if (this->m_BasisType==Superclass::SPF || this->m_BasisType==Superclass::DSPF)
    {
    typedef SphericalPolarFourierRadialGenerator<double> SPFGenerator;
    int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
    int n_b_ra = this->m_RadialRank+1;
    int n_b = n_b_sh * n_b_ra;
    this->m_SPFToFeatureTransform = MatrixPointer (new MatrixType(n_b_sh,n_b) );
    this->m_SPFToFeatureTransform->Fill(0.0);
    if (this->m_BasisType==Superclass::SPF)
      {
      for ( int n = 0; n <= this->m_RadialRank; n += 1 ) 
        {
        double temp;
        double first_temp = 2.0 / std::pow((double)this->m_BasisScale,(double)1.5) * utl::Factorial(n) / utl::Gamma(n+1.5);
        first_temp = std::sqrt(first_temp);
        switch ( this->m_ODFOrder )
          {
        case -1 :
            {
            // std::cout << "order 0, ODF by Tuch, Funk-Radon Transform" << std::endl;
            temp = first_temp * std::exp( -1.0*qValue*qValue / (2*this->m_BasisScale) ) * utl::Lagurre(n, 0.5, qValue*qValue/this->m_BasisScale);
            for ( int j = 0; j < n_b_sh; j += 1 ) 
              (*this->m_SPFToFeatureTransform)(j,n*n_b_sh+j) = m_P(j) * temp;
            break;
            }
        case 0 :
            {
            // std::cout << "order 0, ODF by Tuch, integrate in a disk or a plane" << std::endl;
            double second_temp = 0;
            if ( C>0 )
              {
              // formula (7)
              for ( int i = 0; i <= n; i += 1 ) 
                second_temp += utl::Binomial(n+0.5,n-i) * std::pow((double)-2.0,(double)i) * utl::GammaLower(i+1,0.5*C) / utl::Factorial(i);
              }
            else
              {
              // formula (6)
              for ( int i = 0; i <= n; i += 1 ) 
                {
                double tmp = utl::Binomial(i-0.5,i);
                second_temp += (n-i)%2==0 ? tmp : -tmp;
                }
              }
            second_temp *= this->m_BasisScale;
            temp = first_temp*second_temp;
            //**************************************************************

            for ( int j = 0; j < n_b_sh; j += 1 ) 
              (*this->m_SPFToFeatureTransform)(j,n*n_b_sh+j) = m_P(j) * temp; 
            break;
            }
        case 2 :
            {
            // std::cout << "order 2, OPDF (maginal pdf), integrate in a disk or a plane" << std::endl;
            // NOTE: -1 is in transformBasis
            double second_temp = 0;
            if ( C>0 )
              {
              // formula (15)
              for ( int i = 1; i <= n; i += 1 ) 
                {
                double tmp = utl::Binomial(n+0.5,n-i) * std::pow((double)2.0,(double)i) * utl::GammaLower(i,0.5*C) / utl::Factorial(i);
                second_temp += (i%2==0 ? tmp : -tmp);  
                }
              }
            else
              {
              // formula (16)
              for ( int i = 1; i <= n; i += 1 ) 
                {
                double tmp = utl::Binomial(n+0.5,n-i) * std::pow((double)2.0,(double)i) / i;
                second_temp += (i%2==0 ? tmp : -tmp);  
                }
              }
            second_temp *= 0.5;
            temp = first_temp*second_temp;
            //**************************************************************

            for ( int j = 0; j < n_b_sh; j += 1 ) 
              (*this->m_SPFToFeatureTransform)(j,n*n_b_sh+j) = -1.0 * m_P(j) * m_L(j) * temp;
            break;
            }
        default :
          utlGlobalException(true,"wrong type");
          break;
          }
        }

      }
    else
      {
      utlGlobalException(true, "TODO");
      }
    }
  if (this->GetDebug())
    utl::PrintUtlMatrix(*this->m_SPFToFeatureTransform, "m_SPFToFeatureTransform");
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
ODFFromSPFImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_P = m_P;
  rval->m_L = m_L;
  rval->m_ODFOrder = m_ODFOrder;
  rval->m_BMax = m_BMax;
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
ODFFromSPFImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData ( )
{
  // utlShowPosition(true);
  this->SetSPFIEstimator();

  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  m_P.ReSize(n_b_sh);
  m_L.ReSize(n_b_sh);
  int j = 0;
  double diag;
  for(int l = 0; l <= this->m_SHRank; l+=2)
    {
    diag = 2*M_PI* utl::LegendrePolynomialAt0(l);
    for(int m = -l; m <= l; m++) 
      {
      m_L(j) = -l*(l+1);
      m_P(j) = diag;
      j++;
      }
    }
  // utl::PrintUtlVector(m_L,"m_L");
  // utl::PrintUtlVector(m_P,"m_P");
  Superclass::BeforeThreadedGenerateData();
}

template <class TInputImage, class TOutputImage>
void
ODFFromSPFImageFilter<TInputImage,TOutputImage>
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
  
  unsigned int outputDim = outputPtr->GetNumberOfComponentsPerPixel();;
  outputPixel.SetSize(outputDim);
  unsigned int inputDim = inputPtr->GetNumberOfComponentsPerPixel();
  inputPixel.SetSize(inputDim);

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  VectorType spfVec, result;
  
  Pointer selfClone = this->Clone();

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

        if (!IsImageEmpty(this->m_ScaleImage))
          selfClone->SetBasisScale(this->m_ScaleImage->GetPixel(inputIndex));
        if (selfClone->m_SPFToFeatureTransform->Size()==0)
          selfClone->ComputeSPFToFeatureTransform();

        spfVec = utl::VariableLengthVectorToUtlVector(inputPixel);
        switch ( this->m_ODFOrder )
          {
        case -1 :
            {
            // if (this->GetDebug())
            //   std::cout << "order 0, ODF by Tuch, Funk-Radon Transform" << std::endl;
            utlException(this->m_BMax<=0,"for Funk-Radon transform, bvalue can not be negative");
            result = (*selfClone->m_SPFToFeatureTransform) * spfVec;
            double normalizefactor = 1.0 / (std::sqrt(4*M_PI)*result[0]);
            result %= normalizefactor;
            break;
            }
        case 0 :
            {
            // if (this->GetDebug())
            //   std::cout << "order 0, ODF by Tuch, integrate in a disk or a plane" << std::endl;
            result = (*selfClone->m_SPFToFeatureTransform) * spfVec;
            double normalizefactor = 1.0 / (std::sqrt(4*M_PI)*result[0]);
            result %= normalizefactor;
            break;
            }
        case 2 :
            {
            // if (this->GetDebug())
            //   std::cout << "order 2, OPDF (maginal pdf), integrate in a disk or a plane" << std::endl;
            result = (*selfClone->m_SPFToFeatureTransform) * spfVec / (8*M_PI*M_PI);
            result[0] += 1/std::sqrt(4*M_PI);
            break;
            }
        default :
          utlGlobalException(true,"Wrong type");
          break;
          }
        if (this->GetDebug())
          utl::PrintUtlVector(result, "odf");

        outputPixel = utl::UtlVectorToVariableLengthVector(result);

        }
      else
        outputPixel.Fill(0.0);
      }
    else
      outputPixel.Fill(0.0);

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
ODFFromSPFImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent
    << ", m_ODFOrder = " << this->m_ODFOrder
    << ", m_BMax = " << this->m_BMax
    << std::endl;
}

}

#endif 



