/**
 *       @file  itkProfileFromSPFImageFilter.hxx
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


#ifndef __itkProfileFromSPFImageFilter_hxx
#define __itkProfileFromSPFImageFilter_hxx

#include <itkProgressReporter.h>
#include "itkProfileFromSPFImageFilter.h"
#include "itkSphericalPolarFourierGenerator.h"
#include "utl.h"
#include "itkSphericalPolarFourierImageFilter.h"

namespace itk
{  

template< class TInputImage, class TOutputImage >
void
ProfileFromSPFImageFilter< TInputImage, TOutputImage >
::GenerateOutputInformation()
{
  utlShowPosition(this->GetDebug());
  Superclass::GenerateOutputInformation();
  typename TOutputImage::Pointer outputPtr = this->GetOutput();
  unsigned int numberOfComponentsPerPixel = 0;
  if (this->m_Orientations->Rows()==0)
    {
    // output SH coefficients
    numberOfComponentsPerPixel = utl::RankToDimSH(this->GetSHRank());  
    }
  else
    {
    // output samples
    numberOfComponentsPerPixel = this->m_Orientations->Rows();
    utlGlobalException(m_RadiusVector->size()>0 && m_RadiusVector->size()!=this->m_Orientations->Rows(), "wrong size of m_RadiusVector or m_Orientations");
    }
  outputPtr->SetNumberOfComponentsPerPixel(numberOfComponentsPerPixel);
}

template< class TInputImage, class TOutputImage >
void
ProfileFromSPFImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters() const
{
  utlShowPosition(this->GetDebug());
  Superclass::VerifyInputParameters();
  utlGlobalException(this->m_Radius<0 && this->m_RadiusVector->size()==0, "need to set radius or radiusVector");
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
ProfileFromSPFImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  rval->m_Radius = m_Radius;
  rval->m_RadiusVector = m_RadiusVector;
  return loPtr;
}

template< class TInputImage, class TOutputImage >
void
ProfileFromSPFImageFilter< TInputImage, TOutputImage >
::ComputeSPFToFeatureTransform()
{
  utlShowPosition(this->GetDebug());
  if (this->m_Orientations->Rows()==0)
    {
    // output SH coefficients
    utlException(this->m_Radius<0, "need to set radius");
    utlException(this->m_SHRank<0 || this->m_RadialRank<0, "need to set this->m_SHRank and this->m_RadialRank");

    // if it is in q-space, then convert b value to q value
    double radius = this->m_Radius;
    if ((this->m_IsInQSpace && !this->m_IsFourier) || (!this->m_IsInQSpace && this->m_IsFourier))
      radius = std::sqrt(this->m_Radius/(4*M_PI*M_PI*this->m_Tau));

    // std::cout << "radius = " << radius << std::endl << std::flush;
    if (this->m_BasisType==Superclass::SPF || this->m_BasisType==Superclass::DSPF)
      {
      typedef SphericalPolarFourierRadialGenerator<double> SPFGenerator;
      int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
      int n_b_ra = this->m_RadialRank+1;
      int n_b = n_b_sh * n_b_ra;
      typename SPFGenerator::Pointer spf = SPFGenerator::New();

      // NOTE: if use ReSize, different threads modify the same data block, which is not thread safe. 
      // this->m_SPFToFeatureTransform->ReSize(n_b_sh,n_b);
      this->m_SPFToFeatureTransform = MatrixPointer (new MatrixType(n_b_sh,n_b) );
      this->m_SPFToFeatureTransform->Fill(0.0);

      if ( (this->m_BasisType==Superclass::SPF && !this->m_IsFourier) || (this->m_BasisType==Superclass::DSPF && this->m_IsFourier) )
        {
        spf->SetSPFType(SPFGenerator::SPF);
        spf->SetScale(this->m_BasisScale);
        for ( int n = 0; n <= this->m_RadialRank; n += 1 ) 
          {
          spf->SetN(n);
          double spfValue = spf->Evaluate(radius);
          int jj=0; 
          for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
            {
            for ( int m = -l; m <= l; m += 1 ) 
              {
              (*this->m_SPFToFeatureTransform)(jj,n*n_b_sh+jj) = spfValue; 
              // utlPrintVar4(true, n,l,m, spf->Evaluate(radius));
              jj++;
              }
            }
          }
        }
      else
        {
        spf->SetSPFType(SPFGenerator::DSPF);
        spf->SetScale(this->m_BasisScale);
        for ( int n = 0; n <= this->m_RadialRank; n += 1 ) 
          {
          spf->SetN(n);
          int jj = 0;
          for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
            {
            spf->SetL(l);
            double spfValue = spf->Evaluate(radius);
            for ( int m = -l; m <= l; m += 1 ) 
              {
              (*this->m_SPFToFeatureTransform)(jj,n*n_b_sh+jj) = spfValue; 
              jj++;
              }
            }
          }
        }
      }
    }
  else
    {
    // output samples
    if (m_RadiusVector->size()==0)
      {
      utlException(this->m_Radius<0, "need to set m_Radius or m_RadiusVector");
      utl::MatchBVectorAndGradientMatrix(this->m_Radius, *m_RadiusVector, *this->m_Orientations);
      }
    utlException(m_RadiusVector->size()>0 && m_RadiusVector->size()!=this->m_Orientations->Rows(), "wrong size of m_RadiusVector or m_Orientations");
    this->m_SPFIEstimator->SetBasisScale(this->m_BasisScale);
    // std::cout << "this->m_BasisScale = " << this->m_BasisScale << std::endl << std::flush;
    this->m_SPFIEstimator->ComputeBasisMatrix();
    this->m_SPFToFeatureTransform = this->m_SPFIEstimator->GetBasisMatrix();
    // this->m_SPFIEstimator->Print(std::cout<<"m_SPFIEstimator : ");
    }

  if (this->GetDebug())
    utl::PrintUtlMatrix(*this->m_SPFToFeatureTransform, "this->m_SPFToFeatureTransform");
}

template <class TInputImage, class TOutputImage>
void
ProfileFromSPFImageFilter<TInputImage,TOutputImage>
::BeforeThreadedGenerateData()
{
  utlShowPosition(this->GetDebug());
  // this->Print(std::cout<<"this=");
  this->SetSPFIEstimator();

  if (this->m_Orientations->Rows()>0)
    {
    this->m_SPFIEstimator->GetSamplingSchemeQSpace()->SetOrientationsSpherical(this->m_Orientations);
    if (m_RadiusVector->size()==0)
      {
      utlException(this->m_Radius<0, "need to set m_Radius or m_RadiusVector");
      utl::MatchBVectorAndGradientMatrix(this->m_Radius, *m_RadiusVector, *this->m_Orientations);
      }
    this->m_SPFIEstimator->GetSamplingSchemeQSpace()->SetTau(this->m_Tau);
    this->m_SPFIEstimator->SetSHRank(this->m_SHRank);
    this->m_SPFIEstimator->SetRadialRank(this->m_RadialRank);
    this->m_SPFIEstimator->SetBasisScale(this->m_BasisScale);
    this->m_SPFIEstimator->SetScaleImage(this->m_ScaleImage);
    this->m_SPFIEstimator->SetDebug(this->GetDebug());

    bool needToConvertToQValue = (this->m_IsInQSpace && !this->m_IsFourier) || (!this->m_IsInQSpace && this->m_IsFourier);
    if (needToConvertToQValue)
      this->m_SPFIEstimator->GetSamplingSchemeQSpace()->SetBVector(this->m_RadiusVector);
    else
      this->m_SPFIEstimator->GetSamplingSchemeQSpace()->SetRadiusVector(this->m_RadiusVector);
    }

    
  // this->m_SPFIEstimator->Print(std::cout<<"m_SPFIEstimator: ");
  
  Superclass::BeforeThreadedGenerateData();
}

template <class TInputImage, class TOutputImage>
void
ProfileFromSPFImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData( const typename TOutputImage::RegionType &outputRegionForThread, ThreadIdType threadId)
{
  utlShowPosition(this->GetDebug());
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
  
  unsigned int outputDim = outputPtr->GetNumberOfComponentsPerPixel();
  outputPixel.SetSize(outputDim);
  unsigned int inputDim = inputPtr->GetNumberOfComponentsPerPixel();
  inputPixel.SetSize(inputDim);

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  VectorType spfVec, temp;

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
          std::cout << "inputIndex = " << inputIndex << std::endl << std::flush;
        if (!IsImageEmpty(this->m_ScaleImage))
          selfClone->SetBasisScale(this->m_ScaleImage->GetPixel(inputIndex));
        if (selfClone->m_SPFToFeatureTransform->Size()==0)
          selfClone->ComputeSPFToFeatureTransform();
        spfVec = utl::VariableLengthVectorToUtlVector(inputPixel);
        // utl::PrintUtlVector(spfVec, "spfVec");
        // utl::PrintUtlMatrix(*selfClone->m_SPFToFeatureTransform, "*selfClone->m_SPFToFeatureTransform");

        utl::ProductUtlMv(*selfClone->m_SPFToFeatureTransform, spfVec, temp);
        outputPixel = utl::UtlVectorToVariableLengthVector(temp);
        if (this->GetDebug())
          {
          if (!IsImageEmpty(this->m_ScaleImage))
            std::cout << "scale = " << this->m_ScaleImage->GetPixel(inputIndex) << std::endl << std::flush;
          itk::PrintVariableLengthVector(inputPixel, "inputPixel");
          itk::PrintVariableLengthVector(outputPixel, "outputPixel");
          }
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
ProfileFromSPFImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent
    << ", m_Radius = " << this->m_Radius
    << std::endl;
  utl::PrintVector(*m_RadiusVector, "m_RadiusVector", " ", os<<indent);
}

}


#endif 


