/**
 *       @file  itkUnaryFunctorVectorImageFilter.hxx
 *      @brief  
 *     Created  "09-11-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef itkUnaryFunctorVectorImageFilter_hxx
#define itkUnaryFunctorVectorImageFilter_hxx

#include "itkUnaryFunctorVectorImageFilter.h"
#include "utlCoreMacro.h"
#include "itkVectorImageRegionIteratorWithIndex.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage, typename TFunction, class TMaskImage  >
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage >
::UnaryFunctorVectorImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);
}

template< typename TInputImage, typename TOutputImage, typename TFunction, class TMaskImage >
void
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage >
::GenerateOutputInformation()
{
  utlVLogPosition(LOG_DEBUG);
  typename TInputImage::ConstPointer inputImage = this->GetInput();
  OutputImagePointer outImage = this->GetOutput();

  // copy information
  itk::CopyImageInformation(inputImage, outImage);
  utlException(this->m_VectorAxis<0, "need to set non-negative axis");

  int vecInputSize = itk::GetVectorImageVectorSize(inputImage);
  int outVectorSize = this->m_Functor.GetOutputDimension(vecInputSize);
  int outDim= (this->m_VectorAxis==3) ? outVectorSize : vecInputSize;
  if (this->m_VectorAxis!=3)
    {
    OutputImageRegionType regionTmp = outImage->GetLargestPossibleRegion();
    OutputImageSizeType sizeTmp = regionTmp.GetSize();
    sizeTmp[this->m_VectorAxis] = this->m_Functor.GetOutputDimension(sizeTmp[this->m_VectorAxis]);
    regionTmp.SetSize(sizeTmp);
    outImage->SetRegions(regionTmp);
    }
  utlSAException(outDim<=0)(outDim).msg("wrong outDim");
  itk::SetVectorImageVectorSize(outImage, outDim);
}

template< typename TInputImage, typename TOutputImage, typename TFunction, class TMaskImage >
void
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage >
::BeforeThreadedGenerateData()
{
  itkShowPositionThreadedLogger(utl::IsLogDebug(this->m_LogLevel));

  this->VerifyInputParameters();

  if (utl::IsLogDebug(this->m_LogLevel) && this->GetNumberOfThreads()>1)
    this->CreateLoggerVector();

  typename TInputImage::ConstPointer inputImage = this->GetInput();
  std::vector<int> size = itk::GetVectorImageFullSize(inputImage);
  this->m_Functor.VerifyInputParameters(size[this->m_VectorAxis]);  
  this->m_Functor.Initialize();  
}

/**
 * ThreadedGenerateData Performs the pixel-wise addition
 */
template< typename TInputImage, typename TOutputImage, typename TFunction, class TMaskImage  >
void
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction, TMaskImage >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  itkShowPositionThreadedLogger(utl::IsLogDebug(this->m_LogLevel));
  InputImagePointer inputImage = const_cast<InputImageType *>(this->GetInput());
  OutputImagePointer outImage = this->GetOutput();

  Pointer selfClone = this->Clone();
  selfClone->m_ThreadID = threadId;
  std::string threadIDStr =  selfClone->ThreadIDToString();

  InputImageRegionType regionInput;
  itk::CopyImageRegion(outputRegionForThread, regionInput);
  itk::VectorImageRegionIteratorWithIndex<InputImageType> it(inputImage, regionInput, this->m_VectorAxis);
  // itk::VectorImageRegionIteratorWithIndex<InputImageType> it(inputImage, outputRegionForThread, m_VectorAxis);
        
  if (utl::IsLogDebug(this->m_LogLevel))
    {
    std::ostringstream msg;
    selfClone->Print(msg<< threadIDStr << "selfClone = \n");
    this->WriteLogger(msg.str());
    }

  VectorImageRegionIteratorWithIndex<MaskImageType> maskIt;
  int maskDim = 1;
  if (this->IsMaskUsed())
    {
    typename MaskImageType::RegionType regionTmp;
    itk::CopyImageRegion(outputRegionForThread, regionTmp, 1);
    maskIt = VectorImageRegionIteratorWithIndex<MaskImageType>(this->m_MaskImage, regionTmp, this->m_VectorAxis);
    maskDim = itk::GetVectorImageVectorSize(this->m_MaskImage);
    }

  int vecInputSize = itk::GetVectorImageVectorSize(inputImage);
  int outDim = itk::GetVectorImageVectorSize(outImage, this->m_VectorAxis);
  utlSAGlobalException(maskDim>1 && maskDim!=vecInputSize)(maskDim)(vecInputSize).msg("4D mask have a wrong 4-th dimension.");

  itk::VectorImageRegionIteratorWithIndex<TOutputImage> outIt(outImage, outputRegionForThread, this->m_VectorAxis);   

  InputImageIndexType index;
  VariableLengthVector<double> inPixel, outPixel, maskPixel;
  utl::Vector<double> inVec, outVec(outDim);
  outPixel.SetSize(outDim);
  outPixel.Fill(0.0);
  for (it.GoToBegin(), maskIt.GoToBegin(), outIt.GoToBegin(); 
    !it.IsAtEnd(); 
    ++it, ++maskIt, ++outIt)
    {

    if (this->IsMaskUsed() && maskDim==1)
      {
      maskIt.GetVector(maskPixel,0);
      if (maskPixel.GetSquaredNorm()==0)
        {
        outPixel.Fill(0.0); 
        for ( int i = 0; i < (this->m_VectorAxis==3?1:vecInputSize); ++i ) 
          {
          outIt.SetVector(outPixel, i); 
          }
        continue; 
        }
      }

    index = it.GetIndex();
    // if (IsOutsideBox(index, realBox))
    //   { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }


    for ( int i = 0; i < (this->m_VectorAxis==3?1:vecInputSize); ++i ) 
      {

      if (this->IsMaskUsed() && maskDim>1)
        {
        maskIt.GetVector(maskPixel, i);
        if (maskPixel.GetSquaredNorm()==0)
          {
          outPixel.Fill(0.0); 
          outIt.SetVector(outPixel, i); 
          continue; 
          }
        }

      it.GetVector(inPixel, i);

      inVec = utl::VariableLengthVectorToUtlVector(inPixel);
      if (utl::IsLogDebug(this->m_LogLevel))
        {
        std::ostringstream msg;
        msg << "\n" << threadIDStr << "index = " << index << std::endl << std::flush;
        itk::PrintVariableLengthVector(inPixel, "inPixel", " ", msg << threadIDStr);
        this->WriteLogger(msg.str());
        }

      outVec = selfClone->m_Functor(inVec);
      outPixel = utl::UtlVectorToVariableLengthVector(outVec);

      if (utl::IsLogDebug(this->m_LogLevel))
        {
        std::ostringstream msg;
        itk::PrintVariableLengthVector(outPixel, "outPixel", " ", msg << threadIDStr);
        this->WriteLogger(msg.str());
        }
      outIt.SetVector(outPixel,i);
      }
    }
}

} // end namespace itk

#endif


