/**
 *       @file  itkFunctorFromStringImageFilter.hxx
 *      @brief  
 *     Created  "04-06-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef itkFunctorFromStringImageFilter_hxx
#define itkFunctorFromStringImageFilter_hxx

#include "itkFunctorFromStringImageFilter.h"
#include "utlCoreMacro.h"
#include "itkVectorImageRegionIteratorWithIndex.h"
#include "utlITK.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage, class TMaskImage  >
FunctorFromStringImageFilter< TInputImage, TOutputImage, TMaskImage >
::FunctorFromStringImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);
}

template< typename TInputImage, typename TOutputImage, class TMaskImage >
void
FunctorFromStringImageFilter< TInputImage, TOutputImage, TMaskImage >
::GenerateOutputInformation()
{
  itkShowPositionThreadedLogger(utl::IsLogDebug(this->m_LogLevel));
  typename TInputImage::ConstPointer inputImage = this->GetInput();
  typename TInputImage::Pointer inputTmp = TInputImage::New();
  OutputImagePointer outImage = this->GetOutput();

  // copy information
  itk::CopyImageInformation(inputImage, outImage);

  int numberOfInputs = this->GetNumberOfInputs();

  std::vector<std::vector<int> > size4d;
  for ( int i = 0; i < numberOfInputs; ++i ) 
    {
    inputTmp->Graft(this->GetInput(i));
    std::vector<int> tmpVec = itk::GetVectorImageFullSize(inputTmp);
    size4d.push_back(tmpVec);

    if (i!=0)
      {
      for ( int kk = 0; kk < 4; ++kk ) 
        {
        utlSAGlobalException(size4d[0][kk]!=size4d[i][kk])
          (this->m_VectorAxis)(i)(kk).msg("wrong 4d size (x,y,z,t)");
        }
      }
    }
}

template< typename TInputImage, typename TOutputImage, class TMaskImage >
void
FunctorFromStringImageFilter< TInputImage, TOutputImage, TMaskImage >
::BeforeThreadedGenerateData()
{
  itkShowPositionThreadedLogger(utl::IsLogDebug(this->m_LogLevel));

  this->VerifyInputParameters();

  if (utl::IsLogDebug(this->m_LogLevel) && this->GetNumberOfThreads()>1)
    this->CreateLoggerVector();
  
  utlException(this->m_VectorAxis!=3, "m_VectorAxis should be 3, for element-wise function");

  // typename TInputImage::ConstPointer inputImage = this->GetInput();
  // std::vector<int> size = itk::GetVectorImageFullSize(inputImage);
  // this->m_Functor.VerifyInputParameters(size[this->m_VectorAxis]);  
  // this->m_Functor.Initialize();  
}

/**
 * ThreadedGenerateData Performs the pixel-wise addition
 */
template< typename TInputImage, typename TOutputImage, class TMaskImage  >
void
FunctorFromStringImageFilter< TInputImage, TOutputImage, TMaskImage >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  itkShowPositionThreadedLogger(utl::IsLogDebug(this->m_LogLevel));

  InputImagePointer inputImage = const_cast<InputImageType *>(this->GetInput());
  OutputImagePointer outImage = this->GetOutput();
  int numberOfInputs = this->GetNumberOfInputs();
  
  // outImage->Print(std::cout << "outImage =\n");

  Pointer selfClone = this->Clone();
  selfClone->m_ThreadID = threadId;
  std::string threadIDStr =  selfClone->ThreadIDToString();

  std::vector<double> x(numberOfInputs);
  exprtk::symbol_table<double> symbol_table;
  symbol_table.add_vector("x",x);
  symbol_table.add_constants();

  exprtk::expression<double> expression;
  expression.register_symbol_table(symbol_table);

  exprtk::parser<double> parser;
  // parser.compile(m_Expression, expression);

  InputImageRegionType regionInput;
  itk::CopyImageRegion(outputRegionForThread, regionInput);

  std::vector<VectorImageRegionIteratorWithIndex<InputImageType> > inputItVec;
  for ( int i = 0; i < numberOfInputs; ++i ) 
    {
    InputImagePointer inputTmp = const_cast<InputImageType *>(this->GetInput(i));

    VectorImageRegionIteratorWithIndex<InputImageType> it(inputTmp, regionInput, this->m_VectorAxis);
    inputItVec.push_back(it);
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

  if (utl::IsLogDebug(this->m_LogLevel))
    {
    std::ostringstream msg;
    msg << threadIDStr << "outputRegionForThread = " << outputRegionForThread << std::endl << std::flush;
    msg << threadIDStr << "regionInput = " << regionInput << std::endl << std::flush;
    this->WriteLogger(msg.str());
    }

  int vecInputSize = itk::GetVectorImageVectorSize(inputImage);
  int outDim = itk::GetVectorImageVectorSize(outImage, this->m_VectorAxis);
  utlSAGlobalException(maskDim>1 && maskDim!=vecInputSize)(maskDim)(vecInputSize).msg("4D mask have a wrong 4-th dimension.");


  itk::VectorImageRegionIteratorWithIndex<TOutputImage> outIt(outImage, outputRegionForThread, this->m_VectorAxis);   

  OutputImageIndexType index;
  VariableLengthVector<double> inPixel, outPixel, maskPixel;
  std::vector<utl::Vector<double> > pixelVec(numberOfInputs);
  std::vector<double> inStdVec(numberOfInputs);
  utl::Vector<double> outVec(outDim);
  outPixel.SetSize(outDim);
  outPixel.Fill(0.0);
    
  for ( int n = 0; n < numberOfInputs; n += 1 ) 
    inputItVec[n].GoToBegin();
  outIt.GoToBegin();
  maskIt.GoToBegin();

  while ( !outIt.IsAtEnd() ) 
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
        for ( int n = 0; n < numberOfInputs; n += 1 ) 
          ++inputItVec[n];
        ++outIt;
        ++maskIt;
        continue; 
        }
      }
    
    index = outIt.GetIndex();
    
    for ( int i = 0; i < (this->m_VectorAxis==3?1:vecInputSize); ++i ) 
      {

      if (this->IsMaskUsed() && maskDim>1)
        {
        maskIt.GetVector(maskPixel, i);
        if (maskPixel.GetSquaredNorm()==0)
          {
          outPixel.Fill(0.0); 
          outIt.SetVector(outPixel, i); 
          for ( int n = 0; n < numberOfInputs; n += 1 ) 
            ++inputItVec[n];
          ++outIt;
          ++maskIt;
          continue; 
          }
        }

      for ( int n = 0; n < numberOfInputs; ++n ) 
        {
        inputItVec[n].GetVector(inPixel, i);
        pixelVec[n] = utl::VariableLengthVectorToUtlVector(inPixel);
        }

      if (utl::IsLogDebug(this->m_LogLevel))
        {
        std::ostringstream msg;
        msg << "\n" << threadIDStr << "index = " << index << ", i=" << i << std::endl << std::flush;
        for ( int n = 0; n < numberOfInputs; ++n ) 
          utl::PrintUtlVector(pixelVec[n], "pixelVec[n]", " ", msg << threadIDStr);
        this->WriteLogger(msg.str());
        }

      for ( int v = 0; v < inPixel.Size(); ++v ) 
        {
        for ( int n = 0; n < numberOfInputs; n += 1 ) 
          {
          inStdVec[n] = pixelVec[n][v];
          }
        symbol_table.get_vector("x")->operator=(inStdVec);
        parser.compile(m_Expression,expression);
        outVec[v] = expression.value();
        }

      // outVec = selfClone->m_Functor(pixelVec);
      outPixel = utl::UtlVectorToVariableLengthVector(outVec);

      if (utl::IsLogDebug(this->m_LogLevel))
        {
        std::ostringstream msg;
        itk::PrintVariableLengthVector(outPixel, "outPixel", " ", msg << threadIDStr);
        this->WriteLogger(msg.str());
        }
      outIt.SetVector(outPixel,i);
      }
    
    for ( int n = 0; n < numberOfInputs; n += 1 ) 
      ++inputItVec[n];
    ++outIt;
    ++maskIt;
    }

}

} // end namespace itk

#endif




