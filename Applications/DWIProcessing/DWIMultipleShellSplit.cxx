/**
 *       @file  DWIMultipleShellSplit.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-31-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#include "utl.h"
#include "DWIMultipleShellSplitCLP.h"

/**
 * \brief  Split mutiple shell DWI data into several single shell data.
 */
int 
main (int argc, char const* argv[])
{
  
  PARSE_ARGS;

  typedef float PixelType;
  typedef itk::Image<PixelType, 4> InputImageType;
  typedef InputImageType OutputImageType;
  typedef itk::Image<PixelType, 3> B0ImageType;
  typedef itk::Image<PixelType, 3> MaskImageType;
  
  InputImageType::Pointer inputImage;
  B0ImageType::Pointer b0Image;
  MaskImageType::Pointer mask;
  OutputImageType::Pointer outputImage = OutputImageType::New();

  if (_MaskArg.isSet())
    {
    itk::ReadImage<MaskImageType>(_Mask, mask);
    }

  std::vector<PixelType> bVector, bVectorSeparated;
  utl::ReadVector(_BValueFile, bVector);

  std::vector<std::vector<int> > separateIndex;
  separateIndex = utl::SeparateVector(bVector, bVectorSeparated, (PixelType)_Threshold);
  int b0Index=-1;
  for ( int i = 0; i < bVectorSeparated.size(); i += 1 ) 
    {
    if (bVectorSeparated[i]<_Threshold)
      b0Index = i;
    }
  
  std::string fileNoExt, ext, file;
      
  if (_OutputBValueFileArg.isSet())
    {
    utl::GetFileExtension(_OutputBValueFile, ext, fileNoExt);
    for ( int m = 0; m < separateIndex.size(); m += 1 ) 
      {
      std::vector<PixelType> bVectorSingleShell = utl::SelectVector(bVector, separateIndex[m]);
      file = fileNoExt + "_b" + utl::ConvertNumberToString(bVectorSeparated[m]) + "." + ext;
      utl::SaveVector(bVectorSingleShell, file);
      }
    }

  utlException(!_GradientFileArg.isSet() && _OutputGradientFileArg.isSet(), "Need to set the gradients file");
  utlException(_GradientFileArg.isSet() && !_OutputGradientFileArg.isSet(), "Need to set the gradients file");
  if (_GradientFileArg.isSet())
    {
    std::vector<std::vector<std::string> > gradVector, gradVectorSingleShell;
    utl::ReadLinesFirstlineCheck(_GradientFile, gradVector);
    utlException(gradVector.size()!=bVector.size(), "wrong size for the gradient file");
    utl::GetFileExtension(_OutputGradientFile, ext, fileNoExt);
    for ( int m = 0; m < separateIndex.size(); m += 1 ) 
      {
      gradVectorSingleShell = utl::SelectVector(gradVector, separateIndex[m]);
      file = fileNoExt + "_b" + utl::ConvertNumberToString(bVectorSeparated[m]) + "." + ext;
      utl::Save2DVector(gradVectorSingleShell, file);
      }
    }

 
  itk::ReadImage<InputImageType>(_InputFile, inputImage);
  
  InputImageType::PixelType    inputImagePixelValue;
  InputImageType::IndexType    inputImagePixelIndex;
  InputImageType::RegionType inRegion = inputImage->GetLargestPossibleRegion();
  InputImageType::SizeType inputImageSize = inRegion.GetSize();
  unsigned int numberOfVolumes = inputImageSize[3];
  utlException(b0Index<0 && numberOfVolumes!=bVector.size()+1, "the first volume should be b0 image if there is no 0 in b values");
  utlException(b0Index>=0 && numberOfVolumes!=bVector.size(), "wrong size!, bVector.size()="<<bVector.size()<<", numberOfVolumes="<<numberOfVolumes);
  // utlException(numberOfVolumes!=bVector.size() && numberOfVolumes!=bVector.size()+1, "wrong size!, bVector.size()="<<bVector.size()<<", numberOfVolumes="<<numberOfVolumes);
  
  B0ImageType::IndexType b0ImagePixelIndex;
  MaskImageType::IndexType maskPixelIndex;
  utlException(_Setb0InEachShellArg.isSet() && _B0NormalizationArg.isSet(), "--b0InEachShell and --b0normalization can not be set simultaneously");
  if (_Setb0InEachShellArg.isSet() || _B0NormalizationArg.isSet() || _B0OutputFileArg.isSet())
    {
    b0Image = B0ImageType::New();

    B0ImageType::PixelType b0ImagePixelValue;
    itk::CopyImageInformation<InputImageType, B0ImageType>(inputImage, b0Image);
    b0Image->Allocate();
    b0Image->FillBuffer(0);

    for (unsigned int j = 0; j < inputImageSize[1]; j++)
      for (unsigned int i = 0; i < inputImageSize[0]; i++)
        for (unsigned int k = 0; k < inputImageSize[2]; k++)
          {
          inputImagePixelIndex[0] = i;
          inputImagePixelIndex[1] = j;
          inputImagePixelIndex[2] = k;

          b0ImagePixelIndex[0] = i;
          b0ImagePixelIndex[1] = j;
          b0ImagePixelIndex[2] = k;

          if (b0Index>=0)
            {
            b0ImagePixelValue = 0.0;
            for ( int t = 0; t < separateIndex[b0Index].size(); t += 1 ) 
              {
              inputImagePixelIndex[3] = separateIndex[b0Index][t];
              b0ImagePixelValue += inputImage->GetPixel(inputImagePixelIndex);
              }
            b0Image->SetPixel(b0ImagePixelIndex, b0ImagePixelValue/(PixelType)separateIndex[b0Index].size());
            }
          else
            {
            inputImagePixelIndex[3] = 0;
            b0Image->SetPixel(b0ImagePixelIndex, inputImage->GetPixel(inputImagePixelIndex));
            }
          }
    }

  if (_B0OutputFileArg.isSet())
    itk::SaveImage<B0ImageType>(b0Image, _B0OutputFile, "Write b0 Image to");

  OutputImageType::RegionType  outputImageRegion;
  OutputImageType::SizeType    outputImageSize;
  OutputImageType::PixelType   outputImagePixelValue;
//  OutputImageType::PixelType   OutputImageZeroPixelValue;
  OutputImageType::IndexType   outputImagePixelIndex;
//  OutputImageType::SpacingType OutputImageSpacing;
  
  outputImage->CopyInformation(inputImage);
  outputImageSize[0] = inputImageSize[0];
  outputImageSize[1] = inputImageSize[1];
  outputImageSize[2] = inputImageSize[2];
  // outputImageSize[3] = NumberOfDWIVolumes;
  
  std::cout << "Number of Volumes: " << numberOfVolumes << std::endl;

  int offsetInput = numberOfVolumes==bVector.size() ? 0 : 1;
  int offsetOutput = _Setb0InEachShell? 1 : 0;


  for ( int m = 0; m < separateIndex.size(); m += 1 ) 
    {
    outputImageSize[3] = separateIndex[m].size() + offsetOutput;
    outputImageRegion.SetSize(outputImageSize);
    outputImage->SetRegions(outputImageRegion);
    outputImage->Allocate();

    for (unsigned int i = 0; i < inputImageSize[0]; i++)
      for (unsigned int j = 0; j < inputImageSize[1]; j++)
        for (unsigned int k = 0; k < inputImageSize[2]; k++)
          {
          inputImagePixelIndex[0] = i;
          inputImagePixelIndex[1] = j;
          inputImagePixelIndex[2] = k;

          outputImagePixelIndex[0] = i;
          outputImagePixelIndex[1] = j;
          outputImagePixelIndex[2] = k;

          PixelType mm=1;
          if (_MaskArg.isSet())
            {
            maskPixelIndex[0] = i;
            maskPixelIndex[1] = j;
            maskPixelIndex[2] = k;
            mm = mask->GetPixel(maskPixelIndex);
            }

          PixelType bb=0.0;
          if ( mm>0 && (_B0NormalizationArg.isSet() || _Setb0InEachShellArg.isSet() ) )
            {
            b0ImagePixelIndex[0] = i;
            b0ImagePixelIndex[1] = j;
            b0ImagePixelIndex[2] = k;
            bb = b0Image->GetPixel(b0ImagePixelIndex);
            }

          for ( int n = 0; n < separateIndex[m].size(); n += 1 ) 
            {
            inputImagePixelIndex[3] = separateIndex[m][n] + offsetInput;
            outputImagePixelIndex[3] = n + offsetOutput;

            if (_B0NormalizationArg.isSet())
              {
              if (bb>0 && mm>0)
                outputImage->SetPixel(outputImagePixelIndex, inputImage->GetPixel(inputImagePixelIndex)/bb);
              else
                outputImage->SetPixel(outputImagePixelIndex, 0.0);
              }
            else if (mm>0)
              outputImage->SetPixel(outputImagePixelIndex, inputImage->GetPixel(inputImagePixelIndex));
            else
              outputImage->SetPixel(outputImagePixelIndex, 0.0);
            }

          if (_Setb0InEachShellArg.isSet())
            {
            outputImagePixelIndex[3] = 0;
            outputImage->SetPixel(outputImagePixelIndex, bb);
            }
          }
    utl::GetFileExtension(_OutputFile, ext, fileNoExt);
    file = fileNoExt + "_b" + utl::ConvertNumberToString(bVectorSeparated[m]) + "." + ext;
    itk::SaveImage<OutputImageType>(outputImage, file);
    }

  return 0;
}
