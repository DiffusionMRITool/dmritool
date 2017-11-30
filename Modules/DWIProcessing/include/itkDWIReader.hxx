/**
 *       @file  itkDWIReader.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-11-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkDWIReader_hxx
#define __itkDWIReader_hxx

#include "itkDWIReader.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "utl.h"


namespace itk
{

template <class TPixelType, unsigned int VImageDimension>
DWIReader<TPixelType, VImageDimension>
::DWIReader() 
{
  m_NormalizeDWI = true;
  // m_DWIWithB0 = false;
  m_IsInput4DImage = true;
  m_CorrectDWIValues = true;
  m_ShowWarnings = true;

  m_MaskImage = MaskImageType::New();
  m_B0Image = B0ImageType::New();

  m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
}

template <class TPixelType, unsigned int VImageDimension>
bool
DWIReader<TPixelType, VImageDimension>
::DetermineIsInput4DImage(const std::string& dataStr)
{
  typedef Image<TPixelType,4> DWI4DImageType;
  typename DWI4DImageType::Pointer dwi4DTempImage;
  itk::ReadImageInformation<DWI4DImageType>(dataStr, dwi4DTempImage);
  typename DWI4DImageType::SizeType size = dwi4DTempImage->GetLargestPossibleRegion().GetSize();
  int numberOfDWI = size[3];
  if (numberOfDWI>1)
    return true;  // 4D Image
  else
    {
    typename DWIImageType::Pointer dwiTempImage;
    itk::ReadImageInformation<DWIImageType>(dataStr, dwiTempImage);
    numberOfDWI = dwiTempImage->GetNumberOfComponentsPerPixel();
    if (numberOfDWI>1)
      return false; // VectorImage
    else
      return true;  // 4D Image, one component
    }
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIReader<TPixelType, VImageDimension>
::ReadFromConfigurationFile(const std::string& file)
{
  utlShowPosition(this->GetDebug());

  // read file
  utlGlobalException(m_ConfigurationFile=="", "need to set the configuration file!");
  std::vector<std::vector<std::string> >  stringMatrix;
  utl::ReadLinesFirstlineCheck(m_ConfigurationFile, stringMatrix);

  // check configuration file and dwi size
  std::string dwiPath, dwiFile;
  utl::GetPath(file, dwiPath, dwiFile);
  double bvalue=-1;
  std::string dataStr, gradStr, indexStr, bStr;
  std::vector<std::vector<std::string> > gradStrTempMatrix;
  std::vector<int> numDWINob0, additionalB0;
  std::vector<std::vector<double> > bVec; 
  std::vector<std::vector<std::vector<std::string> > > gradStrMatrix;
  std::vector<std::vector<int> > b0IndexVec;
  std::vector<std::vector<int> > selectIndexVec;

  typedef ImageFileReader<DWIImageType> DWIReaderType;
  typename DWIImageType::Pointer dwiTempImage;
  
  typedef Image<TPixelType, 4>  DWI4DImageType;
  typedef ImageFileReader<DWI4DImageType> DWI4DReaderType;
  typename DWI4DImageType::Pointer dwi4DTempImage;

  bool isB0Set = !IsImageEmpty(m_B0Image);
  double bThresholdSingleShell = m_SamplingSchemeQSpace->GetBThresholdSingleShell();

  // Read Input File Info
  int numberOfDWIsWithoutB0=0;
  int numberOfB0=0;
  std::vector<std::string> ossStr;
  std::string dataStr0;
  for ( int i = 0; i < stringMatrix.size(); i += 1 ) 
    {
    // utlPrintVar1(true, i);
    utlGlobalException(stringMatrix[i].size()!=3 && stringMatrix[i].size()!=4, "wrong size of configuration file");
      
    dataStr = utl::IsFileExist(stringMatrix[i][2]) ? stringMatrix[i][2] : dwiPath+stringMatrix[i][2] ;
    gradStr = utl::IsFileExist(stringMatrix[i][1]) ? stringMatrix[i][1] : dwiPath+stringMatrix[i][1] ;

    if (i==0)
      {
      dataStr0 = utl::IsFileExist(stringMatrix[i][2]) ? stringMatrix[i][2] : dwiPath+stringMatrix[i][2];
      utlGlobalException( !IsImageEmpty(m_MaskImage) && !itk::VerifyImageSize<MaskImageType>(m_MaskImage, dataStr0, true), "wrong size of m_MaskImage and DWI file " << stringMatrix[i][2] << ". Inconsistent information!");
      utlGlobalException( !IsImageEmpty(m_B0Image) && !itk::VerifyImageSize<B0ImageType>(m_B0Image, dataStr0, true), "wrong size of b0Image and DWI file " << stringMatrix[i][2] << ". Inconsistent information! m_B0Image="<<m_B0Image);
      m_IsInput4DImage = DetermineIsInput4DImage(dataStr0);
      }
    else
      {
      utlGlobalException( !itk::VerifyImageInformation(dataStr0, dataStr, true), "DWI files "<< dataStr << " and " << dataStr0 << " have inconsistent information!"); 
      }

    if (stringMatrix[i].size()==4)
      {
      indexStr = utl::IsFileExist(stringMatrix[i][3]) ? stringMatrix[i][3] : dwiPath+stringMatrix[i][3];
      std::vector<int> tmpVec; 
      utl::ReadVector(indexStr,tmpVec);
      selectIndexVec.push_back(tmpVec);
      }

    utl::ReadLinesFirstlineCheck(gradStr, gradStrTempMatrix);
    std::vector<double> bVecTemp; 
    if (utl::IsNumber(stringMatrix[i][0]))
      {
      std::istringstream ( stringMatrix[i][0] ) >> bvalue;
      bVecTemp.clear();
      bVecTemp.assign(gradStrTempMatrix.size(), bvalue);
      }
    else 
      {
      bStr = utl::IsFileExist(stringMatrix[i][0]) ? stringMatrix[i][0] : dwiPath+stringMatrix[i][0];
      utl::ReadVector(bStr,bVecTemp);
      utlGlobalException(gradStrTempMatrix.size()!=gradStrTempMatrix.size(), "wrong size of b values in " << bStr << " and gradients in " << gradStr);
      }
    bVec.push_back(bVecTemp);

    std::vector<int> b0IndexTemp; 
    for ( int j = 0; j < gradStrTempMatrix.size(); j += 1 ) 
      {
      utlGlobalException(gradStrTempMatrix[j].size()!=3, "wrong size of gradients " << gradStr);
      if (stringMatrix[i].size()==4 && !utl::IsInVector(selectIndexVec[i],j))
        continue;
      double gradx = utl::ConvertStringToNumber<double>(gradStrTempMatrix[j][0]);
      double grady = utl::ConvertStringToNumber<double>(gradStrTempMatrix[j][1]);
      double gradz = utl::ConvertStringToNumber<double>(gradStrTempMatrix[j][2]);
      bool gradIsZero = std::abs(gradx)<1e-4 && std::abs(grady)<1e-4 && std::abs(gradz)<1e-4;
      bool bIsZero = std::abs(bVecTemp[j])<(bThresholdSingleShell>0?bThresholdSingleShell:1e-4);
      // std::cout << "bThresholdSingleShell=" << bThresholdSingleShell << ", bIsZero="<<bIsZero << std::endl << std::flush;
      // utlGlobalException(gradIsZero && !bIsZero || !gradIsZero && bIsZero, "Wrong logic in gradients and b values! bValue = "<< bVecTemp[j] << ", grad=" <<"("<<gradx<<","<<grady<<","<<gradz<<").");
      if (gradIsZero || bIsZero)
        b0IndexTemp.push_back(j);
      }
    utlGlobalException(isB0Set && b0IndexTemp.size()>0, "Since b0 is mannually set already, the configuration file should not have b0 image");
    gradStrMatrix.push_back(gradStrTempMatrix);
    b0IndexVec.push_back(b0IndexTemp);

    std::ostringstream oss;
    if (utl::IsNumber(stringMatrix[i][0]))
      oss << "> b="<< bvalue; 
    else 
      oss << "> bStr="<< bStr; 
    oss << ", gradStrTempMatrix.size()=" << gradStrTempMatrix.size() << ", 4D image = " << dataStr;
    if (stringMatrix[i].size()==4)
      oss << ", index=" << indexStr;
    // oss << std::endl << std::flush;
    ossStr.push_back(oss.str());


    int numberOfDWI;
    if (m_IsInput4DImage)
      {
      itk::ReadImageInformation<DWI4DImageType>(dataStr, dwi4DTempImage);
      typename DWI4DImageType::SizeType size = dwi4DTempImage->GetLargestPossibleRegion().GetSize();
      numberOfDWI = size[3];
      }
    else
      {
      itk::ReadImageInformation<DWIImageType>(dataStr, dwiTempImage);
      numberOfDWI = dwiTempImage->GetNumberOfComponentsPerPixel();
      }
    
    utlGlobalException(numberOfDWI<gradStrTempMatrix.size(), "The number of gradients in " << gradStr << " is "<< gradStrTempMatrix.size() << " more than number of DWIs in " << dataStr << "="<<numberOfDWI);
    additionalB0.push_back(numberOfDWI - gradStrTempMatrix.size());
    if (stringMatrix[i].size()==3)
      {
      numDWINob0.push_back(gradStrTempMatrix.size() - b0IndexTemp.size());
      numberOfDWIsWithoutB0 += gradStrTempMatrix.size() - b0IndexTemp.size();
      numberOfB0 += b0IndexTemp.size() + numberOfDWI - gradStrTempMatrix.size();
      }
    else
      {
      numDWINob0.push_back(selectIndexVec[i].size() - b0IndexTemp.size());
      numberOfDWIsWithoutB0 += selectIndexVec[i].size() - b0IndexTemp.size();
      numberOfB0 += b0IndexTemp.size() + numberOfDWI - gradStrTempMatrix.size();
      }

    // numDWInob0 records the start index for i-th dwi file
    if (i>0)
      numDWINob0[i] += numDWINob0[i-1];
    }

  // utlGlobalException( !this->GetB0Image() && numberOfB0==0, "no b0 images in " << file);
  std::cout << numberOfB0 << " b0 images, " << numberOfDWIsWithoutB0 << " DWIs without b0. "<< std::endl << std::flush;
  if (this->GetDebug())
    {
    utl::PrintVector(numDWINob0, "numDWINob0");
    utl::PrintVector(additionalB0, "additionalB0");
    for ( int jj = 0; jj < b0IndexVec.size(); jj += 1 ) 
      utl::PrintVector(b0IndexVec[jj], "b0IndexVec[jj]");
    }

  // Read data in input file 
  typename DWIImageType::Pointer dwiImage = this->GetOutput();
  STDVectorPointer bVector = m_SamplingSchemeQSpace->GetBVector();
  MatrixPointer orientationsCartesian = m_SamplingSchemeQSpace->GetOrientationsCartesian();
  bVector->clear();
  orientationsCartesian->ReSize(numberOfDWIsWithoutB0, 3);
  typename B0ImageType::IndexType  b0Index;
  typename B0ImageType::PixelType  b0Pixel;
  PixelType dwiPixel, dwiZeroPixel;
  dwiPixel.SetSize(numberOfDWIsWithoutB0);
  dwiZeroPixel.SetSize(numberOfDWIsWithoutB0);
  dwiZeroPixel.Fill(0.0);

  for ( int i = 0; i < stringMatrix.size(); i += 1 ) 
    {
    if (stringMatrix.size()==1)
      std::cout << "loading the " << i+1 << "th DWI data from total " << stringMatrix.size() << " DWI data" << std::endl;
    std::cout << ossStr[i] << std::endl << std::flush;
    dataStr = utl::IsFileExist(stringMatrix[i][2]) ? stringMatrix[i][2] : dwiPath+stringMatrix[i][2];

    // b vector, gradients
    // utl::PrintVector(bVec[i], "bVec[i]");
    for ( int j = 0; j < bVec[i].size(); j += 1 ) 
      {
      if (!utl::IsInVector(b0IndexVec[i], j) && ( stringMatrix[i].size()==3 || (stringMatrix[i].size()==4 && utl::IsInVector(selectIndexVec[i],j))) )
        {
        // utlPrintVar1(true, j);
        bVector->push_back(bVec[i][j]);
        double gx = utl::ConvertStringToNumber<double>(gradStrMatrix[i][j][0]);
        double gy = utl::ConvertStringToNumber<double>(gradStrMatrix[i][j][1]);
        double gz = utl::ConvertStringToNumber<double>(gradStrMatrix[i][j][2]);
        double gnorm = std::sqrt(gx*gx+gy*gy+gz*gz);
        (*orientationsCartesian)(bVector->size()-1,0) = gx/gnorm;
        (*orientationsCartesian)(bVector->size()-1,1) = gy/gnorm;
        (*orientationsCartesian)(bVector->size()-1,2) = gz/gnorm;
        // utlPrintVar3(true, (*orientationsCartesian)(bVector->size()-1,0), (*orientationsCartesian)(bVector->size()-1,1), (*orientationsCartesian)(bVector->size()-1,2));
        }
      }

    // read individual dwi
    int numberOfDWI;
    if (m_IsInput4DImage)
      {
      itk::ReadImage<DWI4DImageType>(dataStr, dwi4DTempImage);
      typename DWI4DImageType::SizeType size = dwi4DTempImage->GetLargestPossibleRegion().GetSize();
      numberOfDWI = size[3];
      }
    else
      {
      itk::ReadImage<DWIImageType>(dataStr, dwiTempImage);
      numberOfDWI = dwiTempImage->GetNumberOfComponentsPerPixel();
      }

    // allocate b0, dwi
    if (i==0)
      {
      if (IsImageEmpty(m_B0Image))
        {
        m_B0Image = B0ImageType::New();
        if (m_IsInput4DImage)
          itk::CopyImageInformation<DWI4DImageType, B0ImageType>( dwi4DTempImage, m_B0Image);
        else
          itk::CopyImageInformation<DWIImageType, B0ImageType>( dwiTempImage, m_B0Image);
        m_B0Image->Allocate();
        if (numberOfB0>0)
          m_B0Image->FillBuffer(0);
        else 
          {
          std::cout << "b0=1 is used" << std::endl << std::flush;
          m_B0Image->FillBuffer(1.0);  // assume b0=1.0, if there is no b0 images.
          }
        }
      if (m_IsInput4DImage)
        itk::CopyImageInformation<DWI4DImageType, DWIImageType>( dwi4DTempImage, dwiImage);
      else
        itk::CopyImageInformation<DWIImageType, DWIImageType>( dwiTempImage, dwiImage);
      dwiImage->SetNumberOfComponentsPerPixel( numberOfDWIsWithoutB0 );
      dwiImage->Allocate();
      dwiImage->FillBuffer(dwiZeroPixel);
      }
    
    // concatenate dwi
    ImageRegionIteratorWithIndex<B0ImageType> b0It(m_B0Image, m_B0Image->GetLargestPossibleRegion());
    ImageRegionIteratorWithIndex<DWIImageType> dwiIt(dwiImage, dwiImage->GetLargestPossibleRegion());
    ImageRegionIteratorWithIndex<MaskImageType> maskIt;
    if (!IsImageEmpty(m_MaskImage))
      maskIt = ImageRegionIteratorWithIndex<MaskImageType>(m_MaskImage, m_MaskImage->GetLargestPossibleRegion());
    IndexType dwiIndex;

    if (m_IsInput4DImage)
      {
      typename DWI4DImageType::IndexType dwi4DTempIndex;
      for (b0It.GoToBegin(), dwiIt.GoToBegin(), maskIt.GoToBegin(); 
        !b0It.IsAtEnd(); 
        ++b0It, ++dwiIt, ++maskIt) 
        {
        if (!IsImageEmpty(m_MaskImage) && maskIt.Get()<=1e-10)
          {
          dwiIt.Set(dwiZeroPixel);
          b0It.Set(0.0);
          continue;
          }

        b0Index = b0It.GetIndex();
        if (this->GetDebug())
          std::cout << "index = " << b0Index << std::endl << std::flush;
        for ( int kk = 0; kk < VImageDimension; kk += 1 ) 
          dwi4DTempIndex[kk] = b0Index[kk];

        std::ostringstream ossWarn, ossWarnTmp;
        // additional b0
        double sumPixelb=0;
        for ( int kk = 0; kk < additionalB0[i]; kk += 1 ) 
          {
          dwi4DTempIndex[VImageDimension] = kk;
          double dwi4DTempPixel = dwi4DTempImage->GetPixel(dwi4DTempIndex);
          sumPixelb += dwi4DTempPixel*dwi4DTempPixel;
          if (dwi4DTempPixel<1e-10)
            ossWarnTmp << "zero b0 at "<< dwi4DTempIndex <<", " << "b0=" << dwi4DTempPixel << ", ";
          if (numberOfB0>0)
            {
            b0Pixel = b0It.Get() + dwi4DTempPixel/((double)numberOfB0);
            b0It.Set(b0Pixel);
            }
          }

        // b0 in gradients
        for ( int kk = 0; kk < b0IndexVec[i].size(); kk += 1 ) 
          {
          dwi4DTempIndex[VImageDimension] = b0IndexVec[i][kk]+additionalB0[i];
          double dwi4DTempPixel = dwi4DTempImage->GetPixel(dwi4DTempIndex);
          sumPixelb += dwi4DTempPixel*dwi4DTempPixel;
          if (dwi4DTempPixel<1e-10)
            ossWarnTmp << "zero b0 at "<< dwi4DTempIndex <<", " << "b0=" << dwi4DTempPixel << ", ";
          if (numberOfB0>0)
            {
            b0Pixel = b0It.Get() + dwi4DTempPixel/((double)numberOfB0);
            b0It.Set(b0Pixel);
            }
          }

        if (b0It.Get()<1e-10)
          {
          dwiIt.Set(dwiZeroPixel);
          continue;
          }

        // dwi without b0
        double sumPixeldwi=0;
        dwiPixel = dwiIt.Get();
        int jj=0;
        for ( int kk = additionalB0[i]; kk < numberOfDWI; kk += 1 ) 
          {
          if (!utl::IsInVector(b0IndexVec[i], kk-additionalB0[i]) && ( stringMatrix[i].size()==3 || (stringMatrix[i].size()==4 && utl::IsInVector(selectIndexVec[i],kk-additionalB0[i]))) )
            {
            dwi4DTempIndex[VImageDimension] = kk;
            double dwi4DTempPixel = dwi4DTempImage->GetPixel(dwi4DTempIndex);
            sumPixeldwi += dwi4DTempPixel*dwi4DTempPixel;
            if (dwi4DTempPixel<1e-10)
              ossWarnTmp << "zero dwi at "<< dwi4DTempIndex <<", ";
            dwiPixel[jj + (i==0?0:numDWINob0[i-1])] = dwi4DTempPixel;
            jj++;
            }
          }
        if (ossWarnTmp.str()!="" && sumPixelb>0 && sumPixeldwi>0)
          ossWarn << ossWarnTmp.str();

        if (m_ShowWarnings && (b0It.Get()>1e-10 || dwiPixel.GetSquaredNorm()>1e-10) && ossWarn.str()!="")
          {
          // utlWarning(true, "Warning: " << ossWarn.str() << "dwiPixel=" << dwiPixel<<";");
          std::cout << "Warning: " << ossWarn.str() << std::endl << std::flush;
          }

        if (this->GetDebug())
          {
          utlPrintVar(true, b0It.Get());
          itk::PrintVariableLengthVector(dwiPixel, "dwiPixel");
          }

        dwiIt.Set(dwiPixel);

        }

      }
    else
      {
      ImageRegionIteratorWithIndex<DWIImageType> dwiTempIt(dwiTempImage, dwiTempImage->GetLargestPossibleRegion());
      PixelType dwiTempPixel;
      for (b0It.GoToBegin(), dwiIt.GoToBegin(), dwiTempIt.GoToBegin(), maskIt.GoToBegin(); 
        !b0It.IsAtEnd(); 
        ++b0It, ++dwiIt, ++dwiTempIt, ++maskIt) 
        {
        if (!IsImageEmpty(m_MaskImage) && maskIt.Get()<=1e-10)
          {
          dwiIt.Set(dwiZeroPixel);
          b0It.Set(0.0);
          continue;
          }

        b0Index = b0It.GetIndex();
        if (this->GetDebug())
          std::cout << "index = " << b0Index << std::endl << std::flush;

        dwiTempPixel = dwiTempIt.Get();
        std::ostringstream ossWarn;
        // additional b0
        for ( int kk = 0; kk < additionalB0[i]; kk += 1 ) 
          {
          if (dwiTempPixel[kk]<1e-10)
            ossWarn << "zero b0 at ["<< b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<","<<kk <<"], b0=" << dwiTempPixel[kk] << "\n";
          if (numberOfB0>0)
            {
            b0Pixel = b0It.Get() + dwiTempPixel[kk]/((double)numberOfB0);
            b0It.Set(b0Pixel);
            }
          }
        // b0 in gradients
        for ( int kk = 0; kk < b0IndexVec[i].size(); kk += 1 ) 
          {
          if (dwiTempPixel[kk]<1e-10)
            ossWarn << "zero b0 at ["<< b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<","<<kk <<"], b0=" << dwiTempPixel[kk] << "\n";
          if (numberOfB0>0)
            {
            b0Pixel = b0It.Get() + dwiTempPixel[ b0IndexVec[i][kk]+additionalB0[i] ]/((double)numberOfB0);
            b0It.Set(b0Pixel);
            }
          }

        if (b0It.Get()<1e-10)
          {
          dwiIt.Set(dwiZeroPixel);
          continue;
          }

        // dwi without b0
        dwiPixel = dwiIt.Get();
        int jj=0;
        for ( int kk = additionalB0[i]; kk < numberOfDWI; kk += 1 ) 
          {
          if (!utl::IsInVector(b0IndexVec[i], kk-additionalB0[i]) && ( stringMatrix[i].size()==3 || (stringMatrix[i].size()==4 && utl::IsInVector(selectIndexVec[i],kk-additionalB0[i]))) )
            {
            if (dwiTempPixel[kk]<1e-10)
              ossWarn << "zero dwi at ["<< b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<","<<kk <<"], dwi value=" << dwiTempPixel[kk] << "\n";
            dwiPixel[jj + (i==0?0:numDWINob0[i-1])] = dwiTempPixel[kk];
            jj++;
            }
          }

        if (m_ShowWarnings &&  (b0It.Get()>1e-10 && dwiPixel.GetSquaredNorm()>1e-10) && ossWarn.str()!="")
          {
          // utlWarning(true,  "Warning: " << ossWarn.str() << "dwiPixel=" << dwiPixel<<";");
          std::cout << "Warning: " << ossWarn.str() << std::endl << std::flush;
          }

        if (this->GetDebug())
          {
          utlPrintVar(true, b0It.Get());
          itk::PrintVariableLengthVector(dwiPixel, "dwiPixel");
          }
        dwiIt.Set(dwiPixel);
        }

      }

    }

  m_SamplingSchemeQSpace->SetOrientationsCartesian(orientationsCartesian);
  m_SamplingSchemeQSpace->SetBVector(bVector);
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIReader<TPixelType, VImageDimension>
::CorrectDWI()
{
  utlShowPosition(this->GetDebug());
  utlGlobalException(!m_B0Image, "no b0 images");
  typename DWIImageType::Pointer dwiImage = this->GetOutput();

  utlGlobalException(!dwiImage, "no dwi images");
  utlGlobalException(!(itk::VerifyImageSize<DWIImageType, B0ImageType>(dwiImage, m_B0Image, true)), "dwi image and b0 image have different information");

  ImageRegionIteratorWithIndex<B0ImageType> b0It(m_B0Image, m_B0Image->GetLargestPossibleRegion());
  ImageRegionIteratorWithIndex<DWIImageType> dwiIt(dwiImage, dwiImage->GetLargestPossibleRegion());
  ImageRegionIteratorWithIndex<MaskImageType> maskIt;
  if (!IsImageEmpty(m_MaskImage))
    maskIt = ImageRegionIteratorWithIndex<MaskImageType>(this->m_MaskImage, m_MaskImage->GetLargestPossibleRegion());
  typename B0ImageType::IndexType  b0Index;
  typename B0ImageType::PixelType  b0Pixel;
  PixelType dwiPixel;
  int numberOfDWIsWithoutB0 = dwiImage->GetNumberOfComponentsPerPixel();
  dwiPixel.SetSize(numberOfDWIsWithoutB0);
  for (b0It.GoToBegin(), dwiIt.GoToBegin(), maskIt.GoToBegin(); 
    !b0It.IsAtEnd(); 
    ++b0It, ++dwiIt, ++maskIt) 
    {
    if (!IsImageEmpty(m_MaskImage) && maskIt.Get()<=1e-10)
      continue;
    b0Pixel = b0It.Get();
    if (b0Pixel<1e-10)
      {
      dwiPixel = dwiIt.Get();
      dwiPixel.Fill(0.0);
      dwiIt.Set(dwiPixel);
      continue;
      }

    dwiPixel = dwiIt.Get();
    if (dwiPixel.GetSquaredNorm()<1e-10)
      {
      b0It.Set(0.0);
      continue;
      }

    b0Index = b0It.GetIndex();
    double maxData=0.0, minData=std::numeric_limits<double>::max();
    for ( int j = 0; j < numberOfDWIsWithoutB0; j += 1 ) 
      {
      if (dwiPixel[j]<=0.99*b0Pixel && dwiPixel[j]>=1e-8*b0Pixel && dwiPixel[j]>maxData)
        maxData = dwiPixel[j];
      if (dwiPixel[j]<=0.99*b0Pixel && dwiPixel[j]>=1e-8*b0Pixel && dwiPixel[j]<minData)
        minData = dwiPixel[j];
      }
    
    if (maxData < minData)
      {
      if (this->GetDebug())
        {
        std::cout << "Logical ERROR: ("<<b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<"), all dwi values are not in [0, 0.99]*b0. " << std::endl << std::flush;
        b0It.Set(0.0);
        dwiPixel = dwiIt.Get();
        dwiPixel.Fill(0.0);
        dwiIt.Set(dwiPixel);
        continue;
        }
      }

    for ( int j = 0; j < numberOfDWIsWithoutB0; j += 1 ) 
      {
      if (dwiPixel[j]>=(1+1e-8)*b0Pixel)
        {
        if (this->GetDebug())
          {
          std::cout << "Logical ERROR in itk::DWIReader, dwi("<<b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<","<<j<<")=" << dwiPixel[j] << " > " << "b0("<<b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<")=" << b0Pixel
            << " Thus we force data["<<j<<"]="<< maxData <<", which is the maximal value of the other dwi values."<< std::endl;
          }
        dwiPixel[j] = maxData;
        }
      if (dwiPixel[j]<=1e-8*b0Pixel)
        {
        if (this->GetDebug())
          {
          std::cout << "Logical ERROR in itk::DWIReader, dwi("<<b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<","<<j<<")=" << dwiPixel[j] << " < " << "b0("<<b0Index[0]<<","<<b0Index[1]<<","<<b0Index[2]<<")=" << b0Pixel
            << " Thus we force data["<<j<<"]="<< minData <<", which is the minimal value of the other dwi values."<< std::endl;
          }
        dwiPixel[j] = minData;
        }
      }
    dwiIt.Set(dwiPixel);
    }

  this->Modified();
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIReader<TPixelType, VImageDimension>
::NormalizeDWI()
{
  utlShowPosition(this->GetDebug());
  utlGlobalException(!m_B0Image, "no b0 images");
  typename DWIImageType::Pointer dwiImage = this->GetOutput();

  utlGlobalException(!dwiImage, "no dwi images");
  utlGlobalException(!(itk::VerifyImageSize<DWIImageType, B0ImageType>(dwiImage, m_B0Image, true)), "dwi image and b0 image have different information");

  ImageRegionIterator<B0ImageType> b0It(m_B0Image, m_B0Image->GetLargestPossibleRegion());
  ImageRegionIterator<DWIImageType> dwiIt(dwiImage, dwiImage->GetLargestPossibleRegion());
  typename B0ImageType::PixelType  b0Pixel;
  PixelType dwiPixel;
  dwiPixel.SetSize(dwiImage->GetNumberOfComponentsPerPixel());

  for (b0It.GoToBegin(), dwiIt.GoToBegin(); 
    !b0It.IsAtEnd(); 
    ++b0It, ++dwiIt) 
    {
    b0Pixel = b0It.Get();
    if (b0Pixel<1e-10)
      continue;
    dwiPixel = dwiIt.Get();
    dwiPixel /= b0Pixel;
    dwiIt.Set(dwiPixel);
    }
  this->Modified();
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIReader<TPixelType, VImageDimension>
::GenerateData()
{
  if (m_ConfigurationFile!="")
    ReadFromConfigurationFile(m_ConfigurationFile);
  else
    utlGlobalException(true, "unknown format! ");
  
  if (m_SamplingSchemeQSpace->GetBThresholdSingleShell()>0)
    m_SamplingSchemeQSpace->CorrectBValues();

  // correct values if needed
  if (m_CorrectDWIValues)
    CorrectDWI();

  // normalize values if needed
  if (m_NormalizeDWI)
    NormalizeDWI();
}

template <class TPixelType, unsigned int VImageDimension>
void
DWIReader<TPixelType, VImageDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  PrintVar3(true, m_IsInput4DImage, m_NormalizeDWI, m_CorrectDWIValues, os<<indent );
  PrintVar1(true, m_ConfigurationFile, os<<indent);
  os << indent << "m_SamplingSchemeQSpace = " << m_SamplingSchemeQSpace << std::endl << std::flush;
  if (!IsImageEmpty(m_MaskImage))
    os << indent << "m_MaskImage = " << m_MaskImage << std::endl << std::flush;
  if (!IsImageEmpty(m_B0Image))
    os << indent << "m_B0Image = " << m_B0Image << std::endl << std::flush;
}

}

#endif 


