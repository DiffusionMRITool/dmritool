/**
 *       @file  PrintImage.cxx
 *      @brief  
 *     Created  "10-11-2013
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "PrintImageCLP.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorImageRegionIterator.h"
#include "itkVectorImageRegionIteratorWithIndex.h"
#include "itkVectorImageRegionIteratorWithIndexSpecial.h"
#include "itkSpatiallyDenseSparseVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImageFileReader.h"
#include "itkSpatiallyDenseSparseVectorImageFileWriter.h"
#include "itkVariableLengthVectorImageFileReader.h"
#include "itkVariableLengthVectorImageFileWriter.h"
#include "utl.h"

#include "itkMultiVolumeImageToVectorImageFilter.h"
#include "itkVectorImageToMultiVolumeImageFilter.h"

typedef double  ScalarType; 
typedef itk::VectorImage<ScalarType, 3> VectorImageType;
typedef itk::Image<ScalarType, 4> NDImageType;
typedef itk::SpatiallyDenseSparseVectorImage<ScalarType, 3> SparseImageType;
typedef itk::Image<itk::VariableLengthVector<ScalarType>, 3> VariableLengthImageType;

template <class ImageType, class Image2Type>
void 
ConvertImage ( const itk::SmartPointer<ImageType>& image, itk::SmartPointer<Image2Type>& imageOut )
{
  // if (itk::GetImageType(image)==itk::GetImageType(imageOut))
  //   imageOut = image;
}

template<>
void
ConvertImage<NDImageType, VectorImageType>( const itk::SmartPointer<NDImageType>& image, itk::SmartPointer<VectorImageType>& imageOut )
{
  itk::MultiVolumeToVectorImage(image, imageOut);
}

template<>
void
ConvertImage<VectorImageType, NDImageType>( const itk::SmartPointer<VectorImageType>& image, itk::SmartPointer<NDImageType>& imageOut )
{
  itk::VectorToMultiVolumeImage(image, imageOut);
}

  
template <class IndexType>
bool 
IsOutsideBox ( const IndexType& index, const std::vector<int>& box )
{
  return index[0]<box[0] || index[0]>box[1] || index[1]<box[2] || index[1]>box[3] || index[2]<box[4] || index[2]>box[5];
}

template <class ImageType, class ReaderType, class Image0Type=ImageType, class Reader0Type=ReaderType>
int
PrintImage(int argc, char const* argv[])
{
  PARSE_ARGS;
  typename ImageType::Pointer image = ImageType::New();
  itk::ReadImage<ImageType, ReaderType>(_InputImageFile, image);
  
  itk::VectorImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion());
  int vecSize = itk::GetVectorImageVectorSize(image);
  std::vector<int> size = itk::GetVectorImage3DVolumeSize(image);

  utlGlobalException(_Box.size()!=6, "wrong size of _Box");
  std::vector<int> realBox(_Box);
  realBox[0] = realBox[0]<0 ? 0 : (realBox[0]>size[0]-1  ? size[0]-1  : realBox[0]); 
  realBox[1] = realBox[1]==-1 ? size[0]-1 : (realBox[1]<0 ? 0 : (realBox[1]>size[0]-1  ? size[0]-1  : realBox[1]) ); 
  realBox[2] = realBox[2]<0 ? 0 : (realBox[2]>size[1]-1  ? size[1]-1  : realBox[2]); 
  realBox[3] = realBox[3]==-1 ? size[1]-1 : (realBox[3]<0 ? 0 : (realBox[3]>size[1]-1  ? size[1]-1  : realBox[3]) ); 
  realBox[4] = realBox[4]<0 ? 0 : (realBox[4]>size[2]-1  ? size[2]-1  : realBox[4]); 
  realBox[5] = realBox[5]==-1 ? size[2]-1 : (realBox[5]<0 ? 0 : (realBox[5]>size[2]-1  ? size[2]-1  : realBox[5]) ); 

  // utl::PrintVector(realBox, "realBox");

  typedef itk::Image<double, 3> ScalarImageType;
  ScalarImageType::Pointer mask = ScalarImageType::New();
  itk::ImageRegionIteratorWithIndex<ScalarImageType> maskIt;
  if (_MaskImageFile!="")
    {
    itk::ReadImage<ScalarImageType>(_MaskImageFile, mask);
    std::vector<int> sizeMask = itk::GetVectorImage3DVolumeSize(mask);
    utlGlobalException(!utl::IsSameVector(size,sizeMask), "mask image has different size");
    
    maskIt = itk::ImageRegionIteratorWithIndex<ScalarImageType> (mask, mask->GetLargestPossibleRegion());
    }

  utlGlobalException((_DifferencePercent==1 || _DifferencePercent==2) && _BaseImageFile=="", "need to set _BaseImageFile");

  typename ImageType::IndexType index;
  itk::VariableLengthVector<double> vec;

  std::cout << "\nImage option:" << std::endl;
  for (it.GoToBegin(), maskIt.GoToBegin(); 
    !it.IsAtEnd(); 
    ++it, ++maskIt)
    {
    if (_MaskImageFile!="" && maskIt.Get()<=0)
      continue;

    index = it.GetIndex();
    if (IsOutsideBox(index, realBox))
      continue;

    it.GetVector(vec);

    if (!_PrintAllVoxels && vec.GetNorm()<1e-8) 
      continue;
      
    if (_DifferencePercent!=2)
      {
      std::ostringstream oss;
      oss << "Image1 " << index << " ";
      
      if (_OrientationsFile=="" && std::fabs(_Power-1.0)>1e-10)
        {
        // only use power without _OrientationsFile
        for ( int i = 0; i < vec.Size(); ++i ) 
          vec[i] = std::pow(vec[i], _Power);
        }

      utl::PrintVector(vec, vec.Size(), oss.str());
      }
    }
  
  if (_BaseImageFile!="")
    {
    std::cout << "\nImage difference option:" << std::endl;

    typename Image0Type::Pointer image0 = Image0Type::New();
    int baseImageType = itk::GetImageType(_BaseImageFile);
    if (itk::GetImageType(image0)==baseImageType)
      {
      itk::ReadImage<Image0Type, Reader0Type>(_BaseImageFile, image0);
      }
    else
      {
      if (baseImageType==IMAGE_ND)
        {
        typename NDImageType::Pointer imageTmp = NDImageType::New();
        itk::ReadImage<NDImageType>(_BaseImageFile, imageTmp);
        ConvertImage(imageTmp, image0);
        }
      else if (baseImageType==IMAGE_VECTOR)
        {
        typename VectorImageType::Pointer imageTmp = VectorImageType::New();
        itk::ReadImage<VectorImageType>(_BaseImageFile, imageTmp);
        ConvertImage(imageTmp, image0);
        }
      else
        utlGlobalException(true, "wrong image type");
      }

    std::vector<int> size0 = itk::GetVectorImage3DVolumeSize(image0);
    int vecSize0 = itk::GetVectorImageVectorSize(image0);
    utlGlobalException(!utl::IsSameVector(size0,size) || vecSize!=vecSize0, "these two image have different sizes");

    itk::VectorImageRegionIteratorWithIndex<Image0Type> it0(image0, image0->GetLargestPossibleRegion());
    itk::VariableLengthVector<double> vec0, diff;

    for (it.GoToBegin(), it0.GoToBegin(), maskIt.GoToBegin(); 
      !it.IsAtEnd(); 
      ++it, ++it0, ++maskIt)
      {
      if (_MaskImageFile!="" && maskIt.Get()<=0)
        continue;

      index = it.GetIndex();
      if (IsOutsideBox(index, realBox))
        continue;
      
      it.GetVector(vec);
      it0.GetVector(vec0);

      if (_OrientationsFile=="" && std::fabs(_Power-1.0)>1e-10)
        {
        // only use power without _OrientationsFile
        for ( int i = 0; i < vec.Size(); ++i ) 
          {
          vec[i] = std::pow(vec[i],_Power);
          vec0[i] = std::pow(vec0[i],_Power);
          }
        }
      diff = vec-vec0;
      for ( int jj = 0; jj < diff.Size(); ++jj ) 
        diff[jj] = std::fabs(diff[jj]);

      if (!_PrintAllVoxels && (diff.GetNorm()<1e-8 || vec0.GetNorm()<1e-8 || vec.GetNorm()<1e-8)) 
        continue;

      if (_DifferencePercent!=2)
        {
        std::ostringstream oss;
        oss << "difference " << index << " ";

        utl::PrintVector(diff, diff.Size(), oss.str());
        }
      }

    if (_DifferencePercent==1 || _DifferencePercent==2)
      {
      std::cout << "\nImage difference percentage option:" << std::endl;

      std::vector<double> perVec;
      for (it.GoToBegin(), it0.GoToBegin(), maskIt.GoToBegin(); 
        !it.IsAtEnd(); 
        ++it, ++it0, ++maskIt)
        {
        if (_MaskImageFile!="" && maskIt.Get()<=0)
          continue;
          
        index = it.GetIndex();
        if (IsOutsideBox(index, realBox))
          continue;

        it.GetVector(vec);
        it0.GetVector(vec0);

        if (_OrientationsFile=="" && std::fabs(_Power-1.0)>1e-10)
          {
          for ( int i = 0; i < vec.Size(); ++i ) 
            {
            vec[i] = std::pow(vec[i],_Power);
            vec0[i] = std::pow(vec0[i],_Power);
            }
          }
        diff = vec-vec0;
        for ( int jj = 0; jj < diff.Size(); ++jj ) 
          diff[jj] = std::fabs(diff[jj]);

        if (!_PrintAllVoxels && (diff.GetNorm()<1e-8 || vec0.GetNorm()<1e-8 || vec.GetNorm()<1e-8)) 
          continue;
        
        double diffPer = diff.GetNorm()/vec0.GetNorm();
        if (_DifferencePercent==1)
          {
          std::ostringstream oss;
          oss << "difference percentage " << index << " ";
          std::cout << oss.str() << ": " << diffPer << std::endl;
          }
        perVec.push_back(diffPer);
        }
      utl::PrintVector(perVec, "percentage vector");
      }

    }

  if (vecSize==6)
    {
    std::cout << "\nDTI option:" << std::endl;
    itk::VariableLengthVector<double> dt(9);
    for (it.GoToBegin(), maskIt.GoToBegin(); 
      !it.IsAtEnd(); 
      ++it, ++maskIt)
      {
      if (_MaskImageFile!="" && maskIt.Get()<=0)
        continue;
      
      index = it.GetIndex();
      if (IsOutsideBox(index, realBox))
        continue;
      
      it.GetVector(vec);
      
      if (!_PrintAllVoxels && vec.GetNorm()<1e-8) 
        continue;

      if (_DifferencePercent!=2)
        {
        std::ostringstream oss;
        oss << "DTI " << index << " ";
        if (_TensorStorageFormat=="6D_UPPER") utl::ConvertTensor6DTo9D(vec, dt, TENSOR_UPPER_TRIANGULAR);
        if (_TensorStorageFormat=="6D_LOWER") utl::ConvertTensor6DTo9D(vec, dt, TENSOR_LOWER_TRIANGULAR);
        if (_TensorStorageFormat=="6D_EMBED") utl::ConvertTensor6DTo9D(vec, dt, TENSOR_EMBED6D);
        if (_TensorStorageFormat=="6D_DIAGONAL_FIRST") utl::ConvertTensor6DTo9D(vec, dt, TENSOR_DIAGONAL_FIRST);
        if (_TensorStorageFormat=="9D") dt=vec;
        utl::PrintVector(dt, dt.Size(), oss.str());
        }
      }
    }

  if (_OrientationsFile!="")
    {
    std::cout << "\nSpherical function option:" << std::endl;
    typedef utl::NDArray<double,2> MatrixType;
    typedef utl_shared_ptr<MatrixType> MatrixPointer;
    typedef utl::NDArray<double,1> VectorType;

    int shRank = utl::DimToRankSH(vecSize);
    utlPrintVar2(true, vecSize, shRank);
    MatrixPointer grad = utl::ReadGrad<double>(_OrientationsFile, DIRECTION_NODUPLICATE, CARTESIAN_TO_SPHERICAL);
    MatrixPointer basisMatrix = utl::ComputeSHMatrix(shRank, *grad, SPHERICAL_TO_SPHERICAL);

    VectorType sf, sh(vecSize);
    for (it.GoToBegin(), maskIt.GoToBegin(); 
      !it.IsAtEnd(); 
      ++it, ++maskIt)
      {
      if (_MaskImageFile!="" && maskIt.Get()<=0)
        continue;
      
      index = it.GetIndex();
      if (IsOutsideBox(index, realBox))
        continue;
      
      it.GetVector(vec);
      
      if (!_PrintAllVoxels && vec.GetNorm()<1e-8) 
        continue;

      utl::VectorToVector(vec, sh, vecSize);
      std::ostringstream oss;
      oss << "SF " << index << " ";
      sf = (*basisMatrix) * sh;

      if (std::fabs(_Power-1.0)>1e-10)
        {
        for ( int i = 0; i < sf.Size(); ++i ) 
          sf[i] = std::pow(sf[i], _Power);
        }
      utl::PrintVector(sf, sf.Size(), oss.str());
      }
    }
  return 0;
}

/**
 * \brief  Print image
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;


  int inputType = itk::GetImageType(_InputImageFile);

  if (inputType==IMAGE_VECTOR)
    return PrintImage<VectorImageType, itk::ImageFileReader<VectorImageType>>(argc, argv);
  else if (inputType==IMAGE_ND)
    return PrintImage<NDImageType, itk::ImageFileReader<NDImageType>>(argc, argv);
  else if (inputType==IMAGE_SPARSE)
    return PrintImage<SparseImageType, itk::SpatiallyDenseSparseVectorImageFileReader<SparseImageType>>(argc, argv);
  else if (inputType==IMAGE_VARIABLELENGTH)
    return PrintImage<VariableLengthImageType, itk::VariableLengthVectorImageFileReader<VariableLengthImageType>>(argc, argv);

  // if (!_BaseImageFileArg.isSet())
  //   {
  //   if (isInputVectorImage)
  //     return PrintImage<VectorImageType, VectorImageType>(argc, argv);
  //   else 
  //     return PrintImage<NDImageType, VectorImageType>(argc, argv);
  //   }
  // else
  //   {
  //   isBaseVectorImage = itk::IsVectorImage(_BaseImageFile);
  //   if (isInputVectorImage && isBaseVectorImage)
  //     return PrintImage<VectorImageType, VectorImageType>(argc, argv);
  //   else if (!isInputVectorImage && isBaseVectorImage)
  //     return PrintImage<NDImageType, VectorImageType>(argc, argv);
  //   else if (isInputVectorImage && !isBaseVectorImage)
  //     return PrintImage<VectorImageType, NDImageType>(argc, argv);
  //   else if (!isInputVectorImage && !isBaseVectorImage)
  //     return PrintImage<NDImageType, NDImageType>(argc, argv);
  //   }
  
  return 0;
}
