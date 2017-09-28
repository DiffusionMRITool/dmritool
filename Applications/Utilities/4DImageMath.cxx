/**
 *       @file  4DImageMath.cxx
 *      @brief  
 *     Created  "07-12-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "4DImageMathCLP.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkPermuteAxesImageFilter.h"

#include "itkVectorImageRegionIterator.h"
#include "itkVectorImageRegionIteratorWithIndex.h"
#include "utl.h"
#include "itkFunctors.h"
#include "itkMultiVolumeImageToVectorImageFilter.h"
#include "itkVectorImageToMultiVolumeImageFilter.h"

#include "itkUnaryFunctorVectorImageFilter.h"
#include "itkMultiVariableFunctorVectorImageFilter.h"

#include "itkFunctorFromStringImageFilter.h"
  
typedef double  ScalarType; 
typedef itk::Image<ScalarType, 4> NDImageType;
typedef itk::VectorImage<ScalarType, 3> VectorImageType;

template <class IndexType>
bool 
IsOutsideBox ( const IndexType& index, const std::vector<int>& box )
{
  return index[0]<box[0] || index[0]>box[1] || index[1]<box[2] || index[1]>box[3] || index[2]<box[4] || index[2]>box[5];
}


enum
{
  // no op
  OP_NULL=0,
  // binary op (scalar)
  OP_ADD,
  OP_MINUS,
  OP_MULTIPLY,
  OP_DIVIDE,
  OP_MAX,
  OP_MIN,
  // unary op (scalar)
  OP_FUNC,
  OP_ABS,
  OP_EXP,
  OP_LOG,
  OP_SQUARE,
  OP_SQRT,
  // unary op (vector)
  OP_NORM,
  OP_MEAN,
  OP_MEDIAN,
  OP_SUM,
  OP_MAX_AXIS,
  OP_MIN_AXIS,
  // binary op (vector)
  OP_DOTPRODUCT,
  // multiple op
  OP_COMPOSE,
  // complex op
  OP_SHRED,
  OP_CROP,
};
  
int _Operation = OP_NULL;
std::vector<int> _realBox;
std::vector<int> _size;
int _numberOfThreads = -1;
std::string _axis;
int _axisNumber=-1;

void
SetOperationWithChecking( int &op, int value )
{
  if ( op==OP_NULL )
    op = value;
  else if ( op==value )
    ;
  else
    utlGlobalException(true, "Only one type of operation is allowed!");
}

#define _SetOperationWithChecking(op, value, numberOfInputsCond)                                   \
do                                                                                                 \
  {                                                                                                \
  SetOperationWithChecking(op, value);                                                             \
  utlSAGlobalException(!(numberOfInputsCond))(numberOfInputsCond).msg("wrong number of inputs!");  \
  } while ( 0 );


void 
GetImageFiles(const std::vector<std::string>& imageVec, std::string& inImage0, std::string& outImage)
{
  utlGlobalException(imageVec.size()<2, "need to set at least 2 image files (at least one input image, and one output image)");
  inImage0 = imageVec[0];
  outImage = imageVec.back();
}

int GetNumberFromAxis(std::string axis)
{
  if (axis=="X" || axis=="x" || axis=="0") return 0;
  else if (axis=="Y" || axis=="y" || axis=="1" ) return 1;
  else if (axis=="Z" || axis=="z" || axis=="2" ) return 2;
  else if (axis=="T" || axis=="t" || axis=="3" ) return 3;
  else
    utlGlobalException(true, "wrong logic");
}

template <class ImageType, class Image2Type>
void 
ConvertImage ( const itk::SmartPointer<ImageType>& image, itk::SmartPointer<Image2Type>& imageOut )
{
  imageOut = image;
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


// template <class ImageType>
// void
// PermuteImage ( const itk::SmartPointer<ImageType>& image, const std::string axis, itk::SmartPointer<NDImageType>& imageOut )
// {
//   imageOut = NDImageType::New();
//   typedef itk::PermuteAxesImageFilter <NDImageType> PermuteAxesImageFilterType;
//   PermuteAxesImageFilterType::Pointer permuteAxesFilter = PermuteAxesImageFilterType::New();
//   itk::FixedArray<unsigned int, 4> order;
//   if (axis=="X")
//     order[0]=3, order[1]=1, order[2]=2, order[3]=0;
//   else if (axis=="Y")
//     order[0]=0, order[1]=3, order[2]=2, order[3]=1;
//   else if (axis=="Z")
//     order[0]=0, order[1]=1, order[2]=3, order[3]=2;
//   else
//     utlGlobalException(true, "axis must be X or Y or Z");
//   permuteAxesFilter->SetOrder(order);
//   auto imageTmp = NDImageType::New();
//   ConvertImage(image, imageTmp);
//   permuteAxesFilter->SetInput(imageTmp);
//   permuteAxesFilter->Update();
//   imageOut = permuteAxesFilter->GetOutput();
// }

#define __SetImage(ImageType, fileName, imageName, imageItName, axis)                                                       \
  typename ImageType::Pointer imageName = ImageType::New();                                                                 \
  itk::VectorImageRegionIteratorWithIndex<ImageType> imageItName;                                                           \
  if (fileName!="")                                                                                                         \
    {                                                                                                                       \
    itk::ReadImage<ImageType>(fileName, imageName);                                                                         \
    std::vector<int> sizeTmp = itk::GetVectorImage3DVolumeSize(imageName);                                                  \
    utlGlobalException(!utl::IsSameVector(_size,sizeTmp), "images have different size");                                     \
    imageItName = itk::VectorImageRegionIteratorWithIndex<ImageType> (imageName, imageName->GetLargestPossibleRegion(),axis);    \
    }                                                                                                                   


#define __SetImageOrScalar(ImageType, fileName, imageName, imageItName, valueName)                                          \
  typename ImageType::Pointer imageName = ImageType::New();                                                                 \
  itk::VectorImageRegionIteratorWithIndex<ImageType> imageItName;                                                           \
  double valueName=0;                                                                                                       \
  if (fileName!="")                                                                                                         \
    {                                                                                                                       \
    if (utl::IsNumber(fileName))                                                                                            \
      {                                                                                                                     \
      valueName = utl::ConvertStringToNumber<double>(fileName);                                                             \
      }                                                                                                                     \
    else                                                                                                                    \
      {                                                                                                                     \
      itk::ReadImage<ImageType>(fileName, imageName);                                                                       \
      std::vector<int> sizeTmp = itk::GetVectorImage3DVolumeSize(imageName);                                                \
      utlGlobalException(!utl::IsSameVector(_size,sizeTmp), "images have different size");                                   \
      imageItName = itk::VectorImageRegionIteratorWithIndex<ImageType> (imageName, imageName->GetLargestPossibleRegion());  \
      }                                                                                                                     \
    }                                                                                                                   


// #define __XYZVectorToScalarFunctor(argName, opName, axisStr, funcName)                                                                 \
//   do {                                                                                                                                 \
//   if (argName##Arg.isSet())                                                                                                            \
//     {                                                                                                                                  \
//     typename NDImageType::Pointer outImagePermute = NDImageType::New();                                                                \
//     typename NDImageType::Pointer outImagePermute2 = NDImageType::New();                                                               \
//     SetOperationWithChecking(_Operation, opName);                                                                                       \
//     PermuteImage(image, axisStr, imagePermute);                                                                                        \
//     UnaryVectorOPImage<NDImageType, NDImageType>(imagePermute, outImagePermute, _MaskImageFile, funcName);                           \
//     PermuteImage(outImagePermute, axisStr, outImagePermute2);                                                                          \
//     ConvertImage(outImagePermute2, outImage);                                                                                          \
//     }                                                                                                                                  \
//   } while(0)



#define __FunctionFromStringOpImage(funcStr)                                                                                              \
do                                                                                                                                        \
  {                                                                                                                                       \
  itk::FunctorFromStringOPImage(imageVec, outImage, funcStr, mask, _numberOfThreads);                                                     \
  } while ( 0 );


#define __UnaryScalarFunctor(argNoAxis, opNoAxis, funcName)                                                                              \
  if (argNoAxis##Arg.isSet())                                                                                                            \
    {                                                                                                                                    \
    _SetOperationWithChecking(_Operation, OP_##opNoAxis, imageVec.size()==1);                                                            \
    utl::Functor::ScalarFunctorWrapper<funcName > func;                                                                                  \
    itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, func, mask, _numberOfThreads);                                     \
    }


#define __XYZTUnaryVectorFunctor(axisNum, argNoAxis, opNoAxis, funcName)                                                                 \
  if (argNoAxis##Arg.isSet())                                                                                                            \
  {                                                                                                                                      \
  _SetOperationWithChecking(_Operation, OP_##opNoAxis, imageVec.size()==1);                                                              \
  funcName func;                                                                                                                         \
  itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, func, mask, _numberOfThreads, axisNum);                              \
  }


#define __XYZTUnaryVectorFunctorWithArguments(axisNum, argNoAxis, opNoAxis, funcName)                                                    \
  if (argNoAxis##Arg.isSet())                                                                                                            \
    {                                                                                                                                    \
    _SetOperationWithChecking(_Operation, OP_##opNoAxis, imageVec.size()==1);                                                            \
    funcName func;                                                                                                                       \
    func.SetArguments(argNoAxis);                                                                                                        \
    itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, func, mask, _numberOfThreads, axisNum);                            \
    }                                                                                                                                   


#define __XYZTMultiVectorFunctor(axisNum, argNoAxis, opNoAxis, funcName)                                                                 \
  if (argNoAxis##Arg.isSet())                                                                                                            \
    {                                                                                                                                    \
    _SetOperationWithChecking(_Operation, OP_##opNoAxis, imageVec.size()>=2);                                                            \
    funcName func;                                                                                                                       \
    itk::MultiVariableVectorOPImage<ImageType, ImageOutType>(imageVec, outImage, func, mask, _numberOfThreads, axisNum);                 \
    }                                                                                                                                    


#define __XYZTMultiVectorFunctor_FixedSize(axisNum, argNoAxis, opNoAxis, funcName, inSize)                                               \
  if (argNoAxis##Arg.isSet())                                                                                                            \
    {                                                                                                                                    \
    _SetOperationWithChecking(_Operation, OP_##opNoAxis, imageVec.size()==inSize);                                                       \
    funcName func;                                                                                                                       \
    itk::MultiVariableVectorOPImage<ImageType, ImageOutType>(imageVec, outImage, func, mask, _numberOfThreads, axisNum);                 \
    }                                                                                                                                    


// #define __XYZTBinaryVectorFunctor(axisNum, argNoAxis, opNoAxis, funcName)                                                                \
//   if (argNoAxis##Arg.isSet())                                                                                                            \
//     {                                                                                                                                    \
//     _SetOperationWithChecking(_Operation, OP_##opNoAxis, imageVec.size()==2);                                                                                 \
//     utlGlobalException(imageVec.size()!=2, "need 2 input images");                                                                       \
//     BinaryVectorOPImage<ImageType, ImageOutType>(image, imageVec[1], outImage, _MaskImageFile, funcName, axisNum);                       \
//     }                                                                                                                                    

// [>*  Unary operations. Output image has the same size as the input image. <]  
// template <class ImageType, class ImageOutType, class OpFunctor>
// void
// UnaryOPImage(const itk::SmartPointer<ImageType>& image, itk::SmartPointer<ImageOutType>& outImage, std::string _MaskImageFile, const OpFunctor& func)
// {
//   itk::VectorImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion());

//   typedef itk::Image<double, 4> ScalarImageType;
//   __SetImage(ScalarImageType, _MaskImageFile, mask, maskIt, -1);
  
//   outImage  = ImageOutType::New();       
//   itk::CopyImageInformation(image, outImage);
//   outImage->Allocate();
//   itk::VectorImageRegionIteratorWithIndex<ImageOutType> outIt(outImage, outImage->GetLargestPossibleRegion());   

//   typename ImageType::IndexType index;
//   itk::VariableLengthVector<double> inPixel, outPixel;
//   outPixel.SetSize(itk::GetVectorImageVectorSize(image));
//   for (it.GoToBegin(), maskIt.GoToBegin(), outIt.GoToBegin(); 
//     !it.IsAtEnd(); 
//     ++it, ++maskIt, ++outIt)
//     {
//     if (_MaskImageFile!="" && maskIt.Get()<=0)
//       { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }

//     index = it.GetIndex();
//     if (IsOutsideBox(index, _realBox))
//       { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }

//     it.GetVector(inPixel);

//     for ( int i = 0; i < inPixel.GetSize(); ++i ) 
//       outPixel[i] = func(inPixel[i]);

//     if (utl::IsLogDebug())
//       {
//       std::cout << "index = " << index << std::endl << std::flush;
//       itk::PrintVariableLengthVector(inPixel, "inPixel");
//       itk::PrintVariableLengthVector(outPixel, "outPixel");
//       }
//     outIt.SetVector(outPixel);
//     }
// }

/**  Binary operations. Output image has the same size as the input image. */  
template <class ImageType, class ImageOutType, class OpFunctor>
void
BinaryOPImage(const itk::SmartPointer<ImageType>& image, itk::SmartPointer<ImageOutType>& outImage, std::string _MaskImageFile, std::string _opImageFile, const OpFunctor& func)
{
  utlGlobalException(_opImageFile=="", "cannot be empty");
  bool isValue = utl::IsNumber(_opImageFile);

  itk::VectorImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion());

  typedef itk::Image<double, 4> ScalarImageType;
  __SetImage(ScalarImageType, _MaskImageFile, mask, maskIt, -1);

  __SetImageOrScalar(ImageType, _opImageFile, opImage, opIt, opValue);

  int imageVecSize = itk::GetVectorImageVectorSize(image);
  if (!isValue)
    {
    int opImageVecSize = itk::GetVectorImageVectorSize(opImage);
    utlSAGlobalException(opImageVecSize!=1 && opImageVecSize!=imageVecSize)
      (opImageVecSize)(imageVecSize).msg("wrong size of inputImage and opImage");
    }
  
  outImage  = ImageOutType::New();       
  itk::CopyImageInformation(image, outImage);
  outImage->Allocate();
  itk::VectorImageRegionIteratorWithIndex<ImageOutType> outIt(outImage, outImage->GetLargestPossibleRegion());   

  typename ImageType::IndexType index;
  itk::VariableLengthVector<double> inPixel, outPixel, inTmp;
  outPixel.SetSize(imageVecSize);
  for (it.GoToBegin(), maskIt.GoToBegin(), opIt.GoToBegin(), outIt.GoToBegin(); 
    !it.IsAtEnd(); 
    ++it, ++maskIt, ++opIt, ++outIt)
    {
    if (_MaskImageFile!="" && maskIt.Get()<=0)
      { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }

    index = it.GetIndex();
    if (IsOutsideBox(index, _realBox))
      { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }

    it.GetVector(inPixel);

    if (isValue)
      {
      for ( int i = 0; i < inPixel.GetSize(); ++i ) 
        outPixel[i] = func(inPixel[i], opValue);
      }
    else
      {
      opIt.GetVector(inTmp);
      for ( int i = 0; i < inPixel.GetSize(); ++i ) 
        outPixel[i] = inTmp.Size()==1 ? func(inPixel[i], inTmp[0]) : func(inPixel[i], inTmp[i]);
      }

    if (utl::IsLogDebug())
      {
      std::cout << "index = " << index << std::endl << std::flush;
      itk::PrintVariableLengthVector(inPixel, "inPixel");
      itk::PrintVariableLengthVector(outPixel, "outPixel");
      }
    outIt.SetVector(outPixel);
    }

}

// template <class ImageType, class ImageOutType, class OpFunctor>
// void 
// BinaryVectorOPImage(const itk::SmartPointer<ImageType>& image, const itk::SmartPointer<ImageType>& image2, itk::SmartPointer<ImageOutType>& outImage, std::string _MaskImageFile, const OpFunctor& func, int vectorAxis=3)
// {
//   utlGlobalException(!itk::VerifyImageSize(image, image2, false), "the two images should have the same shape.");
//   utlSAAssert(vectorAxis>=0 && vectorAxis<=3)(vectorAxis).msg("wrong vectorAxis");
//   itk::VectorImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion(), vectorAxis);
//   itk::VectorImageRegionIteratorWithIndex<ImageType> it2(image2, image2->GetLargestPossibleRegion(), vectorAxis);

//   typedef itk::Image<double, 4> ScalarImageType;
//   __SetImage(ScalarImageType, _MaskImageFile, mask, maskIt, vectorAxis);

//   int vecInputSize = itk::GetVectorImageVectorSize(image);
  
//   outImage  = ImageOutType::New();       
//   itk::CopyImageInformation(image, outImage);
//   int outVectorSize = func.GetOutputDimension(vecInputSize);
//   int outDim= (vectorAxis==3) ? func.GetOutputDimension(vecInputSize) : vecInputSize;
//   if (vectorAxis!=3)
//     {
//     typename ImageOutType::RegionType regionTmp = outImage->GetLargestPossibleRegion();
//     typename ImageOutType::SizeType sizeTmp = regionTmp.GetSize();
//     sizeTmp[vectorAxis]=func.GetOutputDimension(sizeTmp[vectorAxis]);
//     outVectorSize = func.GetOutputDimension(sizeTmp[vectorAxis]);
//     regionTmp.SetSize(sizeTmp);
//     outImage->SetRegions(regionTmp);
//     }
//   utlPrintVar(utl::IsLogDebug(), vectorAxis, outDim, vecInputSize, func.GetOutputDimension(vecInputSize));
//   utlSAException(outDim<=0)(outDim).msg("wrong outDim");
//   SetVectorImageVectorSize(outImage, outDim);
//   outImage->Allocate();
//   // outImage->Print(std::cout<<"outImage =\n");
//   itk::VectorImageRegionIteratorWithIndex<ImageOutType> outIt(outImage, outImage->GetLargestPossibleRegion(), vectorAxis);   

//   typename ImageType::IndexType index;
//   itk::VariableLengthVector<double> inPixel, inPixel2, outPixel;
//   utl::Vector<double> inVec, inVec2, outVec(outVectorSize);
//   outPixel.SetSize(outVectorSize);
//   outPixel.Fill(0.0);
//   for (it.GoToBegin(), it2.GoToBegin(), maskIt.GoToBegin(), outIt.GoToBegin(); 
//     !it.IsAtEnd(); 
//     ++it, ++it2, ++maskIt, ++outIt)
//     {
//     if (_MaskImageFile!="" && maskIt.Get()<=0)
//       { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }

//     index = it.GetIndex();
//     if (IsOutsideBox(index, _realBox))
//       { outPixel.Fill(0.0); outIt.SetVector(outPixel); continue; }

//     for ( int i = 0; i < (vectorAxis==3?1:vecInputSize); ++i ) 
//       {
//       it.GetVector(inPixel, i);
//       it2.GetVector(inPixel2, i);

//       inVec = utl::VariableLengthVectorToUtlVector(inPixel);
//       inVec2 = utl::VariableLengthVectorToUtlVector(inPixel2);
//       outVec = func(inVec, inVec2);
//       outPixel = utl::UtlVectorToVariableLengthVector(outVec);

//       if (utl::IsLogDebug())
//         {
//         std::cout << "index = " << index << ", i = " << i << std::endl << std::flush;
//         itk::PrintVariableLengthVector(inPixel, "inPixel");
//         itk::PrintVariableLengthVector(inPixel2, "inPixel2");
//         itk::PrintVariableLengthVector(outPixel, "outPixel");
//         }
//       outIt.SetVector(outPixel,i);
//       }
//     }
// }

template <class ImageType, class ImageOutType>
int
ImageMath(int argc, char const* argv[])
{
  PARSE_ARGS;
  utl::LogLevel = _Verbose;
  _numberOfThreads = _NumberOfThreads; 

  if (_AxisArg.isSet())
    {
    _axis = utl::StringToUpperCase(_Axis);
    _axisNumber = GetNumberFromAxis(_axis);
    }
  
  std::string _InputImageFile, _OutputImageFile;
  GetImageFiles(_ImageFiles, _InputImageFile, _OutputImageFile);
    
  typedef itk::Image<double,4> Scalar4DImageType;
  typename Scalar4DImageType::Pointer mask = Scalar4DImageType::New();
  if (_MaskImageFileArg.isSet())
    {
    itk::ReadImage(_MaskImageFile, mask);
    }

  typedef typename ImageType::Pointer ImagePointer;
  typename ImageType::Pointer image = ImageType::New();
  typename NDImageType::Pointer imagePermute = NDImageType::New();
  typename ImageOutType::Pointer outImage  = ImageOutType::New();       
      
  typedef utl::Functor::VectorUnaryFunctionWrapper<> VectorUnaryFunctionWrapper;
  typedef utl::Functor::VectorMultiVariableFunctionWrapper<>  VectorMultiVariableFunctionWrapper;

  std::vector<ImagePointer> imageVec(_ImageFiles.size()-1);
  for ( int i = 0; i < _ImageFiles.size()-1; ++i ) 
    {
    itk::ReadImage<ImageType>(_ImageFiles[i], imageVec[i]);
    }
  image = imageVec[0];
  _size = itk::GetVectorImage3DVolumeSize(image);
  int dimT = itk::GetVectorImageVectorSize(image);

  utlGlobalException(_Box.size()!=6, "wrong size of _Box");
  _realBox = _Box;
  _realBox[0] = _realBox[0]<0 ? 0 : (_realBox[0]>_size[0]-1  ? _size[0]-1  : _realBox[0]); 
  _realBox[1] = _realBox[1]==-1 ? _size[0]-1 : (_realBox[1]<0 ? 0 : (_realBox[1]>_size[0]-1  ? _size[0]-1  : _realBox[1]) ); 
  _realBox[2] = _realBox[2]<0 ? 0 : (_realBox[2]>_size[1]-1  ? _size[1]-1  : _realBox[2]); 
  _realBox[3] = _realBox[3]==-1 ? _size[1]-1 : (_realBox[3]<0 ? 0 : (_realBox[3]>_size[1]-1  ? _size[1]-1  : _realBox[3]) ); 
  _realBox[4] = _realBox[4]<0 ? 0 : (_realBox[4]>_size[2]-1  ? _size[2]-1  : _realBox[4]); 
  _realBox[5] = _realBox[5]==-1 ? _size[2]-1 : (_realBox[5]<0 ? 0 : (_realBox[5]>_size[2]-1  ? _size[2]-1  : _realBox[5]) ); 

  // binary 
  if (_AddImageFileArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_ADD, imageVec.size()==1);
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, _AddImageFile, std::plus<double>());
    }
  if (_MinusImageFileArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_MINUS, imageVec.size()==1);
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, _MinusImageFile, std::minus<double>());
    }
  if (_MultiplyImageFileArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_MULTIPLY, imageVec.size()==1);
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, _MultiplyImageFile, std::multiplies<double>());
    }
  if (_DivideImageFileArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_DIVIDE, imageVec.size()==1);
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, _DivideImageFile, std::divides<double>());
    }
  if (_MaxImageFileArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_MAX, imageVec.size()==1);
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, _MaxImageFile, utl::Functor::Max<double>());
    }
  if (_MinImageFileArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_MIN, imageVec.size()==1);
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, _MinImageFile, utl::Functor::Min<double>());
    }


  // math expression
  if (_FunctorArg.isSet())
    {
    SetOperationWithChecking(_Operation, OP_FUNC);
    utlPrintVar(true, _Functor);
     __FunctionFromStringOpImage(_Functor);
    }

  // unary
  __UnaryScalarFunctor(_Abs, ABS, utl::Functor::Abs<double>);
  __UnaryScalarFunctor(_Exp, EXP, utl::Functor::Exp<double>);
  __UnaryScalarFunctor(_Log, LOG, utl::Functor::Log<double>);
  __UnaryScalarFunctor(_Square, SQUARE, utl::Functor::Square<double>);
  __UnaryScalarFunctor(_Sqrt, SQRT, utl::Functor::Sqrt<double>);

  // unary vector
  if (_NormArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_NORM, imageVec.size()==1);
    if (_Norm=="L2") itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, utl::Functor::TwoNorm<utl::Vector<double> >(), mask, _numberOfThreads, _axisNumber);
    if (_Norm=="L1") itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, utl::Functor::OneNorm<utl::Vector<double> >(), mask, _numberOfThreads, _axisNumber);
    if (_Norm=="L0") itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, utl::Functor::ZeroNorm<utl::Vector<double> >(), mask, _numberOfThreads, _axisNumber);
    if (_Norm=="INF") itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, outImage, utl::Functor::InfNorm<utl::Vector<double> >(), mask, _numberOfThreads, _axisNumber);
    }
  __XYZTUnaryVectorFunctor(_axisNumber, _Mean, MEAN, utl::Functor::Mean<utl::Vector<double> >);
  __XYZTUnaryVectorFunctor(_axisNumber, _Median, MEDIAN, utl::Functor::Median<utl::Vector<double> >);
  __XYZTUnaryVectorFunctor(_axisNumber, _Sum, SUM, utl::Functor::Sum<utl::Vector<double> >);
  __XYZTUnaryVectorFunctor(_axisNumber, _Max_Axis, MAX_AXIS, utl::Functor::MaxValue<utl::Vector<double> >);
  __XYZTUnaryVectorFunctor(_axisNumber, _Min_Axis, MIN_AXIS, utl::Functor::MinValue<utl::Vector<double> >);

  __XYZTUnaryVectorFunctorWithArguments(_axisNumber, _Shred, SHRED, utl::Functor::Shred<utl::Vector<double> >);

  __XYZTMultiVectorFunctor(_axisNumber, _Compose, COMPOSE, utl::Functor::Compose<utl::Vector<double> >);
    
  // __XYZTBinaryVectorFunctor(_axisNumber, _DotProduct, DOTPRODUCT, utl::Functor::DotProduct<utl::Vector<double> >());
  __XYZTMultiVectorFunctor_FixedSize(_axisNumber, _DotProduct, DOTPRODUCT, utl::Functor::DotProduct<utl::Vector<double> >, 2);

  if (_CropArg.isSet())
    {
    _SetOperationWithChecking(_Operation, OP_CROP, imageVec.size()==1);
    std::vector<int> cropBox(8,-1);
    utlGlobalException(_Crop.size()>8, "wrong size of _Crop, should be no more than 8");
    for ( int i = 0; i < _Crop.size(); ++i ) 
      cropBox[i] = _Crop[i];
    cropBox[0] = cropBox[0]<0 ? 0 : (cropBox[0]>_size[0]-1  ? _size[0]-1  : cropBox[0]); 
    cropBox[1] = cropBox[1]==-1 ? _size[0]-1 : (cropBox[1]<0 ? 0 : (cropBox[1]>_size[0]-1  ? _size[0]-1  : cropBox[1]) ); 
    cropBox[2] = cropBox[2]<0 ? 0 : (cropBox[2]>_size[1]-1  ? _size[1]-1  : cropBox[2]); 
    cropBox[3] = cropBox[3]==-1 ? _size[1]-1 : (cropBox[3]<0 ? 0 : (cropBox[3]>_size[1]-1  ? _size[1]-1  : cropBox[3]) ); 
    cropBox[4] = cropBox[4]<0 ? 0 : (cropBox[4]>_size[2]-1  ? _size[2]-1  : cropBox[4]); 
    cropBox[5] = cropBox[5]==-1 ? _size[2]-1 : (cropBox[5]<0 ? 0 : (cropBox[5]>_size[2]-1  ? _size[2]-1  : cropBox[5]) ); 
    cropBox[6] = cropBox[6]<0 ? 0 : (cropBox[6]>dimT-1  ? dimT-1  : cropBox[6]); 
    cropBox[7] = cropBox[7]==-1 ? dimT-1 : (cropBox[7]<0 ? 0 : (cropBox[7]>dimT-1  ? dimT-1  : cropBox[7]) ); 

    utlPrintVar(true, dimT);
    utl::PrintVector(_Crop, "_Crop");
    utl::PrintVector(cropBox, "cropBox");
    if (cropBox[6]==0 && cropBox[7]==dimT-1)
      {
      // no crop in t-axis
      typename ImageType::Pointer imageTmp = ImageType::New();
      typename ImageType::IndexType desiredStart;
      typename ImageType::SizeType desiredSize;

      for ( int i = 0; i < ImageType::ImageDimension; ++i ) 
        desiredStart[i] = cropBox[2*i];
      for ( int i = 0; i < ImageType::ImageDimension; ++i ) 
        desiredSize[i] = cropBox[2*i+1]-cropBox[2*i]+1;

      typename ImageType::RegionType desiredRegion(desiredStart, desiredSize);

      typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetRegionOfInterest(desiredRegion);
      filter->SetInput(image);
      filter->Update();
      imageTmp = filter->GetOutput();
      ConvertImage(imageTmp, outImage);
      }
    else
      {
      // crop in t-axis
      utl::Functor::Shred<utl::Vector<double> > op;
      op.SetOffset(cropBox[6]);
      op.SetChunkSize(cropBox[7]-cropBox[6]+1);
      op.SetSpace(dimT);
      typename ImageOutType::Pointer imageTmp = ImageOutType::New();
      itk::UnaryVectorOPImage<ImageType, ImageOutType>(image, imageTmp, op, mask, _numberOfThreads);

      if (cropBox[0]==0 && cropBox[1]==_size[0]-1 && cropBox[2]==0 && cropBox[3]==_size[1]-1 && cropBox[4]==0 && cropBox[5]==_size[2]-1)
        {
        // no crop in x/y/z axis
        outImage = imageTmp;
        }
      else
        {
        // crop in x/y/z axis
        typename ImageOutType::IndexType desiredStart;
        typename ImageOutType::SizeType desiredSize;

        for ( int i = 0; i < ImageOutType::ImageDimension; ++i ) 
          desiredStart[i] = cropBox[2*i];
        for ( int i = 0; i < ImageOutType::ImageDimension; ++i ) 
          desiredSize[i] = cropBox[2*i+1]-cropBox[2*i]+1;

        // correct t-axis, since crop in t-axis has been done already.
        if (ImageOutType::ImageDimension==4)
          {
          desiredStart[3] = 0;
          desiredSize[3] = cropBox[7]-cropBox[6]+1;
          }

        typename ImageOutType::RegionType desiredRegion(desiredStart, desiredSize);

        typedef itk::RegionOfInterestImageFilter< ImageOutType, ImageOutType > FilterType;
        typename FilterType::Pointer filter = FilterType::New();
        filter->SetRegionOfInterest(desiredRegion);
        filter->SetInput(imageTmp);
        filter->Update();
        outImage = filter->GetOutput();
        }
      }
    }

  if (_Operation!=OP_NULL)
    itk::SaveImage(outImage, _OutputImageFile);
  else if (_ImageFiles.size()==2)
    {
    // if no operation is set, then conversion of image formats (with mask if provided)
    BinaryOPImage<ImageType, ImageOutType>(image, outImage, _MaskImageFile, utl::ConvertNumberToString(1.0), std::multiplies<double>());
    // __FunctionFromStringOpImage("x");
    itk::SaveImage(outImage, _OutputImageFile);
    }
  else
    utlGlobalException(true, "wrong operator or inputs");

  return 0;
}

/**
 * \brief  4DImage math. It works for both itk::Image<double,4> and itk::VectorImage<double,3>
 */
int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;

  // utlGlobalException(!_OutputImageFileArg.isSet(), "need to set output -o");
  std::string _InputImageFile, _OutputImageFile;
  GetImageFiles(_ImageFiles, _InputImageFile, _OutputImageFile);
    
  if (!_OutputFormatArg.isSet() || (_OutputFormatArg.isSet() && _OutputFormat=="NONE"))
    {
    if (itk::IsVectorImage(_InputImageFile) || itk::Is3DImage(_InputImageFile))
      return ImageMath<VectorImageType, VectorImageType>(argc, argv);
    else
      return ImageMath<NDImageType, NDImageType>(argc, argv);
    }
  else
    {
    if (itk::IsVectorImage(_InputImageFile) && _OutputFormat=="VECTOR")
      return ImageMath<VectorImageType, VectorImageType>(argc, argv);
    else if (itk::IsVectorImage(_InputImageFile) && _OutputFormat=="4D")
      return ImageMath<VectorImageType, NDImageType>(argc, argv);
    else if (!itk::IsVectorImage(_InputImageFile) && _OutputFormat=="VECTOR")
      return ImageMath<NDImageType, VectorImageType>(argc, argv);
    else if (!itk::IsVectorImage(_InputImageFile) && _OutputFormat=="4D")
      return ImageMath<NDImageType, NDImageType>(argc, argv);
    }
  
  return 0;
}

