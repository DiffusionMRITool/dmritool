/**
 *       @file  itkUnaryFunctorVectorImageFilterTest.cxx
 *      @brief  
 *     Created  "09-12-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "itkUnaryFunctorVectorImageFilterTestCLP.h"
#include "itkUnaryFunctorVectorImageFilter.h"

#include "itkMultiVariableFunctorVectorImageFilter.h"

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
  OP_ABS,
  OP_EXP,
  OP_LOG,
  OP_SQUARE,
  OP_SQRT,
  // unary op (vector)
  OP_TNORM,
  OP_TMEAN,
  OP_TMEDIAN,
  OP_XMEDIAN,
  OP_YMEDIAN,
  OP_ZMEDIAN,
  OP_TSUM,
  OP_TMAX,
  OP_TMIN,
  OP_XMEAN,
  OP_YMEAN,
  OP_ZMEAN,
  OP_XNORM,
  OP_YNORM,
  OP_ZNORM,
  OP_XSUM,
  OP_YSUM,
  OP_ZSUM,
  // multiple op
  OP_TCOMPOSE,
  OP_XCOMPOSE,
  OP_YCOMPOSE,
  OP_ZCOMPOSE,
  // binary op (vector)
  OP_TDOTPRODUCT,
  OP_XDOTPRODUCT,
  OP_YDOTPRODUCT,
  OP_ZDOTPRODUCT,
  // complex op
  OP_TSHRED,
  OP_XSHRED,
  OP_YSHRED,
  OP_ZSHRED,
  OP_CROP,
};

int Operation = OP_NULL;
std::vector<int> size;
bool debug=false;
int numberOfThreads = -1;


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

#define __SetImage(ImageType, fileName, imageName, imageItName, axis)                                                       \
  typename ImageType::Pointer imageName = ImageType::New();                                                                 \
  itk::VectorImageRegionIteratorWithIndex<ImageType> imageItName;                                                           \
  if (fileName!="")                                                                                                         \
    {                                                                                                                       \
    itk::ReadImage<ImageType>(fileName, imageName);                                                                         \
    std::vector<int> sizeTmp = itk::GetVectorImage3DVolumeSize(imageName);                                                  \
    utlGlobalException(!utl::IsSameVector(size,sizeTmp), "images have different size");                                     \
    imageItName = itk::VectorImageRegionIteratorWithIndex<ImageType> (imageName, imageName->GetLargestPossibleRegion(),axis);    \
    }                                                                                                                   


#define __XYZTUnaryVectorFunctorWithArguments(argNoAxis, opNoAxis, funcName)                                                             \
  if (_X##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##X##opNoAxis);                                                                               \
    funcName.SetArguments(_X##argNoAxis);                                                                                                \
    outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, funcName,0);                                           \
    }                                                                                                                                    \
  if (_Y##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##Y##opNoAxis);                                                                               \
    funcName.SetArguments(_Y##argNoAxis);                                                                                                \
    outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, funcName,1);                                           \
    }                                                                                                                                    \
  if (_Z##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##Z##opNoAxis);                                                                               \
    funcName.SetArguments(_Z##argNoAxis);                                                                                                \
    outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, funcName,2);                                           \
    }                                                                                                                                    \
  if (_T##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##T##opNoAxis);                                                                               \
    funcName.SetArguments(_T##argNoAxis);                                                                                                \
    outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, funcName,3);                                           \
    }                                                                                                                                    \


#define __XYZTMultiVectorFunctor(argNoAxis, opNoAxis, funcName)                                                                          \
  if (_X##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##X##opNoAxis);                                                                               \
    outImage = MultiVectorOPImage<ImageType, ImageOutType>(image, _X##argNoAxis, _MaskImageFile, funcName,0);                            \
    }                                                                                                                                    \
  if (_Y##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##Y##opNoAxis);                                                                               \
    outImage = MultiVectorOPImage<ImageType, ImageOutType>(image, _Y##argNoAxis, _MaskImageFile, funcName,1);                            \
    }                                                                                                                                    \
  if (_Z##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##Z##opNoAxis);                                                                               \
    outImage = MultiVectorOPImage<ImageType, ImageOutType>(image, _Z##argNoAxis, _MaskImageFile, funcName,2);                            \
    }                                                                                                                                    \
  if (_T##argNoAxis##Arg.isSet())                                                                                                        \
    {                                                                                                                                    \
    SetOperationWithChecking(Operation, OP_##T##opNoAxis);                                                                               \
    outImage = MultiVectorOPImage<ImageType, ImageOutType>(image, _T##argNoAxis, _MaskImageFile, funcName,3);                            \
    }                                                                                                                                    \

template <class ImageType, class ImageOutType, class OpFunctor>
typename ImageOutType::Pointer
UnaryVectorOPImage(const itk::SmartPointer<ImageType>& image, std::string _MaskImageFile, const OpFunctor& func, int vectorAxis=3)
{
  utlSAAssert(vectorAxis>=0 && vectorAxis<=3)(vectorAxis).msg("wrong vectorAxis");
  itk::VectorImageRegionIteratorWithIndex<ImageType> it(image, image->GetLargestPossibleRegion(), vectorAxis);

  typedef itk::Image<double, 4> ScalarImageType;
  __SetImage(ScalarImageType, _MaskImageFile, mask, maskIt, vectorAxis);


  typedef itk::UnaryFunctorVectorImageFilter<ImageType, ImageOutType, OpFunctor, ScalarImageType> UnaryFunctorFilterType;
  typename UnaryFunctorFilterType::Pointer filter = UnaryFunctorFilterType::New();

  if (_MaskImageFile!="")
    filter->SetMaskImage(mask);
  filter->SetInput(image);
  filter->SetVectorAxis(vectorAxis);
  filter->SetFunctor(func);
  filter->SetDebug(debug);
  utlPrintVar(true, numberOfThreads);

  utl::InitializeThreadedLibraries(numberOfThreads);
  if (numberOfThreads>0)
    filter->SetNumberOfThreads(numberOfThreads);

  filter->Update();

  typename ImageOutType::Pointer outImage = filter->GetOutput();       


  return outImage;
}

template <class ImageType, class ImageOutType, class OpFunctor>
typename ImageOutType::Pointer
MultiVectorOPImage(const itk::SmartPointer<ImageType>& image, const std::vector<std::string>& inputStrVec, std::string _MaskImageFile, const OpFunctor& func, int vectorAxis=3)
{
  utlSAAssert(vectorAxis>=0 && vectorAxis<=3)(vectorAxis).msg("wrong vectorAxis");

  typedef itk::MultiVariableFunctorVectorImageFilter<ImageType, ImageOutType, OpFunctor> FunctorFilterType;
  typename FunctorFilterType::Pointer filter = FunctorFilterType::New();

  if (_MaskImageFile!="")
    filter->SetMaskImage(_MaskImageFile);
  filter->SetInput(image);
  typename ImageType::Pointer imageTmp = ImageType::New();
  for ( int i = 0; i < inputStrVec.size(); ++i ) 
    {
    itk::ReadImage(inputStrVec[i], imageTmp);
    filter->SetInput(i+1, imageTmp);
    }
  filter->SetVectorAxis(vectorAxis);
  filter->SetFunctor(func);
  filter->SetDebug(debug);
  if (numberOfThreads>0)
    filter->SetNumberOfThreads(numberOfThreads);

  std::cout << "p1" << std::endl << std::flush;
  filter->Update();
  std::cout << "p2" << std::endl << std::flush;

  typename ImageOutType::Pointer outImage = filter->GetOutput();       
  return outImage;
}

int 
main (int argc, char const* argv[])
{
  // GenerateCLP
  PARSE_ARGS;
  
  debug = _Debug;
  numberOfThreads = _NumberOfThreads; 


  utlGlobalException(!_OutputImageFileArg.isSet(), "need to set output -o");

  // typedef itk::VectorImage<double, 3> ImageType;
  typedef itk::Image<double, 4> ImageType;
  // typedef itk::VectorImage<double, 3> ImageOutType;
  typedef itk::Image<double, 4> ImageOutType;
  
  typename ImageType::Pointer image = ImageType::New();
  typename ImageOutType::Pointer outImage  = ImageOutType::New();       

  itk::ReadImage<ImageType>(_InputImageFile, image);
  size = itk::GetVectorImage3DVolumeSize(image);
  int dimT = itk::GetVectorImageVectorSize(image);

    {
    utl::Functor::Compose<utl::Vector<double> > op;
    __XYZTMultiVectorFunctor(Compose, COMPOSE, op);
    }

  // unary
  // if (_SquareArg.isSet())
  //   {
  //   SetOperationWithChecking(Operation, OP_SQUARE);
  //   outImage = UnaryOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::Square<double>());
  //   }

  // unary vector
  if (_TNormArg.isSet())
    {
    SetOperationWithChecking(Operation, OP_TNORM);
    if (_TNorm=="L2") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::TwoNorm<utl::Vector<double> >());
    if (_TNorm=="L1") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::OneNorm<utl::Vector<double> >());
    if (_TNorm=="L0") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::ZeroNorm<utl::Vector<double> >());
    if (_TNorm=="INF") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::InfNorm<utl::Vector<double> >());
    }
  if (_XNormArg.isSet())
    {
    SetOperationWithChecking(Operation, OP_XNORM);
    // if (_XNorm=="L2") __XYZVectorToScalarFunctor(_XNorm, OP_XNORM, "X", utl::Functor::TwoNorm<utl::Vector<double> >());
    if (_XNorm=="L2") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::TwoNorm<utl::Vector<double> >(), 0);
    if (_XNorm=="L1") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::OneNorm<utl::Vector<double> >(), 0);
    if (_XNorm=="L0") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::ZeroNorm<utl::Vector<double> >(), 0);;
    if (_XNorm=="INF") outImage = UnaryVectorOPImage<ImageType, ImageOutType>(image, _MaskImageFile, utl::Functor::InfNorm<utl::Vector<double> >(), 0);
    }

    {
    utl::Functor::Shred<utl::Vector<double> > op;
    __XYZTUnaryVectorFunctorWithArguments(Shred, SHRED, op);
    }


  itk::SaveImage(outImage, _OutputImageFile);
  
  return 0;
}
