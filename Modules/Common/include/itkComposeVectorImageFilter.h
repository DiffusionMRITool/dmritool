/**
 *       @file  itkComposeVectorImageFilter.h
 *      @brief  
 *     Created  "06-27-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkComposeVectorImageFilter_h
#define __itkComposeVectorImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorImage.h"
#include "itkImageRegionConstIterator.h"
#include <vector>
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "utlCoreMacro.h"

namespace itk
{
/** \class ComposeVectorImageFilter
 * \brief ComposeVectorImageFilter combine several vector images into a vector image
 *
 * \par Inputs and Usage
 * \code
 *    filter->SetInput( 0, image0 );
 *    filter->SetInput( 1, image1 );
 *    ...
 *    filter->Update();
 *    itk::VectorImage< PixelType, dimension >::Pointer = filter->GetOutput();
 * \endcode
 * All input images are expected to have the same template parameters and have
 * the same size and origin.
 *
 * \sa VectorImage
 * \sa VectorIndexSelectionCastImageFilter
 * \ingroup ITKImageCompose
 *
 * \wiki
 * \wikiexample{VectorImages/ImageToVectorImageFilter,Create a vector image from a collection of scalar images}
 * \wikiexample{ImageProcessing/Compose3DCovariantVectorImageFilter,Compose a vector image (with 3 components) from three scalar images}
 * \wikiexample{SpectralAnalysis/RealAndImaginaryToComplexImageFilter,Convert a real image and an imaginary image to a complex image}
 * \endwiki
 */

template< typename TInputImage=VectorImage<double,3>, typename TOutputImage=VectorImage<typename TInputImage::PixelType, TInputImage::ImageDimension> >
class ComposeVectorImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  typedef ComposeVectorImageFilter                         Self;
  typedef SmartPointer< Self >                             Pointer;
  typedef SmartPointer< const Self >                       ConstPointer;
  typedef ImageToImageFilter< TInputImage, TOutputImage >  Superclass;
  itkNewMacro(Self);
  itkTypeMacro(ComposeVectorImageFilter, ImageToImageFilter);

  itkStaticConstMacro(Dimension, unsigned int, TInputImage::ImageDimension);

  typedef TInputImage                          InputImageType;
  typedef TOutputImage                         OutputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename InputImageType::RegionType  RegionType;

  void SetInput1(const InputImageType *image1);
  void SetInput2(const InputImageType *image2);
  void SetInput3(const InputImageType *image3);


protected:
  ComposeVectorImageFilter();

  virtual void GenerateOutputInformation(void) ITK_OVERRIDE;

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

  virtual void ThreadedGenerateData(const RegionType & outputRegionForThread, ThreadIdType) ITK_OVERRIDE;

private:
  ComposeVectorImageFilter(const Self &);
  void operator=(const Self &);


  // we have to specialize the code for complex, because it provides no operator[]
  // method
  typedef ImageRegionConstIterator< InputImageType > InputIteratorType;
  typedef std::vector< InputIteratorType >           InputIteratorContainerType;

  // template<typename T>
  // void ComputeOutputPixel(std::complex<T> & pix, InputIteratorContainerType & inputItContainer )
  //   {
  //   pix = std::complex<T>(inputItContainer[0].Get(), inputItContainer[1].Get());
  //   ++( inputItContainer[0] );
  //   ++( inputItContainer[1] );
  //   }
  template<typename TPixel>
  void ComputeOutputPixel(TPixel & pix, InputIteratorContainerType & inputItContainer)
    {
    int j=0;
    for ( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
      {
      typename InputImageType::PixelType pixel = inputItContainer[i].Get();
      for ( int k = 0; k < pixel.Size(); ++k ) 
        {
        pix[j] = pixel[k];
        j++;
        }
      ++( inputItContainer[i] );
      }
    }
};


// template <class ImageType, class FilterType >
// SmartPointer<ImageType>
// VectorImageChannelFilterApply ( const SmartPointer<ImageType>& input, const SmartPointer<FilterType>& filter )
// {
//   typedef typename FilterType::InputImageType FilterInputImageType;
//   typedef typename FilterType::OutputImageType FilterOutputImageType;
//   typedef itk::ComposeVectorImageFilter<FilterOutputImageType, ImageType> ComposeImageFilterType;
//   typename ComposeImageFilterType::Pointer composeImageFilter = ComposeImageFilterType::New();

//   typedef itk::VectorIndexSelectionCastImageFilter<ImageType, FilterInputImageType> IndexSelectionType;

//   int N = input->GetNumberOfComponentsPerPixel();
//   input->Print(std::cout<<"input=\n");
//   for ( int i = 0; i < N; ++i ) 
//     {
//     typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
//     indexSelectionFilter->SetIndex(i);
//     indexSelectionFilter->SetInput(input);
//     indexSelectionFilter->Update();
    
//     typename FilterType::Pointer filterCopy = filter->Clone();
//     filterCopy->SetInput(indexSelectionFilter->GetOutput());
//     filterCopy->Update();
//     typename FilterOutputImageType::Pointer out = filterCopy->GetOutput();

//     utlPrintVar3(true, i, N, out->GetNumberOfComponentsPerPixel());

//     composeImageFilter->SetInput(i, out);
//     }

//   composeImageFilter->Update();
//   return composeImageFilter->GetOutput();
// }

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkComposeVectorImageFilter_hxx)
#include "itkComposeVectorImageFilter.hxx"
#endif

#endif

