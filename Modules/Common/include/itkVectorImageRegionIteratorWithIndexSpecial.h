/**
 *       @file  itkVectorImageRegionIteratorWithIndexSpecial.h
 *      @brief  
 *     Created  "11-12-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkVectorImageRegionIteratorWithIndexSpecial_h
#define __itkVectorImageRegionIteratorWithIndexSpecial_h

#include "itkVectorImageRegionIteratorWithIndex.h"
#include "itkSpatiallyDenseSparseVectorImage.h"


namespace itk
{

template <typename TPixel, unsigned int VImageDimension>
class VectorImageRegionIteratorWithIndex<SpatiallyDenseSparseVectorImage<TPixel, VImageDimension> > : public ImageRegionIteratorWithIndex<SpatiallyDenseSparseVectorImage<TPixel, VImageDimension> >
{

public:
  typedef VectorImageRegionIteratorWithIndex                Self;
  typedef ImageRegionIteratorWithIndex< SpatiallyDenseSparseVectorImage<TPixel, VImageDimension>  > Superclass;

  /** Types inherited from the Superclass */
  typedef typename Superclass::IndexType             IndexType;
  typedef typename Superclass::SizeType              SizeType;
  typedef typename Superclass::OffsetType            OffsetType;
  typedef typename Superclass::RegionType            RegionType;
  typedef typename Superclass::ImageType             ImageType;
  typedef typename Superclass::PixelContainer        PixelContainer;
  typedef typename Superclass::PixelContainerPointer PixelContainerPointer;
  typedef typename Superclass::InternalPixelType     InternalPixelType;
  typedef typename Superclass::PixelType             PixelType;
  typedef typename Superclass::AccessorType          AccessorType;

  VectorImageRegionIteratorWithIndex() : Superclass()
  {}

  VectorImageRegionIteratorWithIndex(ImageType *ptr, const RegionType & region) : Superclass(ptr, region)
  {
  }

  VectorImageRegionIteratorWithIndex(const ImageIteratorWithIndex< ImageType > & it) : Superclass(it)
    {
    }

  void GetVector(PixelType& vec) const
  { 
  vec = Superclass::Get();
  }
  
  void SetVector(const PixelType & value)
  { 
  Superclass::Set(value);
  }
};

}

#endif 
