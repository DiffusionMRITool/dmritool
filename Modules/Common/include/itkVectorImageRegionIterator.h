/**
 *       @file  itkVectorImageRegionIterator.h
 *      @brief  
 *     Created  "06-10-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef itkVectorImageRegionIterator_h
#define itkVectorImageRegionIterator_h

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "utlITKConceptChecking.h"
#include "itkVectorImage.h"
#include "itkVariableLengthVector.h"

namespace itk
{

/**
 * \class VectorImageRegionIterator
 * \brief A multi-dimensional iterator templated over image type.
 *
 * \ingroup ImageIterators
 * \ingroup ITKCommon
 *
 * \sa itkVectorImageRegionIteratorWithIndexTest.cxx 
 */

template< typename TImage >
class VectorImageRegionIterator : public ImageRegionIterator< TImage >
{
public:
  /** Standard class typedefs. */
  typedef VectorImageRegionIterator Self;

  /** Dimension of the image the iterator walks.  This constant is needed so
   * functions that are templated over image iterator type (as opposed to
   * being templated over pixel type and dimension) can have compile time
   * access to the dimension of the image that the iterator walks. */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int,
                      TImage::ImageDimension);

  /** Define the superclass */
  typedef ImageRegionIterator< TImage > Superclass;
  

  /** Inherit types from the superclass */
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

  typedef VariableLengthVector< InternalPixelType >  PixelVectorType;
  
#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( SameTypeCheck,
                   ( itk::Concept::SameType< ImageType, VectorImage<InternalPixelType, ImageIteratorDimension> > ) );
  itkConceptMacro( SameTypeCheck2,
                   ( itk::Concept::SameType< PixelType, PixelVectorType> ) );
#endif


  /** Default Constructor. Need to provide a default constructor since we
   * provide a copy constructor. */
  VectorImageRegionIterator() : Superclass() 
  {
  m_VectorStride=0;
  m_VectorSize=0;
  m_VectorAxis=ImageIteratorDimension;
  }

  /** Default Destructor */
  ~VectorImageRegionIterator() {}

  /** Copy Constructor. The copy constructor is provided to make sure the
   * handle to the image is properly reference counted. */
  VectorImageRegionIterator(const Self & it) : Superclass(it) 
  { 
  m_VectorStride = it.m_VectorStride;
  m_VectorSize = it.m_VectorSize;
  m_VectorAxis = it.m_VectorAxis;
  }

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  VectorImageRegionIterator(TImage *ptr, const RegionType & region, int vectorAxis=-1) : Superclass(ptr, region) 
  { 
  // -1 means t axis
  m_VectorAxis = (vectorAxis==-1) ? ImageIteratorDimension : vectorAxis;
  itkAssertOrThrowMacro((m_VectorAxis<=(int)ImageIteratorDimension && m_VectorAxis>=0), "wrong vector axis, m_VectorAxis=" << m_VectorAxis << ", ImageIteratorDimension=" << ImageIteratorDimension);

  if (m_VectorAxis==ImageIteratorDimension)
    {
    // Superclass(ptr, regionInput);
    m_VectorSize = ptr->GetNumberOfComponentsPerPixel();
    m_VectorStride = 1;
    }
  else
    {
    this->SetRegion(region);
    typename RegionType::SizeType regionSize = this->m_Region.GetSize();
    typename RegionType::IndexType regionIndex = this->m_Region.GetIndex();
    if (ptr)
      {
      SizeType size = ptr->GetLargestPossibleRegion().GetSize();
      const typename ImageType::OffsetValueType* offsetTable = ptr->GetOffsetTable();
      m_VectorStride = ptr->GetNumberOfComponentsPerPixel()*offsetTable[m_VectorAxis];
      m_VectorSize = size[m_VectorAxis];
      }
    }
  }

  /** operator= is provided to make sure the handle to the image is properly
   * reference counted. */
  Self & operator=(const Self & it) 
    {
    Superclass::operator=(it);
    m_VectorStride = it.m_VectorStride;
    m_VectorSize = it.m_VectorSize;
    m_VectorAxis = it.m_VectorAxis;
    };
  
  /** Set the region of the image to iterate over. */
  void SetRegion(const RegionType & regionInput)
  {
  if (m_VectorAxis==ImageIteratorDimension)
    Superclass::SetRegion(regionInput);
  else
    {
    typename RegionType::SizeType sizeInput = regionInput.GetSize();
    typename RegionType::IndexType indexInput = regionInput.GetIndex();
    sizeInput[m_VectorAxis] = 1;
    indexInput[m_VectorAxis] = 0;
    RegionType region;
    this->m_Region.SetIndex(indexInput);
    this->m_Region.SetSize(sizeInput);

    if ( region.GetNumberOfPixels() > 0 ) // If region is non-empty
      {
      const RegionType & bufferedRegion = this->m_Image->GetBufferedRegion();
      itkAssertOrThrowMacro( ( bufferedRegion.IsInside(this->m_Region) ),
                             "Region " << this->m_Region << " is outside of buffered region " << bufferedRegion );
      }

    // Compute the start offset
    this->m_Offset = this->m_Image->ComputeOffset( this->m_Region.GetIndex() );
    this->m_BeginOffset = this->m_Offset;

    // Compute the end offset. If any component of m_Region.GetSize()
    // is zero, the region is not valid and we set the EndOffset
    // to be same as BeginOffset so that iterator end condition is met
    // immediately.
    if ( this->m_Region.GetNumberOfPixels() == 0 )
      {
      // region is empty, probably has a size of 0 along one dimension
      this->m_EndOffset = this->m_BeginOffset;
      }
    else
      {
      IndexType ind( this->m_Region.GetIndex() );
      SizeType  size( this->m_Region.GetSize() );
      for ( unsigned int i = 0; i < ImageIteratorDimension; ++i )
        {
        if (i!=m_VectorAxis)
          ind[i] += ( static_cast< IndexValueType >( size[i] ) - 1 );
        else
          ind[i] = 0;
        }
      this->m_EndOffset = this->m_Image->ComputeOffset(ind);
      this->m_EndOffset++;
      }
    }
  }

  /** Set the pixel value */
  void SetVector(const PixelVectorType & value, const int offIndex=0) const
  {
  if (m_VectorAxis==ImageIteratorDimension)
    {
    this->m_PixelAccessorFunctor.Set(*( const_cast< InternalPixelType * >(
          this->m_Buffer + this->m_Offset ) ), value);
    }
  else
    {
    int numberOfComponens = this->m_Image->GetNumberOfComponentsPerPixel();
    InternalPixelType* p = const_cast<InternalPixelType*>(this->m_Buffer + numberOfComponens*this->m_Offset);
    int off = offIndex;
    for ( int i = 0; i < m_VectorSize; i += 1 ) 
      {
      *(p+off) = value[i];
      off += m_VectorStride;
      }
    }
  }
  
  /** Get the pixel value */
  void GetVector(PixelVectorType& vec, const int offIndex=0) const
  { 
  if (m_VectorAxis==ImageIteratorDimension)
    {
    vec = this->m_PixelAccessorFunctor.Get( *( this->m_Buffer + this->m_Offset ) );
    }
  else
    {
    vec.SetSize(m_VectorSize);
    int numberOfComponens = this->m_Image->GetNumberOfComponentsPerPixel();
    InternalPixelType* p = const_cast<InternalPixelType*>(this->m_Buffer + numberOfComponens*this->m_Offset);
    int off = offIndex;
    for ( int i = 0; i < m_VectorSize; i += 1 ) 
      {
      vec[i] = *(p+off);
      off += m_VectorStride;
      }
    }
  }
  


protected:

  /** This constructor is declared protected in order to enforce
    const-correctness */
  VectorImageRegionIterator(const ImageRegionIterator< TImage > & it);
  Self & operator=(const ImageRegionIterator< TImage > & it);

protected: 
  
  OffsetValueType m_VectorStride;
  int m_VectorSize;
  int m_VectorAxis;
};


template <typename TPixel, unsigned int VImageDimension>
class VectorImageRegionIterator<Image<TPixel, VImageDimension> > : public ImageRegionIterator<Image<TPixel, VImageDimension> >
{
public:
  /** Standard class typedefs. */
  typedef VectorImageRegionIterator Self;


  /** Dimension of the image the iterator walks.  This constant is needed so
   * functions that are templated over image iterator type (as opposed to
   * being templated over pixel type and dimension) can have compile time
   * access to the dimension of the image that the iterator walks. */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int,
                      VImageDimension);

  /** Define the superclass */
  typedef ImageRegionIterator< Image<TPixel, VImageDimension> > Superclass;

  /** Inherit types from the superclass */
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
  
  typedef VariableLengthVector< InternalPixelType >  PixelVectorType;

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( SameTypeCheck,
                   ( itk::Concept::SameType< ImageType, Image<InternalPixelType, VImageDimension> > ) );
#endif

  /** Default Constructor. Need to provide a default constructor since we
   * provide a copy constructor. */
  VectorImageRegionIterator() : Superclass() 
  {
  m_VolumeSize=0;
  m_VectorStride=0;
  m_VectorSize=0;
  m_VectorAxis=ImageIteratorDimension-1;
  }

  /** Default Destructor */
  ~VectorImageRegionIterator() {}

  /** Copy Constructor. The copy constructor is provided to make sure the
   * handle to the image is properly reference counted. */
  VectorImageRegionIterator(const Self & it) : Superclass(it) 
  { 
  m_VolumeSize = it.m_VolumeSize;
  m_VectorStride = it.m_VectorStride;
  m_VectorSize = it.m_VectorSize;
  m_VectorAxis = it.m_VectorAxis;
  };

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  VectorImageRegionIterator(ImageType *ptr, const RegionType & region, int vectorAxis=-1) : Superclass(ptr, region) 
  { 
  m_VectorAxis = (vectorAxis==-1) ? Superclass::ImageDimension-1 : vectorAxis;
  itkAssertOrThrowMacro((m_VectorAxis<(int)ImageIteratorDimension && m_VectorAxis>=0), "wrong vector axis, m_VectorAxis=" << m_VectorAxis);

  this->SetRegion(region);
  typename RegionType::SizeType regionSize = this->m_Region.GetSize();
  typename RegionType::IndexType regionIndex = this->m_Region.GetIndex();
  if (ptr)
    {
    SizeType size = ptr->GetLargestPossibleRegion().GetSize();
    const typename ImageType::OffsetValueType* offsetTable = ptr->GetOffsetTable();
    m_VectorStride = offsetTable[m_VectorAxis];
    m_VectorSize = size[m_VectorAxis];
    }
  } 

  /** operator= is provided to make sure the handle to the image is properly
   * reference counted. */
  Self & operator=(const Self & it) 
    {
    Superclass::operator=(it);
    m_VolumeSize = it.m_VolumeSize;
    m_VectorStride = it.m_VectorStride;
    m_VectorSize = it.m_VectorSize;
    m_VectorAxis = it.m_VectorAxis;
    }

  /** Set the region of the image to iterate over. */
  void SetRegion(const RegionType & regionInput)
  {
    typename RegionType::SizeType sizeInput = regionInput.GetSize();
    typename RegionType::IndexType indexInput = regionInput.GetIndex();
    m_VolumeSize=1;
    for ( int i = 0; i < VImageDimension-1; ++i ) 
      m_VolumeSize *= sizeInput[i];

    sizeInput[m_VectorAxis] = 1;
    indexInput[m_VectorAxis] = 0;
    sizeInput[3] = 1;
    indexInput[3] = 0;
    this->m_Region.SetIndex(indexInput);
    this->m_Region.SetSize(sizeInput);

    if ( this->m_Region.GetNumberOfPixels() > 0 ) // If region is non-empty
      {
      const RegionType & bufferedRegion = this->m_Image->GetBufferedRegion();
      itkAssertOrThrowMacro( ( bufferedRegion.IsInside(this->m_Region) ),
                             "Region " << this->m_Region << " is outside of buffered region " << bufferedRegion );
      }

    // Compute the start offset
    this->m_Offset = this->m_Image->ComputeOffset( this->m_Region.GetIndex() );
    this->m_BeginOffset = this->m_Offset;

    // Compute the end offset. If any component of m_Region.GetSize()
    // is zero, the region is not valid and we set the EndOffset
    // to be same as BeginOffset so that iterator end condition is met
    // immediately.
    if ( this->m_Region.GetNumberOfPixels() == 0 )
      {
      // region is empty, probably has a size of 0 along one dimension
      this->m_EndOffset = this->m_BeginOffset;
      }
    else
      {
      IndexType ind( this->m_Region.GetIndex() );
      SizeType  size( this->m_Region.GetSize() );
      for ( unsigned int i = 0; i < VImageDimension; ++i )
        {
        if (i!=m_VectorAxis)
          ind[i] += ( static_cast< IndexValueType >( size[i] ) - 1 );
        else
          ind[i] = 0;
        }
      this->m_EndOffset = this->m_Image->ComputeOffset(ind);
      this->m_EndOffset++;
      }
  }


  /** Set the pixel value */
  void SetVector(const PixelVectorType & value, const int offIndex=0)
  {
  int off=offIndex*m_VolumeSize;
  for ( int i = 0; i < m_VectorSize; ++i ) 
    {
    this->m_PixelAccessorFunctor.Set(
      *( const_cast< InternalPixelType * >( this->m_Buffer ) + this->m_Offset + off), value[i]);
    off += m_VectorStride;
    }
  }
  
  /** Get the pixel value */
  void GetVector(PixelVectorType& vec, const int offIndex=0) const
  { 
  vec.SetSize(m_VectorSize);
  int off=offIndex*m_VolumeSize;
  for ( int i = 0; i < m_VectorSize; ++i ) 
    {
    vec[i] = this->m_PixelAccessorFunctor.Get( *( this->m_Buffer + this->m_Offset + off ) );
    off += m_VectorStride;
    }
  }
  

protected:

  /** This constructor is declared protected in order to enforce
    const-correctness */
  VectorImageRegionIterator(const ImageRegionIterator< ImageType > & it);
  Self & operator=(const ImageRegionIterator< ImageType > & it);

protected: //made protected so other iterators can access

  OffsetValueType m_VolumeSize;
  OffsetValueType m_VectorStride;
  int m_VectorSize;
  int m_VectorAxis;

};


} // end namespace itk


#endif

