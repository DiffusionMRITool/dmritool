/**
 *       @file  itkVectorImageRegionIteratorWithIndex.h
 *      @brief  
 *     Created  "06-19-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef itkVectorImageRegionIteratorWithIndex_h
#define itkVectorImageRegionIteratorWithIndex_h

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "utlITKConceptChecking.h"
#include "itkVectorImage.h"
#include "itkVariableLengthVector.h"
#include "utlCoreMacro.h"
#include "utlITK.h"

namespace itk
{

/**
 * \class VectorImageRegionIteratorWithIndex
 * \brief A multi-dimensional iterator templated over image type. It provides the same interfaces for both itk::Image<T, N+1> and itk::VectorImage<T, N>
 *
 * \author Jian Cheng
 * 
 *  If  image4d is a \c itk::Image<double,4> or \c itk::VectorImage<double,3> object, and k is 0,1,2,3 axis, 
 *  then the iterator could loop over axises except k-axis, and obtain a vector a long k-axis
 *
 * \code
 *
 * itk::VectorImageRegionIteratorWithIndex<NDImageType> it(image4d, image4d->GetLargestPossibleRegion(),k);
 * int sizeVec = itk::GetVectorImageVectorSize(image4d); // number of components in t-axis
 * itk::VariableLengthVector<double> vec; // vector along k-axis
 * for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
 *   {
 *   // i means the coordinate in the last axis (t-axis)
 *   for ( int i = 0; i < (k==3?1:sizeVec); ++i ) 
 *     {
 *     index = it.GetIndex();
 *     it.GetVector(vec, i);
 *
 *     // do somthing for vec ...
 *
 *     it.SetVector(vec, i);
 *     }
 *   }
 *
 * \endcode
 *
 * \ingroup ImageIterators
 * \ingroup ITKCommon
 *
 * \sa itkVectorImageRegionIteratorWithIndexTest.cxx 
 */

template< typename TImage >
class VectorImageRegionIteratorWithIndex : public ImageRegionIteratorWithIndex< TImage >
{
public:
  /** Standard class typedefs. */
  typedef VectorImageRegionIteratorWithIndex Self;

  /** Dimension of the image the iterator walks.  This constant is needed so
   * functions that are templated over image iterator type (as opposed to
   * being templated over pixel type and dimension) can have compile time
   * access to the dimension of the image that the iterator walks. */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int,
                      TImage::ImageDimension);

  /** Define the superclass */
  typedef ImageRegionIteratorWithIndex< TImage > Superclass;
  

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

  typedef ImageRegion<TImage::ImageDimension+1>   NDImageRegionType;
  typedef ImageRegion<TImage::ImageDimension>     VectorImageRegionType;

  typedef VariableLengthVector< InternalPixelType >  PixelVectorType;
  
#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( SameTypeCheck,
                   ( itk::Concept::SameType< ImageType, VectorImage<InternalPixelType, ImageIteratorDimension> > ) );
  itkConceptMacro( SameTypeCheck2,
                   ( itk::Concept::SameType< PixelType, PixelVectorType> ) );
#endif


  /** Default Constructor. Need to provide a default constructor since we
   * provide a copy constructor. */
  VectorImageRegionIteratorWithIndex() : Superclass() 
  {
  m_VectorStride=0;
  m_VectorSize=0;
  m_VectorAxis=ImageIteratorDimension;
  m_BeginBuffer   = ITK_NULLPTR;
  }

  /** Default Destructor */
  ~VectorImageRegionIteratorWithIndex() {}

  /** Copy Constructor. The copy constructor is provided to make sure the
   * handle to the image is properly reference counted. */
  VectorImageRegionIteratorWithIndex(const Self & it) : Superclass(it) 
  { 
  m_VectorStride = it.m_VectorStride;
  m_VectorSize = it.m_VectorSize;
  m_VectorAxis = it.m_VectorAxis;
  m_BeginBuffer = it.m_BeginBuffer;
  };

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  VectorImageRegionIteratorWithIndex(TImage *ptr, const RegionType & regionInput, int vectorAxis=-1) //: Superclass(ptr, regionInput) 
  { 
  Initialize(ptr, regionInput, vectorAxis);
  } 
  
  // VectorImageRegionIteratorWithIndex(TImage *ptr, const NDImageRegionType & regionInput, int vectorAxis=-1) : Superclass(ptr, regionInput) 
  // { 
  // RegionType region;
  // itk::CopyImageRegion(regionInput, region);
  // Initialize(ptr, region, vectorAxis);
  // } 

  /** operator= is provided to make sure the handle to the image is properly
   * reference counted. */
  Self & operator=(const Self & it) 
    {
    Superclass::operator=(it);
    m_VectorStride = it.m_VectorStride;
    m_VectorSize = it.m_VectorSize;
    m_VectorAxis = it.m_VectorAxis;
    m_BeginBuffer = it.m_BeginBuffer;
    return *this;
    };

  /** Set the pixel value. offindex is the coordinate in t-axis, it is only used when m_VectorAxis!=ImageIteratorDimension. */
  void SetVector(const PixelVectorType & value, const int offIndex=0) const
  {
  if (m_VectorAxis==ImageIteratorDimension)
    {
    this->m_PixelAccessorFunctor.Set(*( const_cast< InternalPixelType * >( this->m_Position ) ), value);
    }
  else
    {
    int numberOfComponens = this->m_Image->GetNumberOfComponentsPerPixel();
    InternalPixelType* p = const_cast<InternalPixelType*>(this->m_BeginBuffer + numberOfComponens*(this->m_Position-this->m_Begin));
    int off = offIndex;
    for ( int i = 0; i < m_VectorSize; i += 1 ) 
      {
      *(p+off) = value[i];
      off += m_VectorStride;
      }
    }
  }
  
  /** Get the pixel value. offindex is the coordinate in t-axis, it is only used when m_VectorAxis!=ImageIteratorDimension. */
  void GetVector(PixelVectorType& vec, const int offIndex=0) const
  { 
  if (m_VectorAxis==ImageIteratorDimension)
    {
    vec = this->m_PixelAccessorFunctor.Get( *this->m_Position );
    }
  else
    {
    vec.SetSize(m_VectorSize);
    int numberOfComponens = this->m_Image->GetNumberOfComponentsPerPixel();
    InternalPixelType* p = const_cast<InternalPixelType*>(this->m_BeginBuffer + numberOfComponens*(this->m_Position-this->m_Begin));
    int off = offIndex;
    for ( int i = 0; i < m_VectorSize; i += 1 ) 
      {
      vec[i] = *(p+off);
      off += m_VectorStride;
      }
    }
  }
  


protected:

  void Initialize(TImage *ptr, const RegionType & regionInput, int vectorAxis=-1)
    {
    // -1 means t axis
    m_VectorAxis = (vectorAxis==-1) ? ImageIteratorDimension : vectorAxis;
    itkAssertOrThrowMacro((m_VectorAxis<=(int)ImageIteratorDimension && m_VectorAxis>=0), "wrong vector axis, m_VectorAxis=" << m_VectorAxis << ", ImageIteratorDimension=" << ImageIteratorDimension);
      
    this->m_Image = ptr;
    const InternalPixelType *buffer   = this->m_Image->GetBufferPointer();

    if (m_VectorAxis==ImageIteratorDimension)
      {
      // Superclass(ptr, regionInput);

      this->m_BeginIndex        = regionInput.GetIndex();
      this->m_PositionIndex     = this->m_BeginIndex;
      this->m_Region            = regionInput;

      if ( regionInput.GetNumberOfPixels() > 0 ) // If region is non-empty
        {
        const RegionType & bufferedRegion = this->m_Image->GetBufferedRegion();
        itkAssertOrThrowMacro( ( bufferedRegion.IsInside(this->m_Region) ),
                               "Region " << this->m_Region << " is outside of buffered region " << bufferedRegion );
        }

      std::copy(this->m_Image->GetOffsetTable(),
                this->m_Image->GetOffsetTable()+Superclass::ImageDimension + 1 ,
                this->m_OffsetTable);

      // Compute the start position
      OffsetValueType offs =  this->m_Image->ComputeOffset(this->m_BeginIndex);
      this->m_Begin = buffer + offs;
      this->m_Position = this->m_Begin;

      // Compute the end offset
      this->m_Remaining = false;
      IndexType pastEnd;
      for ( unsigned int i = 0; i < Superclass::ImageDimension; ++i )
        {
        SizeValueType size = regionInput.GetSize()[i];
        if ( size > 0 )
          {
          this->m_Remaining = true;
          }
        this->m_EndIndex[i] = this->m_BeginIndex[i] + static_cast< OffsetValueType >( size );
        pastEnd[i]    = this->m_BeginIndex[i] + static_cast< OffsetValueType >( size ) - 1;
        }
      this->m_End = buffer + this->m_Image->ComputeOffset(pastEnd);

      this->m_PixelAccessor = this->m_Image->GetPixelAccessor();
      this->m_PixelAccessorFunctor.SetPixelAccessor(this->m_PixelAccessor);
      this->m_PixelAccessorFunctor.SetBegin(buffer);

      this->GoToBegin();


      m_VectorSize = ptr->GetNumberOfComponentsPerPixel();
      m_VectorStride = 1;
      m_BeginBuffer = buffer + offs;
      }
    else
      {

      typename RegionType::SizeType sizeInput = regionInput.GetSize();
      typename RegionType::IndexType indexInput = regionInput.GetIndex();
      sizeInput[m_VectorAxis] = 1;
      indexInput[m_VectorAxis] = 0;
      RegionType region;
      region.SetIndex(indexInput);
      region.SetSize(sizeInput);


      this->m_BeginIndex        = region.GetIndex();
      this->m_PositionIndex     = this->m_BeginIndex;
      this->m_Region            = region;

      if ( region.GetNumberOfPixels() > 0 ) // If region is non-empty
        {
        const RegionType & bufferedRegion = this->m_Image->GetBufferedRegion();
        itkAssertOrThrowMacro( ( bufferedRegion.IsInside(this->m_Region) ),
  "Region " << this->m_Region << " is outside of buffered region " << bufferedRegion );
        }

      std::copy(this->m_Image->GetOffsetTable(),
          this->m_Image->GetOffsetTable()+Superclass::ImageDimension + 1 ,
          this->m_OffsetTable);

      // Compute the start position
      OffsetValueType offs =  this->m_Image->ComputeOffset(this->m_BeginIndex);
      this->m_Begin = buffer + offs;
      this->m_Position = this->m_Begin;

      // Compute the end offset
      this->m_Remaining = false;
      IndexType pastEnd;
      for ( unsigned int i = 0; i < Superclass::ImageDimension; ++i )
        {
        if (i!=m_VectorAxis)
          {
          SizeValueType size = region.GetSize()[i];
          if ( size > 0 )
            {
            this->m_Remaining = true;
            }
          this->m_EndIndex[i] = this->m_BeginIndex[i] + static_cast< OffsetValueType >( size );
          pastEnd[i]    = this->m_BeginIndex[i] + static_cast< OffsetValueType >( size ) - 1;
          }
        else
          {
          this->m_EndIndex[m_VectorAxis] = 0;
          pastEnd[i] = 0;
          this->m_Remaining = true;
          }
        }
      this->m_End = buffer + this->m_Image->ComputeOffset(pastEnd);

      this->m_PixelAccessor = this->m_Image->GetPixelAccessor();
      this->m_PixelAccessorFunctor.SetPixelAccessor(this->m_PixelAccessor);
      this->m_PixelAccessorFunctor.SetBegin(buffer);

      this->GoToBegin();

      if (ptr)
        {
        SizeType size = ptr->GetLargestPossibleRegion().GetSize();
        const typename ImageType::OffsetValueType* offsetTable = ptr->GetOffsetTable();
        m_VectorStride = ptr->GetNumberOfComponentsPerPixel()*offsetTable[m_VectorAxis];
        m_VectorSize = size[m_VectorAxis];
        m_BeginBuffer = buffer + offs*m_VectorStride;
        }
      }
    }

  /** This constructor is declared protected in order to enforce
    const-correctness */
  VectorImageRegionIteratorWithIndex(const ImageRegionIteratorWithIndex< TImage > & it);
  Self & operator=(const ImageRegionIteratorWithIndex< TImage > & it);

protected: 

  const InternalPixelType *m_BeginBuffer;

  /** stride in vector array  */
  OffsetValueType m_VectorStride;
  /** size alone m_VectorAxis axis  */
  int m_VectorSize;
  /** the axis which vector arrays need to be extract  */
  int m_VectorAxis;
};


/** For Image<TPixel, 4>  */
template <typename TPixel, unsigned int VImageDimension>
class VectorImageRegionIteratorWithIndex<Image<TPixel, VImageDimension> > : public ImageRegionIteratorWithIndex<Image<TPixel, VImageDimension> >
{
public:
  /** Standard class typedefs. */
  typedef VectorImageRegionIteratorWithIndex Self;


  /** Dimension of the image the iterator walks.  This constant is needed so
   * functions that are templated over image iterator type (as opposed to
   * being templated over pixel type and dimension) can have compile time
   * access to the dimension of the image that the iterator walks. */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int,
                      VImageDimension);

  /** Define the superclass */
  typedef ImageRegionIteratorWithIndex< Image<TPixel, VImageDimension> > Superclass;

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
  
  typedef ImageRegion<VImageDimension>   NDImageRegionType;
  typedef ImageRegion<VImageDimension-1>     VectorImageRegionType;
  
  typedef VariableLengthVector< InternalPixelType >  PixelVectorType;

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro( SameTypeCheck,
                   ( itk::Concept::SameType< ImageType, Image<InternalPixelType, VImageDimension> > ) );
#endif

  /** Default Constructor. Need to provide a default constructor since we
   * provide a copy constructor. */
  VectorImageRegionIteratorWithIndex() : Superclass() 
  {
  m_VolumeSize=0;
  m_VectorStride=0;
  m_VectorSize=0;
  m_VectorAxis=ImageIteratorDimension-1;
  }

  /** Default Destructor */
  ~VectorImageRegionIteratorWithIndex() {}

  /** Copy Constructor. The copy constructor is provided to make sure the
   * handle to the image is properly reference counted. */
  VectorImageRegionIteratorWithIndex(const Self & it) : Superclass(it) 
  { 
  m_VolumeSize = it.m_VolumeSize;
  m_VectorStride = it.m_VectorStride;
  m_VectorSize = it.m_VectorSize;
  m_VectorAxis = it.m_VectorAxis;
  };

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  VectorImageRegionIteratorWithIndex(ImageType *ptr, const RegionType & regionInput, int vectorAxis=-1)  
  {  
  Initialize(ptr, regionInput, vectorAxis);
  } 
  
  // VectorImageRegionIteratorWithIndex(ImageType *ptr, const VectorImageRegionType & regionInput, int vectorAxis=-1)
  // { 
  // RegionType region;
  // itk::CopyImageRegion(regionInput, region);
  // Initialize(ptr, region, vectorAxis);
  // } 

  /** operator= is provided to make sure the handle to the image is properly
   * reference counted. */
  Self & operator=(const Self & it) 
    {
    Superclass::operator=(it);
    m_VolumeSize = it.m_VolumeSize;
    m_VectorStride = it.m_VectorStride;
    m_VectorSize = it.m_VectorSize;
    m_VectorAxis = it.m_VectorAxis;
    return *this;
    };


  /** Set the pixel value. offindex is the coordinate in t-axis, it is only used when m_VectorAxis!=ImageIteratorDimension-1. */
  void SetVector(const PixelVectorType & value, const int offIndex=0)
  {
  int off=offIndex*m_VolumeSize;
  for ( int i = 0; i < m_VectorSize; ++i ) 
    {
    this->m_PixelAccessorFunctor.Set(
      *( const_cast< InternalPixelType * >( this->m_Position + off) ), value[i]);
    off += m_VectorStride;
    }
  }
  
  /** Get the pixel value. offindex is the coordinate in t-axis, it is only used when m_VectorAxis!=ImageIteratorDimension-1. */
  void GetVector(PixelVectorType& vec, const int offIndex=0) const
  { 
  vec.SetSize(m_VectorSize);
  int off=offIndex*m_VolumeSize;
  for ( int i = 0; i < m_VectorSize; ++i ) 
    {
    vec[i] = this->m_PixelAccessorFunctor.Get( *( this->m_Position + off ) );
    off += m_VectorStride;
    }
  }
  

protected:
  
  void Initialize(ImageType *ptr, const RegionType & regionInput, int vectorAxis=-1)  
    {

    this->m_Image = ptr;

    const InternalPixelType *buffer   = this->m_Image->GetBufferPointer();

    m_VectorAxis = (vectorAxis==-1) ? Superclass::ImageDimension-1 : vectorAxis;
    itkAssertOrThrowMacro((m_VectorAxis<=(int)ImageIteratorDimension && m_VectorAxis>=0), "wrong vector axis, m_VectorAxis=" << m_VectorAxis);
    
    typename RegionType::SizeType size4d = this->m_Image->GetLargestPossibleRegion().GetSize();
    m_VolumeSize=1;
    for ( int i = 0; i < VImageDimension-1; ++i ) 
      m_VolumeSize *= size4d[i];

    typename RegionType::SizeType sizeInput = regionInput.GetSize();
    typename RegionType::IndexType indexInput = regionInput.GetIndex();

    if (m_VectorAxis<ImageIteratorDimension)
      {
      // set size as 1 in m_VectorAxis and the last t-axis 
      sizeInput[m_VectorAxis] = 1;
      indexInput[m_VectorAxis] = 0;
      sizeInput[ImageIteratorDimension-1] = 1;
      indexInput[ImageIteratorDimension-1] = 0;
      }
    RegionType region;
    region.SetIndex(indexInput);
    region.SetSize(sizeInput);


    this->m_BeginIndex        = region.GetIndex();
    this->m_PositionIndex     = this->m_BeginIndex;
    this->m_Region            = region;

    if ( region.GetNumberOfPixels() > 0 ) // If region is non-empty
      {
      const RegionType & bufferedRegion = this->m_Image->GetBufferedRegion();
      itkAssertOrThrowMacro( ( bufferedRegion.IsInside(this->m_Region) ),
  "Region " << this->m_Region << " is outside of buffered region " << bufferedRegion );
      }

    std::copy(this->m_Image->GetOffsetTable(),
  this->m_Image->GetOffsetTable()+Superclass::ImageDimension + 1 ,
  this->m_OffsetTable);

    // Compute the start position
    OffsetValueType offs =  this->m_Image->ComputeOffset(this->m_BeginIndex);
    this->m_Begin = buffer + offs;
    this->m_Position = this->m_Begin;

    // Compute the end offset
    this->m_Remaining = false;
    IndexType pastEnd;
    for ( unsigned int i = 0; i < Superclass::ImageDimension; ++i )
      {
      if (i!=m_VectorAxis)
        {
        SizeValueType size = region.GetSize()[i];
        if ( size > 0 )
          {
          this->m_Remaining = true;
          }
        this->m_EndIndex[i] = this->m_BeginIndex[i] + static_cast< OffsetValueType >( size );
        pastEnd[i]    = this->m_BeginIndex[i] + static_cast< OffsetValueType >( size ) - 1;
        }
      else
        {
        this->m_EndIndex[m_VectorAxis] = 0;
        pastEnd[i] = 0;
        this->m_Remaining = true;
        }
      }
    this->m_End = buffer + this->m_Image->ComputeOffset(pastEnd);

    this->m_PixelAccessor = this->m_Image->GetPixelAccessor();
    this->m_PixelAccessorFunctor.SetPixelAccessor(this->m_PixelAccessor);
    this->m_PixelAccessorFunctor.SetBegin(buffer);

    this->GoToBegin();

    if (ptr)
      {
      SizeType size = ptr->GetLargestPossibleRegion().GetSize();
      const typename ImageType::OffsetValueType* offsetTable = ptr->GetOffsetTable();
      if (m_VectorAxis<ImageIteratorDimension)
        {
        m_VectorStride = offsetTable[m_VectorAxis];
        m_VectorSize = size[m_VectorAxis];
        }
      else
        {
        m_VectorStride = 0;
        m_VectorSize = 1;
        }
      }

    }

  /** This constructor is declared protected in order to enforce
    const-correctness */
  VectorImageRegionIteratorWithIndex(const ImageRegionIteratorWithIndex< ImageType > & it);
  Self & operator=(const ImageRegionIteratorWithIndex< ImageType > & it);

protected: 

  OffsetValueType m_VolumeSize;
  OffsetValueType m_VectorStride;
  int m_VectorSize;
  int m_VectorAxis;

};

/** For Image<VariableLengthVector<TPixel>, 3>, it works as ImageRegionIteratorWithIndex  */
template <typename TPixel>
class VectorImageRegionIteratorWithIndex<Image<VariableLengthVector<TPixel>, 3 > > : public ImageRegionIteratorWithIndex<Image<VariableLengthVector<TPixel>,3> >
{
public:
  typedef VectorImageRegionIteratorWithIndex                Self;
  typedef ImageRegionIteratorWithIndex< Image<VariableLengthVector<TPixel>, 3> > Superclass;

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


} // end namespace itk


#endif
