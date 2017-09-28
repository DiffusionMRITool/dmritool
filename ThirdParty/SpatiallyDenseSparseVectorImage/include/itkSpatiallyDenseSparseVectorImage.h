/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImage_h
#define __itkSpatiallyDenseSparseVectorImage_h

#include "itkImage.h"
#include "itkSpatiallyDenseSparseVectorImagePixelAccessor.h"
#include "itkSpatiallyDenseSparseVectorImagePixelAccessorFunctor.h"
#include "itkSpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor.h"
#include "itkNeighborhoodAccessorFunctor.h"
#include "itkVariableLengthVector.h"
#include "itkSparseVector.h"


namespace itk
{

/** \class SpatiallyDenseSparseVectorImage
 * \brief An n-dimensional vector image with a sparse memory model.
 *
 * The elements for this image are stored in a hash table, catering to very
 * large images with a small number of relevant pixels. The image is spatially
 * dense, allowing fast random access and thread-safe operations.
 *
 * \ingroup ITKSpatiallyDenseSparseVectorImage
 *
 */
template <class TValueType, unsigned int VImageDimension, typename TKeyType = unsigned long>
class ITK_EXPORT SpatiallyDenseSparseVectorImage :
    public ImageBase< VImageDimension >
{
public:
  /** Standard class typedefs */
  typedef SpatiallyDenseSparseVectorImage  Self;
  typedef ImageBase< VImageDimension >     Superclass;
  typedef SmartPointer<Self>               Pointer;
  typedef SmartPointer<const Self>         ConstPointer;
  typedef WeakPointer<const Self>          ConstWeakPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SpatiallyDenseSparseVectorImage, Image);

  /** Pixel typedef support. */
  typedef VariableLengthVector< TValueType > PixelType;

  /** This is the actual pixel type contained in the buffer. */
  typedef TValueType ValueType;
  typedef TKeyType KeyType;
  typedef SparseVector<TValueType, TKeyType> InternalPixelType;

  /** Typedef alias for PixelType */
  typedef InternalPixelType IOPixelType;

  /** Dimension of the image.  This constant is used by functions that are
   * templated over image type (as opposed to being templated over pixel type
   * and dimension) when they need compile time access to the dimension of
   * the image. */
  itkStaticConstMacro( ImageDimension, unsigned int, VImageDimension );

  /** Container used to store pixels in the image. */
  typedef ImportImageContainer< SizeValueType, InternalPixelType > PixelContainer;

  /** A pointer to the pixel container. */
  typedef typename PixelContainer::Pointer PixelContainerPointer;
  typedef typename PixelContainer::ConstPointer PixelContainerConstPointer;

  /** Map used to store data. */
  typedef typename InternalPixelType::InternalDataType PixelMapType;
  typedef typename PixelMapType::iterator PixelMapIterator;
  typedef typename PixelMapType::const_iterator PixelMapConstIterator;

  /** Index typedef support. An index is used to access pixel values. */
  typedef typename Superclass::IndexType IndexType;

  /** Offset typedef support. An offset is used to access pixel values. */
  typedef typename Superclass::OffsetType OffsetType;

  /** Size typedef support. A size is used to define region bounds. */
  typedef typename Superclass::SizeType SizeType;

  /** Direction typedef support. A matrix of direction cosines. */
  typedef typename Superclass::DirectionType DirectionType;

  /** Region typedef support. A region is used to specify a subset of an image. */
  typedef typename Superclass::RegionType RegionType;

  /** Spacing typedef support.  Spacing holds the size of a pixel.  The
   * spacing is the geometric distance between image samples. */
  typedef typename Superclass::SpacingType SpacingType;

  /** Origin typedef support.  The origin is the geometric coordinates
   * of the index (0,0). */
  typedef typename Superclass::PointType PointType;

  /** Offset typedef (relative position between indices) */
  typedef typename Superclass::OffsetValueType OffsetValueType;

  typedef unsigned long VectorLengthType;

  /** Accessor type that convert data between internal and external
   *  representations. */
  typedef SpatiallyDenseSparseVectorImagePixelAccessor<
    ValueType, KeyType > AccessorType;

  /** Tyepdef for the functor used to access pixels.*/
  typedef SpatiallyDenseSparseVectorImagePixelAccessorFunctor< Self >
    AccessorFunctorType;

  /** Tyepdef for the functor used to access a neighborhood of pixel pointers.*/
  typedef SpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor< Self >
    NeighborhoodAccessorFunctorType;

  /** Allocate the image memory. The size of the image must
   * already be set, e.g. by calling SetRegions(). */
//  void Allocate();
  void Allocate(bool UseDefaultConstructor = false) ITK_OVERRIDE;

  /** Convenience methods to set the LargestPossibleRegion,
   *  BufferedRegion and RequestedRegion. Allocate must still be called.
   */
//  using Superclass::SetRegions;
  void SetRegions(const RegionType & region)
    {
    this->SetLargestPossibleRegion(region);
    this->SetBufferedRegion(region);
    this->SetRequestedRegion(region);
    }

  void SetRegions(const SizeType & size)
    {
    RegionType region; region.SetSize(size);
    this->SetLargestPossibleRegion(region);
    this->SetBufferedRegion(region);
    this->SetRequestedRegion(region);
    }

  /** Buffered region has no meaning for sparse images.
   *  This method does nothing except for computing offset table
   * \sa ImageRegion, SetLargestPossibleRegion(), SetRequestedRegion() */
  virtual void SetBufferedRegion(const RegionType &region)
    {
    this->ComputeOffsetTable();
    this->Modified();
    }

  /** Buffered region has no meaning for sparse images.
   *  This method always returns the largest possible region.
   * \sa ImageRegion, SetLargestPossibleRegion(), SetRequestedRegion() */
  virtual const RegionType& GetBufferedRegion() const
  { return this->GetLargestPossibleRegion(); }

  /** Restore the data object to its initial state. This means releasing
   * memory. */
  virtual void Initialize();

  OffsetValueType ComputeOffset(const IndexType &ind) const
  {
    OffsetValueType offset=0;
    const OffsetValueType *offsetTable = this->GetOffsetTable();

    for (int i=ImageDimension-1; i > 0; i--)
      {
      offset += ind[i]*offsetTable[i];
      }
    offset += ind[0];

    return offset;
  }

  /** Fill the image buffer with a value.  Be sure to call Allocate()
   * first. */
  void FillBuffer(const PixelType& value);

  /** \brief Set a pixel value.
   *
   * Allocate() needs to have been called first -- for efficiency,
   * this function does not check that the image has actually been
   * allocated yet. */
   void SetPixel( const IndexType &index, const PixelType &value )
    {
    PixelMapType *map = (this->GetInternalPixel(index)).GetDataPointer();

    map->clear();

    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
      if ( value[i] != 0 )
        {
        (*map)[i] = value[i];
        }
      }
    }

  /** \brief Get a pixel (read only version).
   *
   * For efficiency, this function does not check that the
   * image has actually been allocated yet. Note that the method returns a
   * pixel on the stack. */
  const PixelType GetPixel(const IndexType &index) const
    {
    PixelType pixel;
    pixel.SetSize(m_VectorLength);

    const PixelMapType *map = (this->GetInternalPixel(index)).GetDataPointer();

    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
      typename PixelMapType::const_iterator it = map->find( i );

      if ( it == map->end() )
        {
        pixel[i] = m_FillBufferValue[i];
        }
      else
        {
        pixel[i] = it->second;
        }
      }

    return pixel;
    }

  /** \brief Get a reference to a pixel (e.g. for editing).
   *
   * For efficiency, this function does not check that the
   * image has actually been allocated yet. */
  PixelType GetPixel(const IndexType &index)
  {
    PixelType pixel;
    pixel.SetSize(m_VectorLength);

    const PixelMapType *map = (this->GetInternalPixel(index)).GetDataPointer();

    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
      typename PixelMapType::const_iterator it = map->find( i );

      if ( it == map->end() )
        {
        pixel[i] = m_FillBufferValue[i];
        }
      else
        {
        pixel[i] = it->second;
        }
      }

    return pixel;
  }

  /** \brief Get an internal pixel (read only version).
   *
   * For efficiency, this function does not check that the
   * image has actually been allocated yet. Note that the method returns a
   * pixel on the stack. */
  const InternalPixelType & GetInternalPixel(const IndexType &index) const
    {
    OffsetValueType offset = this->ComputeOffset(index);
    return (*m_Container)[offset];
    }

  /** \brief Get an internal pixel (e.g. for editing).
   *
   * For efficiency, this function does not check that the
   * image has actually been allocated yet. Note that the method returns a
   * pixel on the stack. */
  InternalPixelType & GetInternalPixel(const IndexType &index)
    {
    OffsetValueType offset = this->ComputeOffset(index);
    return (*m_Container)[offset];
    }

  /** \brief Access a pixel. This version cannot be an lvalue because the pixel
   * is converted on the fly to a VariableLengthVector.
   *
   * For efficiency, this function does not check that the
   * image has actually been allocated yet. */
  PixelType operator[](const IndexType &index)
     { return this->GetPixel(index); }

  /** \brief Access a pixel.
   *
   * For efficiency, this function does not check that the
   * image has actually been allocated yet. */
  PixelType operator[](const IndexType &index) const
     { return this->GetPixel(index); }

  /** Return a pointer to the container. */
  PixelContainer* GetPixelContainer()
  {
    return m_Container.GetPointer();
  }

  /** Return a pointer to the container. */
  const PixelContainer* GetPixelContainer() const
  {
    return m_Container.GetPointer();
  }

  /** Sparse images do not have buffer. This method always returns 0 */
  InternalPixelType * GetBufferPointer()
  {
    return m_Container->GetBufferPointer();
  }
  const InternalPixelType * GetBufferPointer() const
  {
    return m_Container->GetBufferPointer();
  }

  /** Set the container to use. Note that this does not cause the
   * DataObject to be modified. */
  void SetPixelContainer( PixelContainer *container );

  /** Graft the data and information from one image to another. This
   * is a convenient method to setup a second image with all the meta
   * information of another image and use the same pixel
   * container. Note that this method is different than just using two
   * SmartPointers to the same image since separate DataObjects are
   * still maintained. This method is similar to
   * ImageSource::GraftOutput(). The implementation in ImageBase
   * simply calls CopyInformation() and copies the region ivars.
   * The implementation here refers to the superclass' implementation
   * and then copies over the pixel container. */
  virtual void Graft(const DataObject *data);

  /** Return the Pixel Accessor object */
  AccessorType GetPixelAccessor( void )
  {
    return AccessorType(m_FillBufferValue, m_VectorLength);
  }

  /** Return the Pixel Accesor object */
  const AccessorType GetPixelAccessor( void ) const
  {
    return AccessorType(m_FillBufferValue, m_VectorLength);
  }

  /** Return the NeighborhoodAccessor functor */
  NeighborhoodAccessorFunctorType GetNeighborhoodAccessor()
  {
    return NeighborhoodAccessorFunctorType(m_FillBufferValue, m_VectorLength);
  }

  /** Return the NeighborhoodAccessor functor */
  const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
  {
    return NeighborhoodAccessorFunctorType(m_FillBufferValue, m_VectorLength);
  }

  /** Set/Get macros for the length of each vector in the vector image */
  itkSetMacro(VectorLength, VectorLengthType);
  itkGetConstReferenceMacro(VectorLength, VectorLengthType);

  /** Get/Set the number of components each pixel has, ie the VectorLength */
  virtual unsigned int GetNumberOfComponentsPerPixel() const;

  virtual void SetNumberOfComponentsPerPixel(unsigned int n);

protected:
  SpatiallyDenseSparseVectorImage();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  virtual ~SpatiallyDenseSparseVectorImage() {};

private:
  SpatiallyDenseSparseVectorImage( const Self & ); // purposely not implementated
  void operator=(const Self&); //purposely not implemented

  /** Memory for the map containing the pixel data. */
  PixelContainerPointer m_Container;
  PixelType m_FillBufferValue;

  /** Length of the "vector pixel" */
  VectorLengthType m_VectorLength;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatiallyDenseSparseVectorImage.hxx"
#endif

#endif
