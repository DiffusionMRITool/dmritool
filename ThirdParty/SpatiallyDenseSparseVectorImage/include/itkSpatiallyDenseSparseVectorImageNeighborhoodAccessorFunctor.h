/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image Neighborhood Accessor Functor

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor_h
#define __itkSpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor_h

#include "itkImageBoundaryCondition.h"
#include "itkNeighborhood.h"
#include "itkImageBase.h"
#include "itkNumericTraits.h"


namespace itk
{

/** \class SpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor
 * \brief Provides accessor interfaces to Access pixels and is meant to be
 * used on pointers to pixels held by the Neighborhood class.
 *
 * A typical user should not need to use this class. The class is internally
 * used by the neighborhood iterators.
 *
 * \author Pew-Thian Yap (ptyap@med.unc.edu)
 */
template< class TImage >
class SpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor
{
public:
  typedef TImage                                           ImageType;
  typedef typename ImageType::PixelType                    PixelType;
  typedef typename ImageType::InternalPixelType            InternalPixelType;
  typedef typename ImageType::OffsetType                   OffsetType;
  typedef unsigned int                                     VectorLengthType;

  typedef Neighborhood< InternalPixelType *, TImage::ImageDimension>
    NeighborhoodType;

  typedef ImageBoundaryCondition< ImageType > const
                          *ImageBoundaryConditionConstPointerType;

  SpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor(
    PixelType fillBufferValue , VectorLengthType length)
    : m_FillBufferValue( fillBufferValue ), m_VectorLength( length ) {};
  SpatiallyDenseSparseVectorImageNeighborhoodAccessorFunctor()
    : m_FillBufferValue( NumericTraits<PixelType>::Zero, m_VectorLength( 0 ) ) {};

  /** Set the pointer index to the start of the buffer.
   * This must be set by the iterators to the starting location of the buffer.
   * Typically a neighborhood iterator iterating on a neighborhood of an Image,
   * say \c image will set this in its constructor. For instance:
   *
   * \code
   * ConstNeighborhoodIterator( radius, image, )
   *   {
   *   ...
   *   m_NeighborhoodAccessorFunctor.SetBegin( image->GetBufferPointer() );
   *   }
   * \endcode
   */
  inline void SetBegin( const InternalPixelType * begin )  // NOTE: begin is always 0
    { this->m_Begin = const_cast< InternalPixelType * >( begin ); }

  /** Method to dereference a pixel pointer. This is used from the
   * ConstNeighborhoodIterator as the equivalent operation to (*it).
   * This method should be preferred over the former (*it) notation.
   * The reason is that dereferencing a pointer to a location of
   * VectorImage pixel involves a different operation that simply
   * dereferencing the pointer. Here a PixelType (array of InternalPixelType s)
   * is created on the stack and returned.  */
  inline PixelType Get( const InternalPixelType *pixelPointer ) const
    {
    PixelType pixel;
    pixel.SetSize(m_VectorLength);

    const typename InternalPixelType::InternalDataType *map =
        pixelPointer->GetDataPointer();

    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
      typename InternalPixelType::InternalDataType::const_iterator it = map->find( i );

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

  /** Method to set the pixel value at a certain pixel pointer */
  inline void Set( InternalPixelType* &pixelPointer, const PixelType &p ) const
    {
    typename InternalPixelType::InternalDataType *map =
        pixelPointer->GetDataPointer();

    map->clear();

    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
      if ( p[i] != 0 )
        {
        (*map)[i] = p[i];
        }
      }
    }

  inline PixelType BoundaryCondition(
      const OffsetType& point_index,
      const OffsetType &boundary_offset,
      const NeighborhoodType *data,
      const ImageBoundaryConditionConstPointerType boundaryCondition) const
    {
    return boundaryCondition->operator()(point_index, boundary_offset, data, *this);
    }

  /** Required for some filters to compile. */
  void SetVectorLength( VectorLengthType length )
    {
    m_VectorLength = length;
    }

  /** Required for some filters to compile. */
  VectorLengthType GetVectorLength()
    {
    return m_VectorLength;
    }

private:
  PixelType m_FillBufferValue;
  InternalPixelType *m_Begin;  // Begin of the buffer, always 0

  VectorLengthType m_VectorLength;

};

} // end namespace itk
#endif
