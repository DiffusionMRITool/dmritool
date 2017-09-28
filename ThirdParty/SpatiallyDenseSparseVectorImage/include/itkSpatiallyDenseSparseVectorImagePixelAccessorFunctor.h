/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image Pixel Accessor Functor

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImagePixelAccessorFunctor_h
#define __itkSpatiallyDenseSparseVectorImagePixelAccessorFunctor_h

#include "itkMacro.h"


namespace itk
{
/** \class SpatiallyDenseSparseVectorImagePixelAccessorFunctor
 * \brief Provides accessor interfaces to Access pixels and is meant to be
 * used by iterators.
 *
 * A typical user should not need to use this class. The class is internally
 * used by the neighborhood iterators.
 *
 * The pixel accessor is set with the SetPixelAccessor method. This accessor is
 * meant to be used only for SpatiallyDenseSparseVectorImage and not for Image.
 *
 * \sa DefaultSliceContiguousPixelAccessor
 * \sa DefaultPixelAccessor
 * \sa DefaultPixelAccessorFunctor
 *
 * \ingroup ITKSpatiallyDenseSparseVectorImage
 *
 */
template <class TImageType >
class ITK_EXPORT SpatiallyDenseSparseVectorImagePixelAccessorFunctor
{
public:

  typedef TImageType                                   ImageType;
  typedef typename ImageType::InternalPixelType        InternalPixelType;
  typedef typename ImageType::PixelType                ExternalPixelType;
  typedef typename ImageType::AccessorType             PixelAccessorType;
  typedef unsigned int                                 VectorLengthType;

  /** Set the PixelAccessor. This is set at construction time by the image iterators.
   * The type PixelAccessorType is obtained from the ImageType over which the iterators
   * are templated.
   * */
  inline void SetPixelAccessor( PixelAccessorType& accessor )
  {
    m_PixelAccessor = accessor;
  }

  /** Set the pointer index to the start of the buffer. */
  inline void SetBegin( const InternalPixelType * begin ) // NOTE: begin is always 0
  {
    this->m_Begin = const_cast< InternalPixelType * >( begin );
  }

  /** Set output using the value in input */
  inline void Set( InternalPixelType & output, const ExternalPixelType &input ) const
  {
    m_PixelAccessor.Set( output, input, (&output)-m_Begin ); // NOTE: begin is always 0
  }

  /** Get the value from input */
  inline ExternalPixelType Get( const InternalPixelType &input ) const
  {
    return m_PixelAccessor.Get( input, (&input)-m_Begin ); // NOTE: begin is always 0
  }

  /** Required for some filters to compile. */
  static void SetVectorLength( ImageType *image, VectorLengthType length)
  {
    image->SetVectorLength(length);
  }

  /** Required for some filters to compile. */
  static VectorLengthType GetVectorLength( const ImageType *image)
  {
    return image->GetVectorLength();
  }

private:
  PixelAccessorType m_PixelAccessor; // The pixel accessor
  InternalPixelType *m_Begin; // Begin of the buffer, always 0
};

}

#endif
