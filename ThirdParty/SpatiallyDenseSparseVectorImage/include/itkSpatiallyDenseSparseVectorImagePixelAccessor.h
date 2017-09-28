/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image Pixel Accessor

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkSpatiallyDenseSparseVectorImagePixelAccessor_h
#define __itkSpatiallyDenseSparseVectorImagePixelAccessor_h

#include "itkMacro.h"
#include "itkVariableLengthVector.h"
#include "itkSparseVector.h"
#include "itkIntTypes.h"


namespace itk
{

/** \class SpatiallyDenseSparseVectorImagePixelAccessor
 * \brief Give access to partial aspects of a type
 *
 * SpatiallyDenseSparseVectorImagePixelAccessor is specifically meant to provide
 * SpatiallyDenseSparseVectorImage
 * with the same \c DefaultPixelAccessor interface that
 * DefaultPixelAccessor provides to Image.
 *
 * SpatiallyDenseSparseVectorImagePixelAccessor is templated over an internal
 * type and an external type representation. This class encapsulates a
 * customized convertion between the internal and external
 * type representations.
 *
 * \ingroup ITKSpatiallyDenseSparseVectorImage
 *
 */
template <class TValueType, class TKeyType>
class ITK_EXPORT SpatiallyDenseSparseVectorImagePixelAccessor
{
public:

 /** External typedef. It defines the external aspect
   * that this class will exhibit. */
  typedef VariableLengthVector<TValueType> ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data. */
  typedef SparseVector<TValueType, TKeyType> InternalType;

  typedef unsigned long VectorLengthType;

  /** Set output using the value in input */
  inline void Set(InternalType & output,
                  const ExternalType & input,
                  const SizeValueType offset) const
  {
    output.Clear();
    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
      if ( input[i] != 0 )
        {
        output[i] = input[i];
        }
      }
  }

  /** Get the value from input */
  inline ExternalType Get(const InternalType & input,
                          const SizeValueType offset) const
  {
    ExternalType pixel;
    pixel.SetSize(m_VectorLength);

    for ( VectorLengthType i = 0; i < m_VectorLength; i++ )
      {
//      typename InternalType::const_iterator it = input.find( i );

//      if ( it == input.end() )
//        {
//        pixel[i] = m_FillBufferValue[i];
//        }
//      else
//        {
//        pixel[i] = it->second;
//        }
      pixel[i] = input[i];
      }

    return pixel;
  }

  /** Set the length of each vector in the VectorImage */
  void SetVectorLength(VectorLengthType length)
  {
    m_VectorLength = length;
  }

  /** Get Vector lengths */
  VectorLengthType GetVectorLength() const { return m_VectorLength; }

  SpatiallyDenseSparseVectorImagePixelAccessor() : m_VectorLength(0) {}

   /** Constructor to initialize slices and image size at construction time */
   SpatiallyDenseSparseVectorImagePixelAccessor(ExternalType fillBufferValue,
                                                VectorLengthType length)
   {
     m_FillBufferValue = fillBufferValue;
     m_VectorLength = length;
   }

  virtual ~SpatiallyDenseSparseVectorImagePixelAccessor() {};

private:
  ExternalType m_FillBufferValue;
  VectorLengthType m_VectorLength;
};

} // end namespace itk


#endif
