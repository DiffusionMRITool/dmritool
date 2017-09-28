/*=========================================================================

 Program:   Spatially Dense Sparse Vector Image

 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef _itkSpatiallyDenseSparseVectorImage_hxx
#define _itkSpatiallyDenseSparseVectorImage_hxx

#include "itkSpatiallyDenseSparseVectorImage.h"
#include "itkProcessObject.h"

namespace itk
{

template<class TValueType, unsigned int VImageDimension, typename TKeyType>
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::SpatiallyDenseSparseVectorImage()
{
  m_Container = PixelContainer::New();
  m_FillBufferValue.SetSize(0);
}

template<class TValueType, unsigned int VImageDimension, typename TKeyType>
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::Allocate(bool UseDefaultConstructor)
{
  SizeValueType num;
  this->ComputeOffsetTable();
  num = this->GetOffsetTable()[VImageDimension];

  m_Container->Reserve(num, UseDefaultConstructor);
}

template<class TValueType, unsigned int VImageDimension, typename TKeyType>
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::Initialize()
{
  //
  // We don't modify ourselves because the "ReleaseData" methods depend upon
  // no modification when initialized.
  //

  // Call the superclass which should initialize the BufferedRegion ivar.
  Superclass::Initialize();

  // Replace the handle to the container. This is the safest thing to do,
  // since the same container can be shared by multiple images (e.g.
  // Grafted outputs and in place filters).
  m_Container = PixelContainer::New();
}

template<class TValueType, unsigned int VImageDimension, typename TKeyType>
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::FillBuffer(const PixelType& value)
{
  m_FillBufferValue = value;

  const SizeValueType numberOfPixels =
    this->GetBufferedRegion().GetNumberOfPixels();

  for (SizeValueType i=0; i<numberOfPixels; i++)
    {
    (*m_Container)[i].Clear();
    }
}

template<class TValueType, unsigned int VImageDimension, typename TKeyType>
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::SetPixelContainer(PixelContainer *container)
{
  if (m_Container != container)
    {
    m_Container = container;
    this->Modified();
    }
}


template<class TValueType, unsigned int VImageDimension, typename TKeyType>
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::Graft(const DataObject *data)
{
  // call the superclass' implementation
  Superclass::Graft( data );

  if ( data )
    {
    // Attempt to cast data to an Image
    const Self *imgData;

    try
      {
      imgData = dynamic_cast< const Self *>( data );
      }
    catch( ... )
      {
      return;
      }

    // Copy from SpatiallyDenseSparseVectorImage< TPixel, VImageDimension >
    if ( imgData )
      {
      // Now copy anything remaining that is needed
      this->SetPixelContainer( const_cast< PixelContainer *>
                                    (imgData->GetPixelContainer()) );
      }
    else
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::SpatiallyDenseSparseVectorImage::Graft() cannot cast "
                         << typeid(data).name() << " to "
                         << typeid(const Self *).name() );
      }
    }
}

template< class TValueType, unsigned int VImageDimension, typename TKeyType >
unsigned int
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::GetNumberOfComponentsPerPixel() const
{
  return this->m_VectorLength;
}

template< class TValueType, unsigned int VImageDimension, typename TKeyType >
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::SetNumberOfComponentsPerPixel(unsigned int n)
{
  this->SetVectorLength( static_cast< VectorLengthType >( n ) );
  if ( m_FillBufferValue.GetNumberOfElements() == 0 )
    {
    m_FillBufferValue.SetSize(n);
    m_FillBufferValue.Fill(0);
    }
}

template<class TValueType, unsigned int VImageDimension, typename TKeyType>
void
SpatiallyDenseSparseVectorImage<TValueType, VImageDimension, TKeyType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "PixelContainer: " << std::endl;
  m_Container->Print(os, indent.GetNextIndent());
// m_Origin and m_Spacing are printed in the Superclass
}

} // end namespace itk

#endif
