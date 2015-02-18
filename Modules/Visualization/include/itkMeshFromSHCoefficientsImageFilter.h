/*=========================================================================

 Program:   Mesh From SH Coefficients Image Filter

 Copyright (c) Pew-Thian Yap, Jian Cheng. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkMeshFromSHCoefficientsImageFilter_h
#define __itkMeshFromSHCoefficientsImageFilter_h


#include "itkMeshFromContinuousSphericalFunctionImageFilter.h"

namespace itk
{

/** \class MeshFromSHCoefficientsImageFilter
 * \brief Compute mesh from SH coefficients.
 *
 *
 * \author  Pew-Thian Yap, Jian Cheng,  UNC-CH
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromSHCoefficientsImageFilter 
  : public MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromSHCoefficientsImageFilter Self;
  typedef MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromSHCoefficientsImageFilter, MeshFromContinuousSphericalFunctionImageFilter );

  //** Convenient typedefs */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::IndexType InputImageIndexType;
  typedef typename InputImageType::SizeType InputImageSizeType;
  typedef typename InputImageType::SizeValueType InputImageSizeValueType;
  typedef typename InputImageType::SpacingType InputImageSpacingType;
  typedef typename InputImageType::PixelType InputImagePixelType;
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename InputImageType::PointType InputImagePointType;

  typedef TOutputMesh OutputMeshPolyDataType;

  typedef typename Superclass::OutputMeshPointsType   OutputMeshPointsType;
  typedef typename Superclass::OutputMeshCellArrayType  OutputMeshCellArrayType;
  typedef typename Superclass::OutputMeshPolyDataPointer  OutputMeshPolyDataPointer;
  typedef typename Superclass::OutputMeshScalarType  OutputMeshScalarType;
  typedef typename Superclass::OutputMeshRGBType  OutputMeshRGBType;
  
  
  /** Set/Get the min-max normalization */
  typedef SphereTessellator<double> SphereTessellatorType;
  typedef typename SphereTessellatorType::Pointer SphereTessellatorPointer;
  typedef typename SphereTessellatorType::BasicShapeType BasicShapeType;


  /** Orientation Matrices */
  typedef typename Superclass::MatrixType                  MatrixType;
  typedef typename Superclass::MatrixPointer               MatrixPointer;
  typedef typename Superclass::VectorType                  VectorType;
  typedef typename Superclass::VectorPointer               VectorPointer;
  typedef typename Superclass::STDVectorType               STDVectorType;
  typedef typename Superclass::STDVectorPointer            STDVectorPointer;

  
  /** It extracts the first m_MaxOrder SH basis, if m_MaxOrder < inputPixelDimension. If it is -1, use all SH coefficients.  */
  itkSetMacro(MaxOrder, int);
  itkGetMacro(MaxOrder, int);
  
  MatrixPointer ComputeBasisMatrix()
    {
    utlShowPosition(this->GetDebug());

    InputImageConstPointer inputPtr = this->GetInput();
    unsigned int maxNumberOfBasis = inputPtr->GetNumberOfComponentsPerPixel();
    if (m_MaxOrder<=0)
      m_MaxOrder = utl::DimToRankSH(maxNumberOfBasis);

    utlGlobalException(this->m_Orientations->Rows()==0, "need to set m_Orientations");
    this->m_BasisMatrix = utl::ComputeSHMatrix(m_MaxOrder, *this->m_Orientations, CARTESIAN_TO_SPHERICAL);
    // utl::PrintVnlMatrix(*this->m_BasisMatrix, "m_BasisMatrix 0");
    return this->m_BasisMatrix;
    }
  
  // bool IsPixelIndexVisible(const ImageRegionConstIteratorWithIndex<InputImageType> inputIt );
  
protected:
  MeshFromSHCoefficientsImageFilter() : Superclass()
    {
    m_MaxOrder = -1;
    }
  ~MeshFromSHCoefficientsImageFilter()
    {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    PrintVar1(true, m_MaxOrder, os<<indent);
    }
  
  // virtual void GenerateData();
  // void ThreadedGenerateData(const typename TInputImage::RegionType& regionForThread,ThreadIdType threadId );
  // void BeforeThreadedGenerateData();

  int m_MaxOrder;
  

private:
  MeshFromSHCoefficientsImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented


};

} // end namespace itk


#endif
