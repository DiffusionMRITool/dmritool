/**
 *       @file  itkMeshFromSphericalFunctionTessellatedSamplesImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromSphericalFunctionTessellatedSamplesImageFilter_h
#define __itkMeshFromSphericalFunctionTessellatedSamplesImageFilter_h

#include "itkMeshFromContinuousSphericalFunctionImageFilter.h"

namespace itk
{

/** \class MeshFromSphericalFunctionTessellatedSamplesImageFilter
 * \brief Compute mesh from SH coefficients.
 *
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromSphericalFunctionTessellatedSamplesImageFilter 
  : public MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromSphericalFunctionTessellatedSamplesImageFilter Self;
  typedef MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromSphericalFunctionTessellatedSamplesImageFilter, MeshFromContinuousSphericalFunctionImageFilter );

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

  itkSetMacro(DataOrientations, MatrixPointer);
  itkGetMacro(DataOrientations, MatrixPointer);
  
  MatrixPointer ComputeBasisMatrix()
    {
    utlShowPosition(this->GetDebug());

    InputImageConstPointer inputPtr = this->GetInput();
    int numberOfSamples = this->m_Orientations->Rows();
    int inputDimension =  inputPtr->GetNumberOfComponentsPerPixel();

    utlSAGlobalException(m_DataOrientations->Rows()!=inputDimension)
      (m_DataOrientations->Rows())(inputDimension).msg("inconsistent size between input orientations and input image with spherical samples (whole sphere or hemisphere).");
    utlSAGlobalException(numberOfSamples!=inputDimension && numberOfSamples!=2*inputDimension)
      (numberOfSamples)(inputDimension).msg("inconsistent size between tessellation and input image with spherical samples (whole sphere or hemisphere).");

    *this->m_BasisMatrix = (*this->m_Orientations) * m_DataOrientations->GetTranspose();
    // utl::PrintUtlMatrix(*this->m_Orientations, "m_Orientations");
    // utl::PrintUtlMatrix(*this->m_DataOrientations, "m_DataOrientations");
    // utl::PrintUtlMatrix(*this->m_BasisMatrix, "m_BasisMatrix");
    for ( int j = 0; j < this->m_BasisMatrix->Columns(); j += 1 ) 
      {
      int num = 0;
      for ( int i = 0; i < this->m_BasisMatrix->Rows(); i += 1 ) 
        {
        if ( std::fabs((*this->m_BasisMatrix)(i,j)-1.0) < 1e-8 || std::fabs((*this->m_BasisMatrix)(i,j)+1.0)<1e-8  )
          {
          (*this->m_BasisMatrix)(i,j) = 1.0;
          num++;
          }
        else
          (*this->m_BasisMatrix)(i,j) = 0;
        }
      utlGlobalException(inputDimension==numberOfSamples && num!=1, "input orientation is different from the stored one.");
      utlGlobalException(inputDimension==2*numberOfSamples && num!=2, "input orientation is different from the stored one.");
      }

    return this->m_BasisMatrix;
    }
  
  
protected:
  MeshFromSphericalFunctionTessellatedSamplesImageFilter() : Superclass(), 
    m_DataOrientations (new MatrixType())
    {
    }
  ~MeshFromSphericalFunctionTessellatedSamplesImageFilter()
    {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    utl::PrintUtlMatrix(*this->m_DataOrientations, "m_DataOrientations", " ", os<<indent);
    }
  
  MatrixPointer m_DataOrientations;

private:
  MeshFromSphericalFunctionTessellatedSamplesImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented


};

} // end namespace itk


#endif
