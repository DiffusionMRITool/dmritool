/**
 *       @file  itkMeshFromDiscreteFiberODFImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromDiscreteFiberODFImageFilter_h
#define __itkMeshFromDiscreteFiberODFImageFilter_h


#include "itkMeshFromContinuousSphericalFunctionImageFilter.h"
#include "itkBasisMatrixGenerator.h"

namespace itk
{

/** \class MeshFromDiscreteFiberODFImageFilter
 * \brief Compute mesh from discrete fiber ODF represented by a set of basis functions.
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromDiscreteFiberODFImageFilter 
  : public MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromDiscreteFiberODFImageFilter Self;
  typedef MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromDiscreteFiberODFImageFilter, MeshFromContinuousSphericalFunctionImageFilter );

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

  typedef typename Superclass::MaskImageType      MaskImageType;
  
  typedef BasisMatrixGenerator<double>  BasisMatrixGeneratorType;
  typedef typename BasisMatrixGeneratorType::Pointer  BasisMatrixGeneratorPointer;
  
  itkSetObjectMacro(BasisMatrixGenerator, BasisMatrixGeneratorType);
  itkGetObjectMacro(BasisMatrixGenerator, BasisMatrixGeneratorType);
  
  itkSetMacro(Radius, double);
  itkGetMacro(Radius, double);
  
  MatrixPointer ComputeBasisMatrix()
    {
    utlShowPosition(this->GetDebug());
    InputImageConstPointer inputPtr = this->GetInput();
    unsigned int inputPixelDimension = inputPtr->GetNumberOfComponentsPerPixel();

    utlGlobalException(!m_BasisMatrixGenerator, "no m_BasisMatrixGenerator");
    utlSAGlobalException(m_BasisMatrixGenerator->GetNumberOfBasis()!=inputPixelDimension)
      (m_BasisMatrixGenerator->GetNumberOfBasis())(inputPixelDimension)
      .msg("no m_BasisMatrixGenerator");

    // this->Print(std::cout<<"this");
    utlGlobalException(this->m_Orientations->Rows()==0, "no m_Orientations");
    if (m_BasisMatrixGenerator->GetOutputType()==BasisMatrixGeneratorType::DWI)
      {
      m_BasisMatrixGenerator->GetSamplingSchemeQSpace()->SetOrientationsCartesian(this->m_Orientations);
      STDVectorPointer bVec = m_BasisMatrixGenerator->GetSamplingSchemeQSpace()->GetBVector();
      utlGlobalException(bVec->size()==0 && m_Radius<0, "need to set m_Radius (b value) for DWI");
      if (bVec->size()==0 && m_Radius>=0)
        utl::MatchBVectorAndGradientMatrix(m_Radius, *bVec, *this->m_Orientations);
      }
    else if (m_BasisMatrixGenerator->GetOutputType()==BasisMatrixGeneratorType::EAP)
      {
      m_BasisMatrixGenerator->GetSamplingSchemeRSpace()->SetOrientationsCartesian(this->m_Orientations);
      STDVectorPointer rVec = m_BasisMatrixGenerator->GetSamplingSchemeRSpace()->GetRadiusVector();
      utlGlobalException(rVec->size()==0 && m_Radius<0, "need to set m_Radius (r value) for EAP");
      if (rVec->size()==0 && m_Radius>=0)
        utl::MatchBVectorAndGradientMatrix(m_Radius, *rVec, *this->m_Orientations);
      }
    else if (m_BasisMatrixGenerator->GetOutputType()==BasisMatrixGeneratorType::ODF)
      {
      m_BasisMatrixGenerator->GetSamplingSchemeRSpace()->SetOrientationsCartesian(this->m_Orientations);
      }

    // m_BasisMatrixGenerator->Print(std::cout<<"m_BasisMatrixGenerator");
    m_BasisMatrixGenerator->Flip(this->m_Flip[0], this->m_Flip[1], this->m_Flip[2]);
    m_BasisMatrixGenerator->ComputeBasisMatrix();
    this->m_BasisMatrix = m_BasisMatrixGenerator->GetBasisMatrix();
    return this->m_BasisMatrix;
    }

protected:
  MeshFromDiscreteFiberODFImageFilter() : Superclass()
    {
    m_BasisMatrixGenerator = NULL;
    m_Radius = -1.0;
    }
  ~MeshFromDiscreteFiberODFImageFilter()
    {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    if (m_BasisMatrixGenerator)
      m_BasisMatrixGenerator->Print(std::cout<<"m_BasisMatrixGenerator");
    }
  
  double m_Radius;
  // virtual void GenerateData();
  // void ThreadedGenerateData(const typename TInputImage::RegionType& regionForThread,ThreadIdType threadId );
  // void BeforeThreadedGenerateData();
  
  
  BasisMatrixGeneratorPointer m_BasisMatrixGenerator;

private:
  MeshFromDiscreteFiberODFImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented


};

} // end namespace itk


#endif
