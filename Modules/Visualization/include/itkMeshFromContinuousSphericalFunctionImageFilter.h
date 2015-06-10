/**
 *       @file  itkMeshFromContinuousSphericalFunctionImageFilter.h
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

#ifndef __itkMeshFromContinuousSphericalFunctionImageFilter_h
#define __itkMeshFromContinuousSphericalFunctionImageFilter_h


#include "itkMeshFromSphericalFunctionImageFilter.h"

namespace itk
{
/** 
*   \defgroup Visualization
*   \brief Visualization of diffusion data 
* */

/** \class MeshFromContinuousSphericalFunctionImageFilter
 * \brief Compute mesh from continuous spherical function represented by a set of basis.
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromContinuousSphericalFunctionImageFilter 
  : public MeshFromSphericalFunctionImageFilter<TInputImage, TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromContinuousSphericalFunctionImageFilter Self;
  typedef MeshFromSphericalFunctionImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromContinuousSphericalFunctionImageFilter, MeshFromSphericalFunctionImageFilter);

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

  
  itkSetMacro(BasicShape, BasicShapeType);
  itkGetMacro(BasicShape, BasicShapeType);

  itkSetMacro(TessellationOrder, unsigned int);
  itkGetMacro(TessellationOrder, unsigned int);
  
  itkSetMacro(BasisMatrix, MatrixPointer);
  itkGetMacro(BasisMatrix, MatrixPointer);

  virtual MatrixPointer ComputeBasisMatrix()
    {
    return MatrixPointer();
    }
  
protected:
  MeshFromContinuousSphericalFunctionImageFilter(); 

  ~MeshFromContinuousSphericalFunctionImageFilter()
    {};

  void ScaleSamples(VectorType& b) const;

  void PrintSelf(std::ostream& os, Indent indent) const;
  
  // virtual void GenerateData();
  void ThreadedGenerateData(const typename TInputImage::RegionType& regionForThread,ThreadIdType threadId );
  void BeforeThreadedGenerateData();

  SphereTessellatorPointer m_SphereTessellator;
  BasicShapeType m_BasicShape;
  unsigned int m_TessellationOrder;

  MatrixPointer m_BasisMatrix;
  

private:
  MeshFromContinuousSphericalFunctionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_MeshFromContinuousSphericalFunctionImageFilter(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT MeshFromContinuousSphericalFunctionImageFilter< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef MeshFromContinuousSphericalFunctionImageFilter< ITK_TEMPLATE_2 x > MeshFromContinuousSphericalFunctionImageFilter##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkMeshFromContinuousSphericalFunctionImageFilter+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkMeshFromContinuousSphericalFunctionImageFilter_hxx)
#include "itkMeshFromContinuousSphericalFunctionImageFilter.hxx"
#endif

#endif
