/**
 *       @file  itkMeshFromSphericalFunctionImageFilter.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-18-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromSphericalFunctionImageFilter_h
#define __itkMeshFromSphericalFunctionImageFilter_h

#include "itkSphereTessellator.h"
#include "itkMeshFromImageImageFilter.h"

namespace itk
{

/** \class MeshFromSphericalFunctionImageFilter
 * \brief Compute mesh from general spherical function.
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromSphericalFunctionImageFilter 
  : public MeshFromImageImageFilter<TInputImage, TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromSphericalFunctionImageFilter Self;
  typedef MeshFromImageImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromSphericalFunctionImageFilter, MeshFromImageImageFilter );

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

  /** Normalization type */
  typedef enum {NONE=0, MIN_MAX, UNIT_MAX, UNIT_INTEGRAL} NormalizationType;
  /** Set/Get whether and how to do normalization */
  itkSetMacro(Normalization, NormalizationType);
  itkGetMacro(Normalization, NormalizationType);

  itkSetMacro(Orientations, MatrixPointer);
  itkGetMacro(Orientations, MatrixPointer);
  
  /** Set/Get pow factor  */
  itkSetMacro(Pow, double);
  itkGetMacro(Pow, double);
  
  /** Set/Get whether to remove negative values  */
  itkSetMacro(RemoveNegativeValues, bool);
  itkGetMacro(RemoveNegativeValues, bool);
  itkBooleanMacro(RemoveNegativeValues);

  unsigned int GetNumberOfOrientations() const
    {
    return m_Orientations->Rows();
    }

  
protected:
  MeshFromSphericalFunctionImageFilter() : Superclass(), 
    m_Orientations(new MatrixType())
  {
  m_Pow = 1.0;
  m_RemoveNegativeValues = false;
  m_Normalization = NONE;
  this->m_RemoveNegativeValues = false;

  this->m_ColorScheme = Superclass::DIRECTION;
  }

  ~MeshFromSphericalFunctionImageFilter()
    {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    PrintVar(true, os<<indent, m_Pow, m_RemoveNegativeValues);
    PrintEnum4(true, m_Normalization, NONE, MIN_MAX, UNIT_MAX, UNIT_INTEGRAL, os<<indent);
    utl::PrintUtlMatrix(*m_Orientations, "m_Orientations", " ", os<<indent);
    }
  
  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();

    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }

    rval->m_Pow = m_Pow;
    rval->m_RemoveNegativeValues = m_RemoveNegativeValues;
    rval->m_Orientations = m_Orientations;
    rval->m_Normalization = m_Normalization;

    return loPtr;
    }
  
  // virtual void GenerateData();
  // void ThreadedGenerateData(const typename TInputImage::RegionType& regionForThread,ThreadIdType threadId );
  // void BeforeThreadedGenerateData();

  bool m_RemoveNegativeValues;

  double m_Pow;

  MatrixPointer m_Orientations;
  
  NormalizationType m_Normalization;

private:
  MeshFromSphericalFunctionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented

};

} // end namespace itk


#endif
