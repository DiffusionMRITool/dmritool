/**
 *       @file  itkMeshFromPeaksImageFilter.h
 *      @brief  
 *     Created  "08-26-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromPeaksImageFilter_h
#define __itkMeshFromPeaksImageFilter_h


#include "itkPoint.h"
#include "itkPointSet.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkTypeTraits.h"
#include "vtkType.h"

#include "utlITKMacro.h"
#include "utlCoreMacro.h"
#include "itkMeshFromImageImageFilter.h"
#include "itkPeakContainerHelper.h"

namespace itk
{

/** \class MeshFromPeaksImageFilter
 * \brief Generate a mesh from given peaks.
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData >
class ITK_EXPORT MeshFromPeaksImageFilter :
  public MeshFromImageImageFilter<TInputImage, TOutputMesh>
{
public:
  typedef MeshFromPeaksImageFilter Self;
  typedef MeshFromImageImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromPeaksImageFilter, MeshSource );

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
  
  /** Orientation Matrices */
  typedef utl::NDArray<double,2>                  MatrixType;
  typedef utl::NDArray<double,1>                  VectorType;
  typedef utl_shared_ptr<MatrixType>              MatrixPointer;
  typedef utl_shared_ptr<VectorType>              VectorPointer;
  typedef std::vector<double>                     STDVectorType;
  typedef utl_shared_ptr<STDVectorType >          STDVectorPointer;

  itkSetGetMacro(PeakType, PeakType);
  itkSetGetMacro(TubeRadius, double);
  itkSetGetMacro(MaxNumberOfPeaks, int);

  void SetColorPeak(double r, double g, double b)
    {
    m_ColorPeak[0]=r, m_ColorPeak[1]=g, m_ColorPeak[2]=b;
    }

  
protected:
  MeshFromPeaksImageFilter();
  ~MeshFromPeaksImageFilter()
    {  };
  
  void VerifyInputParameters() const
    {
    utlSAException(this->m_ColorScheme!=Superclass::FIXED && this->m_ColorScheme!=Superclass::DIRECTION)
      (this->m_ColorScheme).msg("wrong m_ColorScheme");
    }

  typename LightObject::Pointer InternalClone() const;

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    std::string peaktypeStr = PeakContainerHelper::GetString(m_PeakType);
    PrintVar(true, os<<indent, peaktypeStr, m_TubeRadius, m_MaxNumberOfPeaks);
    PrintVar(true, os<<indent, m_ColorPeak[0], m_ColorPeak[1], m_ColorPeak[2]);
    }

  virtual void GenerateData();

  PeakType m_PeakType;

  /** If it is positive, use tubes for peaks. If it is non-positive, use lines. */
  double m_TubeRadius;

  double m_ColorPeak[3];

  /** If m_MaxNumberOfPeaks>0, then only visualize the first m_MaxNumberOfPeaks peaks.  */
  int m_MaxNumberOfPeaks;

private:
  MeshFromPeaksImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_MeshFromPeaksImageFilter(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT MeshFromPeaksImageFilter< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef MeshFromPeaksImageFilter< ITK_TEMPLATE_2 x > MeshFromPeaksImageFilter##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkMeshFromPeaksImageFilter+-.h"
#endif

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshFromPeaksImageFilter.hxx"
#endif

#endif 
