/**
 *       @file  itkMeshFromImageImageFilter.h
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

#ifndef __itkMeshFromImageImageFilter_h
#define __itkMeshFromImageImageFilter_h

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
#include "vtkSmartPointer.h"

#include "utlNDArray.h"
#include "itkImageSource.h"
#include "itkMaskedImageToImageFilter.h"
#include "utlCore.h"

namespace itk
{

/** \class MeshFromImageImageFilter
 * \brief Compute mesh from SH coefficients.
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromImageImageFilter 
  : public MaskedImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromImageImageFilter Self;
  typedef MaskedImageToImageFilter<TInputImage, TInputImage> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromImageImageFilter, MaskedImageToImageFilter );

  /** Orienation dimension */
  static const unsigned int OrientationDimension = 3;

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

  typedef vtkPoints OutputMeshPointsType;
  typedef vtkCellArray OutputMeshCellArrayType;
  // typedef vtkPolyData OutputMeshPolyDataType;
  typedef vtkSmartPointer<OutputMeshPolyDataType> OutputMeshPolyDataPointer;
  typedef vtkDoubleArray OutputMeshScalarType;
  typedef vtkUnsignedCharArray OutputMeshRGBType;

  /** Orientation Matrices */
  typedef utl::NDArray<double,2>                  MatrixType;
  typedef utl::NDArray<double,1>                  VectorType;
  typedef utl_shared_ptr<MatrixType>              MatrixPointer;
  typedef utl_shared_ptr<VectorType>              VectorPointer;
  typedef std::vector<double>                     STDVectorType;
  typedef utl_shared_ptr<STDVectorType >          STDVectorPointer;
  
  typedef  enum {UNKNOWN=0, FIXED, DIRECTION, MAGNITUDE} ColorSchemeType;

  itkSetMacro(ColorScheme, ColorSchemeType);
  itkGetMacro(ColorScheme, ColorSchemeType);
  
  // virtual void SetInput(const InputImageType *image)
  //   {
  //   this->ProcessObject::SetNthInput( 0,const_cast< InputImageType * >( image ) );
  //   }

  // const InputImageType * GetInput(void) const
  //   {
  //   return itkDynamicCastInDebugMode< const TInputImage * >( this->GetPrimaryInput() );
  //   }


  void SetBoxView(const int x0, const int x1, const int y0, const int y1, const int z0, const int z1 )
    {
    m_BoxView[0]=x0;
    m_BoxView[1]=x1;
    m_BoxView[2]=y0;
    m_BoxView[3]=y1;
    m_BoxView[4]=z0;
    m_BoxView[5]=z1;
    }
  std::vector<int> GetBoxView() const
    {
    return m_BoxView;
    }

  void SetSliceView(const int coronal, const int sagittal, const int transverse )
    {
    m_SliceView[0]=coronal, m_SliceView[1]=sagittal, m_SliceView[2]=transverse;
    }
  std::vector<int> GetSliceView() const
    {
    return m_SliceView;
    }

  void SetFlip(const int flipx, const int flipy, const int flipz )
    {
    m_Flip[0]=flipx, m_Flip[1]=flipy, m_Flip[2]=flipz;
    }
  std::vector<int> GetFlip() const
    {
    return m_Flip;
    }

  /** Set/Get scale factor */
  itkSetMacro(Scale, double);
  itkGetMacro(Scale, double);


  /** Access to Mesh */
  OutputMeshPolyDataPointer GetOutput()
  {
    return m_Mesh;
  }
  
protected:
  MeshFromImageImageFilter()
    {
    this->SetNumberOfRequiredInputs(1);
    // TODO: multi-thread
    this->SetNumberOfThreads(1);
    }
  ~MeshFromImageImageFilter()
    { 
    if (m_Mesh) 
      m_Mesh->Delete(); 
    }

  virtual bool IsPixelIndexVisible(const InputImageIndexType& index )
    {
    if (this->IsMaskUsed() && this->m_MaskImage->GetPixel(index)==0)
      return false;

    if ( (m_BoxView[0]>=0 || m_BoxView[1]>=0) && (index[0]<m_BoxView[0] || index[0]>m_BoxView[1]) )
      return false;
    if ( (m_BoxView[2]>=0 || m_BoxView[3]>=0) && (index[1]<m_BoxView[2] || index[1]>m_BoxView[3]) )
      return false;
    if ( (m_BoxView[4]>=0 || m_BoxView[5]>=0) && (index[2]<m_BoxView[4] || index[2]>m_BoxView[5]) )
      return false;

    if ( (m_SliceView[0]>=0 || m_SliceView[1]>=0 || m_SliceView[2]>=0) && index[0]!=m_SliceView[0] && index[1]!=m_SliceView[1] && index[2]!=m_SliceView[2])
      return false;

    return true;
    }

  virtual void VerifyInputParameters() const ITK_OVERRIDE
    {
    }

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE
    {
    PrintVar(true, os<<indent, m_Scale );
    utl::PrintVector(m_BoxView, "m_BoxView");
    utl::PrintVector(m_SliceView, "m_SliceView");
    utl::PrintVector(m_Flip, "m_Flip");
    PrintVar(true, os<<indent, m_ColorScheme);
    m_Mesh->Print(std::cout<<"m_Mesh=");
    }

  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();

    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_Scale = m_Scale;
    rval->m_ColorScheme = m_ColorScheme;

    rval->m_BoxView = m_BoxView;
    rval->m_SliceView = m_SliceView;
    
    rval->m_Flip = m_Flip;

    return loPtr;
    }


  /** the information of the output image is used for ThreadedGenerateData() in subclass */
  virtual void GenerateOutputInformation() ITK_OVERRIDE
    {
    typename TInputImage::Pointer outputPtr = itkDynamicCastInDebugMode< TInputImage * >( this->GetPrimaryOutput() );
    InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput());
    itk::CopyImageInformation(inputPtr, outputPtr);
    outputPtr->SetNumberOfComponentsPerPixel(1);
    }

  /** The output is m_Mesh, do not allocate output image  */
  virtual void AllocateOutputs() ITK_OVERRIDE
    {
    }

  OutputMeshPolyDataPointer m_Mesh = OutputMeshPolyDataType::New();

  double m_Scale=1.0;

  ColorSchemeType m_ColorScheme=UNKNOWN;

  /** 3D box view (xmin, xmax, ymin, ymax, zmin, zmax)  */
  std::vector<int> m_BoxView = std::vector<int>(6,-1);

  /** slice view (x,y,z), 3 orthogonal slices  */
  std::vector<int> m_SliceView = std::vector<int>(3,-1);
  
  /** flips in x/y/z-axis. 0 means no flip, 1 means flip  */
  std::vector<int> m_Flip = std::vector<int>(3, 0);

private:
  MeshFromImageImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented
};

} // end namespace itk


#endif
