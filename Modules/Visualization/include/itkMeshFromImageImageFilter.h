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

// #include "vnl/vnl_matrix.h"
// #include "vnl/vnl_math.h"
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

#include "utlVector.h"
#include "utlMatrix.h"
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

  
  /** Set/Get whether to remove negative values  */
  itkSetMacro(RemoveNegativeValues, bool);
  itkGetMacro(RemoveNegativeValues, bool);
  itkBooleanMacro(RemoveNegativeValues);

  void SetBoxView(const int x0, const int x1, const int y0, const int y1, const int z0, const int z1 )
    {
    m_BoxView[0]=x0;
    m_BoxView[1]=x1;
    m_BoxView[2]=y0;
    m_BoxView[3]=y1;
    m_BoxView[4]=z0;
    m_BoxView[5]=z1;
    }
  std::vector<int> GetBoxView()
    {
    return m_BoxView;
    }

  /** Set/Get scale factor */
  itkSetMacro(Scale, double);
  itkGetMacro(Scale, double);
  
  /** Set/Get pow factor  */
  itkSetMacro(Pow, double);
  itkGetMacro(Pow, double);


  /** Access to Mesh */
  OutputMeshPolyDataPointer GetOutput()
  {
    return m_Mesh;
  }
  
protected:
  MeshFromImageImageFilter()
    {
    this->SetNumberOfRequiredInputs(1);
    m_RemoveNegativeValues = false;
    m_Scale = 1.0;
    m_Pow = 1.0;
    m_ColorScheme = UNKNOWN;
    m_BoxView = std::vector<int>(6,-1);
    m_Mesh = OutputMeshPolyDataType::New();
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
    return true;
    }

  virtual void VerifyInputParameters() const
    {
    }

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    PrintVar3(true, m_Scale, m_Pow, m_RemoveNegativeValues, os<<indent);
    utl::PrintVector(m_BoxView, "m_BoxView");
    PrintVar1(true, m_ColorScheme, os<<indent);
    m_Mesh->Print(std::cout<<"m_Mesh=");
    }

  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();

    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_RemoveNegativeValues = m_RemoveNegativeValues;
    rval->m_Scale = m_Scale;
    rval->m_Pow = m_Pow;
    rval->m_ColorScheme = m_ColorScheme;
    rval->m_BoxView = m_BoxView;
    return loPtr;
    }


  /** the information of the output image is used for ThreadedGenerateData() in subclass */
  virtual void GenerateOutputInformation()
    {
    typename TInputImage::Pointer outputPtr = itkDynamicCastInDebugMode< TInputImage * >( this->GetPrimaryOutput() );
    InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput());
    itk::CopyImageInformation(inputPtr, outputPtr);
    outputPtr->SetNumberOfComponentsPerPixel(1);
    }

  /** The output is m_Mesh, do not allocate output image  */
  virtual void AllocateOutputs()
    {
    }

  OutputMeshPolyDataPointer m_Mesh;

  bool m_RemoveNegativeValues;

  double m_Scale;
  double m_Pow;

  ColorSchemeType m_ColorScheme;

  std::vector<int> m_BoxView;

private:
  MeshFromImageImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented
};

} // end namespace itk


#endif
