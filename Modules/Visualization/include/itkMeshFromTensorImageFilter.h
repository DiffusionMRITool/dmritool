/**
 *       @file  itkMeshFromTensorImageFilter.h
 *      @brief  
 *     Created  "06-09-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromTensorImageFilter_h
#define __itkMeshFromTensorImageFilter_h

#include <vtkPolyData.h>
#include <vtkTensorGlyph.h>
#include <vtkDataArray.h>
#include "itkMeshFromImageImageFilter.h"

#include "itkDiffusionTensor.h"

namespace itk
{

/** \class MeshFromTensorImageFilter
 * \brief Compute mesh from a tensor image. 
 *
 * The tensors are assumed to be stored in TENSOR_UPPER_TRIANGULAR format [xx, xy, xz, yy, yz, zz].
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
template < class TInputImage, class TOutputMesh=vtkPolyData>
class ITK_EXPORT MeshFromTensorImageFilter 
  : public MeshFromImageImageFilter<TInputImage, TOutputMesh>
{
public:
  /** Standard class typedefs. */
  typedef MeshFromTensorImageFilter Self;
  typedef MeshFromImageImageFilter<TInputImage, TOutputMesh> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( MeshFromTensorImageFilter, MeshFromImageImageFilter );

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
  typedef typename Superclass::MatrixType                  MatrixType;
  typedef typename Superclass::MatrixPointer               MatrixPointer;
  typedef typename Superclass::VectorType                  VectorType;
  typedef typename Superclass::VectorPointer               VectorPointer;
  typedef typename Superclass::STDVectorType               STDVectorType;
  typedef typename Superclass::STDVectorPointer            STDVectorPointer;
  
  typedef Image<double,3> ScalarImageType;
  typedef typename ScalarImageType::Pointer ScalarImagePointer;

  typedef vtkTensorGlyph  TensorGlyphType;

  /** enum of the possible shapes for the glyphs */
  typedef enum 
  {
    GLYPH_LINE=0,      //0
    GLYPH_ARROW,       //1
    GLYPH_DISK,        //2
    GLYPH_CYLINDER,    //3
    GLYPH_CUBE,        //4
    /** ellipsoild  */
    GLYPH_SPHERE,      //5
    /** superquadric   */
    GLYPH_SUPERQUADRIC //6
  } GlyphShapeType;


  typedef enum
  {
  /** default color  */
    COLOR_NONE=0, 
  /** by principal direction  */
    COLOR_BY_DIRECTION, 
    /** by fa  */
    COLOR_BY_FA,   
    /** by md  */
    COLOR_BY_MD, 
   /** by input scalar image  */
    COLOR_BY_IMAGE
  } TensorColorSchemeType;
  
  itkSetGetMacro(ScalarImage, ScalarImagePointer);

  itkSetGetMacro(TensorColorScheme, TensorColorSchemeType);

  itkSetGetMacro(ShapeMode, GlyphShapeType);
  
  void SetGlyphResolution(int res);
  itkGetMacro(GlyphResolution, int);

  void SetGlyphShape (GlyphShapeType i);
  
  /** Set a specific glyph shape */
  void SetGlyphShapeToLine();
  
  /** Set a specific glyph shape */
  void SetGlyphShapeToDisk();
    
  /** Set a specific glyph shape */
  void SetGlyphShapeToArrow();
  
  /** Set a specific glyph shape */
  void SetGlyphShapeToCube();
  
  /** Set a specific glyph shape */
  void SetGlyphShapeToCylinder();
  
  /** Set a specific glyph shape */
  void SetGlyphShapeToSphere();
  
  /** Set a specific glyph shape */
  void SetGlyphShapeToSuperquadric();

  
  static double GetScalarCodeFromTensor(const DiffusionTensor<double>& tensor, TensorColorSchemeType colorScheme);

  // static std::vector<double> GetRGBAFromTensorDirection(const DiffusionTensor<double>& tensor);

  
protected:
  MeshFromTensorImageFilter(); 

  ~MeshFromTensorImageFilter()
    {};
  
  virtual void GenerateData() ITK_OVERRIDE;

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE
    {
    }
  
  typename LightObject::Pointer InternalClone() const ITK_OVERRIDE
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();

    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }

    rval->m_TensorColorScheme = m_TensorColorScheme;

    return loPtr;
    }

  void VerifyInputParameters() const ITK_OVERRIDE
    {
    if (m_TensorColorScheme==COLOR_BY_IMAGE)
      {
      utlSAGlobalException(itk::IsImageEmpty(m_ScalarImage))(this->m_TensorColorScheme).msg("need to set m_ScalarImage for COLOR_BY_IMAGE");

      InputImageConstPointer inputPtr = this->GetInput();
      utlSAGlobalException(!itk::VerifyImageSize(inputPtr, m_ScalarImage))
        .msg("m_ScalarImage and input need to have the same size");
      }
    }

  
  TensorColorSchemeType m_TensorColorScheme=COLOR_BY_DIRECTION;

  GlyphShapeType m_ShapeMode = GLYPH_SPHERE;

  vtkSmartPointer<vtkPolyDataAlgorithm>  m_Shape = vtkPolyDataAlgorithm::New();

  vtkSmartPointer<TensorGlyphType> m_Glyph = TensorGlyphType::New();

  // vtkSmartPointer<vtkDataArray>  m_Scalars = vtkDataArray::New();

  /** resolution for glyphs  */
  int m_GlyphResolution=10;

  ScalarImagePointer m_ScalarImage = ScalarImageType::New();


private:
  MeshFromTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);//purposely not implemented

};

}

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkMeshFromTensorImageFilter_hxx)
#include "itkMeshFromTensorImageFilter.hxx"
#endif


#endif 
