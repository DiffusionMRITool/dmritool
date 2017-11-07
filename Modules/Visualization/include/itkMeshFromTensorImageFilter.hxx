/**
 *       @file  itkMeshFromTensorImageFilter.hxx
 *      @brief  
 *     Created  "06-09-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */



#ifndef __itkMeshFromTensorImageFilter_hxx
#define __itkMeshFromTensorImageFilter_hxx

#include "itkMeshFromTensorImageFilter.h"

#include <vtkObjectFactory.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkStructuredPoints.h>
#include <vtkLookupTable.h>
#include <vtkFieldData.h>
#include <vtkArrowSource.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkDiskSource.h>
#include <vtkLineSource.h>
#include <vtkSuperquadricSource.h>
#include <vtkProperty.h>

#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
template <class TInputImage, class TOutputMesh>
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::MeshFromTensorImageFilter() : Superclass()
{

  this->m_Glyph = vtkTensorGlyph::New();
  this->m_Glyph->ClampScalingOn();
  this->m_Glyph->ColorGlyphsOn();

  // default shape is ellipsoid
  this->SetGlyphShapeToSphere();
  this->m_Scale = 1000; // scale for ellipsoids. 
  this->m_Glyph->SetScaleFactor(this->m_Scale);

  // this->m_LUT = 0;
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShape(GlyphShapeType i)
{
  switch(i)
  {
      case GLYPH_LINE:
        this->SetGlyphShapeToLine();
        break;

      case GLYPH_ARROW:
        this->SetGlyphShapeToArrow();
        break;

      case GLYPH_DISK:
        this->SetGlyphShapeToDisk();
        break;

      case GLYPH_CYLINDER:
        this->SetGlyphShapeToCylinder();
        break;

      case GLYPH_CUBE:
        this->SetGlyphShapeToCube();
        break;

      case GLYPH_SPHERE:
        this->SetGlyphShapeToSphere();
        break;

      case GLYPH_SUPERQUADRIC:
        this->SetGlyphShapeToSuperquadric();
        break;

      default:
        utlGlobalException(true, "shape type is not recognized.");
        return;
  }

  this->SetGlyphResolution( this->m_GlyphResolution );
}


template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToLine()
{
  this->m_ShapeMode = GLYPH_LINE;
  this->m_Shape = vtkLineSource::New();
  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToDisk()
{
  this->m_ShapeMode = GLYPH_DISK;
  vtkDiskSource* source = vtkDiskSource::New();
  source->SetInnerRadius (0.0);
  this->m_Shape = source;
  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToArrow()
{
  this->m_ShapeMode = GLYPH_ARROW;
  this->m_Shape = vtkArrowSource::New();
  vtkArrowSource* arrow = vtkArrowSource::SafeDownCast (this->m_Shape);
  // arrow->SetTipLength (0.0);
  // arrow->SetTipRadius (0.0);
  // arrow->SetShaftRadius (0.18);

  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToCube()
{
  this->m_ShapeMode = GLYPH_CUBE;
  this->m_Shape = vtkCubeSource::New();
  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToCylinder()
{
  this->m_ShapeMode = GLYPH_CYLINDER;
  this->m_Shape = vtkCylinderSource::New();
  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToSphere()
{
  this->m_ShapeMode = GLYPH_SPHERE;
  this->m_Shape = vtkSphereSource::New();
  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphShapeToSuperquadric()
{
  this->m_ShapeMode = GLYPH_SUPERQUADRIC;
  this->m_Shape = vtkSuperquadricSource::New();
  vtkSuperquadricSource::SafeDownCast (this->m_Shape)->SetPhiRoundness (0.3);
  vtkSuperquadricSource::SafeDownCast (this->m_Shape)->SetThetaRoundness (0.3);
  this->m_Glyph->SetSourceConnection(this->m_Shape->GetOutputPort());
  this->m_Shape->Delete();
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::SetGlyphResolution (int res)
{

  this->m_GlyphResolution = res;

  vtkArrowSource* arrowSource = 0;
  vtkDiskSource* diskSource = 0;
  vtkCylinderSource* cylinderSource = 0;
  vtkSphereSource* sphereSource = 0;
  vtkSuperquadricSource* quadricSource = 0;

  switch(this->m_ShapeMode)
    {

  case GLYPH_LINE:
    break;

  case GLYPH_ARROW:
    arrowSource = vtkArrowSource::SafeDownCast (this->m_Shape);
    if( arrowSource )
      arrowSource->SetShaftResolution (res);
    break;

  case GLYPH_DISK:
    diskSource = vtkDiskSource::SafeDownCast (this->m_Shape);
    if( diskSource )
      diskSource->SetCircumferentialResolution(res);
    break;

  case GLYPH_CYLINDER:
    cylinderSource = vtkCylinderSource::SafeDownCast (this->m_Shape);
    if( cylinderSource )
      cylinderSource->SetResolution (res);
    break;

  case GLYPH_CUBE:
    break;

  case GLYPH_SPHERE:
    sphereSource = vtkSphereSource::SafeDownCast (this->m_Shape);
    if( sphereSource )
      {
      sphereSource->SetThetaResolution (res);
      sphereSource->SetPhiResolution (res);
      }
    break;

  case GLYPH_SUPERQUADRIC:
    quadricSource = vtkSuperquadricSource::SafeDownCast (this->m_Shape);
    if( quadricSource )
      {
      quadricSource->SetThetaResolution (res);
      quadricSource->SetPhiResolution (res);
      }
    break;

  default:
    break;

    }

  this->m_Glyph->Modified();
}


template <class TInputImage, class TOutputMesh>
double
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::GetScalarCodeFromTensor(const DiffusionTensor<double>& tensor, TensorColorSchemeType colorScheme)
{
  if (colorScheme==COLOR_BY_FA)
    {
    return tensor.GetFA();
    }
  else if (colorScheme==COLOR_BY_MD)
    {
    return tensor.GetMD();
    }
  else if (colorScheme==COLOR_BY_DIRECTION)
    {
    double fa = tensor.GetFA();
    // if fa is close to 0, then the eigenvector is not reliable. Thus, we set ss=0, if fa=0.
    if (fa<0.01)
      return 0;

    VectorType eigenValues(3), v1(3); 
    MatrixType eigenVectors(3,3);

    tensor.GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
    for ( int i = 0; i < 3; ++i ) 
      v1[i] = eigenVectors(i,2);

    double ss = 0;
    utl::RGBToIndex(std::fabs(v1[0]),std::fabs(v1[1]),std::fabs(v1[2]),ss);
    return ss;
    }
  else
    {
    utlGlobalException(true, "Logical error! Wrong colorScheme");
    return 0.0;
    }
}

// template <class TInputImage, class TOutputMesh>
// std::vector<double>
// MeshFromTensorImageFilter<TInputImage, TOutputMesh>
// ::GetRGBAFromTensorDirection(const DiffusionTensor<double>& tensor)
// {

//   // std::vector<double> rgba(3);
//   // double fa = tensor.GetFA();
//   // // if fa is close to 0, then the eigenvector is not reliable. Thus, we set ss=0, if fa=0.
//   // if (fa<0.01)
//   //   {
//   //   rgba[0]=0, rgba[1]=0, rgba[2]=0;
//   //   return rgba;
//   //   }

//   // VectorType eigenValues(3), v1(3); 
//   // MatrixType eigenVectors(3,3);

//   // tensor.GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
//   // for ( int i = 0; i < 3; ++i ) 
//   //   v1[i] = eigenVectors(i,2);


//   // for (unsigned int d=0;d<3;d++)
//   //   {
//   //   rgba[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>(std::fabs(v1[d])*255.0);
//   //   }

//   // double min_val = std::numeric_limits<double>::max();
//   // double max_val = -std::numeric_limits<double>::max();

//   // for (unsigned int d=0;d<3;d++)
//   //   {
//   //   if (rgba[d] > max_val)
//   //     max_val = rgba[d];
//   //   if (rgba[d] < min_val)
//   //     min_val = rgba[d];
//   //   }

//   // // min_val=0;
//   // for (unsigned int d=0;d<3;d++)
//   //   {
//   //   rgba[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>( (rgba[d] - min_val)/(max_val - min_val) * 255.0 );
//   //   }
      

//   double ss = GetScalarCodeFromTensor(tensor, COLOR_BY_DIRECTION);

//   vtkSmartPointer<vtkLookupTable> lut = vtkLookupTable::New();
//   lut->SetHueRange(0.0,1.0);
//   lut->SetRampToLinear();
//   lut->SetTableRange(0,255.0);
//   lut->Build();

//   // unsigned char* tmp = lut->MapValue(ss);
//   std::vector<double> rgba(3);
//   // rgba[0]=tmp[0], rgba[1]=tmp[1], rgba[2]=tmp[2], rgba[3]=tmp[3];
//   double tmp[3];
//   lut->GetColor(ss, tmp);
//   rgba[0]=tmp[0]*255, rgba[1]=tmp[1]*255, rgba[2]=tmp[2]*255;
//   utlPrintVar(true, "a1", ss, rgba[0], rgba[1], rgba[2]);

//   return rgba;
// }

template <class TInputImage, class TOutputMesh>
void 
MeshFromTensorImageFilter<TInputImage, TOutputMesh>
::GenerateData()
{
  utlShowPosition(this->GetDebug());

  this->VerifyInputParameters();

  this->SetGlyphShape(m_ShapeMode);

  InputImageConstPointer inputPtr = this->GetInput();
  InputImagePixelType inputPixel, vec9d;
  inputPixel.SetSize(inputPtr->GetNumberOfComponentsPerPixel());
  vec9d.SetSize(9);

  std::vector<int> size3d = itk::GetVectorImage3DVolumeSize(inputPtr);
  auto spacing = inputPtr->GetSpacing();
  auto origin = inputPtr->GetOrigin();

  vtkSmartPointer<vtkPoints> tens_points = vtkPoints::New();
  vtkSmartPointer<vtkFloatArray> tens_array = vtkFloatArray::New();
  tens_array->SetNumberOfComponents(9);
  tens_array->SetNumberOfTuples(size3d[0]*size3d[1]*size3d[2]);
    
  vtkSmartPointer<vtkFloatArray> scalars = vtkFloatArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetName("scalar_color");

  vtkSmartPointer<vtkUnsignedCharArray> rgbaArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  rgbaArray->SetNumberOfComponents(3);
  rgbaArray->SetName("RGBA_color");
  double rgba[3]={0,0,0};

  InputImagePointType inputPhysicalPoint;

  ImageRegionConstIteratorWithIndex<InputImageType> inputIt(inputPtr, inputPtr->GetLargestPossibleRegion());
  ImageRegionConstIteratorWithIndex<ScalarImageType> scalarIt;
  if (!itk::IsImageEmpty(m_ScalarImage))
    scalarIt = ImageRegionConstIteratorWithIndex<ScalarImageType>(m_ScalarImage, m_ScalarImage->GetLargestPossibleRegion());
  
  DiffusionTensor<double> tensor;

  long off=0;
  for (inputIt.GoToBegin(), scalarIt.GoToBegin(); 
    !inputIt.IsAtEnd(); 
    ++inputIt, ++scalarIt) 
    {

    auto inputIndex = inputIt.GetIndex();

    if (!this->IsPixelIndexVisible(inputIndex))
      {
      continue;
      }

    auto inputVec = inputIt.Get();
    if (inputVec.GetNorm()< 1e-10)
      {
      continue;
      }

    if (this->GetDebug())
      {
      utlPrintVar(true, inputIndex, off);
      itk::PrintVariableLengthVector(inputVec, "inputVec");
      }

    inputPtr->TransformIndexToPhysicalPoint(inputIndex, inputPhysicalPoint);

    tens_points->InsertPoint(off, inputPhysicalPoint.GetDataPointer());

    for ( int i = 0; i < 6; ++i ) 
      tensor[i] = inputVec[i];
    tensor.Flip(this->m_Flip[0], this->m_Flip[1], this->m_Flip[2]);

    utl::ConvertTensor6DTo9D(tensor, vec9d, TENSOR_UPPER_TRIANGULAR);
    tens_array->InsertTuple9(off,vec9d[0],vec9d[1],vec9d[2],vec9d[3],vec9d[4],vec9d[5],vec9d[6],vec9d[7],vec9d[8]);
    
    if (m_TensorColorScheme!=COLOR_NONE)
      {
      double ss=0;
      if (m_TensorColorScheme!=COLOR_BY_IMAGE && m_TensorColorScheme!=COLOR_BY_DIRECTION)
        {
        ss = GetScalarCodeFromTensor(tensor, m_TensorColorScheme);
        scalars->InsertTuple1(off, ss);
        }
      else if (m_TensorColorScheme==COLOR_BY_IMAGE)
        {
        ss = scalarIt.Get();
        scalars->InsertTuple1(off, ss);
        }
      else if (m_TensorColorScheme==COLOR_BY_DIRECTION)
        {
        ss = GetScalarCodeFromTensor(tensor, m_TensorColorScheme);
        scalars->InsertTuple1(off, ss);

        // auto colorVec = GetRGBAFromTensorDirection(tensor);
        // utl::PrintVector(colorVec, "colorVec");
        // for ( int i = 0; i < colorVec.size(); ++i ) 
        //   rgba[i] = colorVec[i];
        // rgbaArray->InsertTuple(off, rgba);
        }

      if (this->GetDebug())
        utlPrintVar(true, ss);

      // if (this->GetDebug())
      //   utlPrintVar(true, ss, rgba[0], rgba[1], rgba[2], rgba[3]);
      }

    off++;
    }

  this->m_Mesh->SetPoints(tens_points);
  this->m_Mesh->GetPointData()->SetTensors(tens_array);
  if (m_TensorColorScheme!=COLOR_NONE)
    {
    // if (m_TensorColorScheme==COLOR_BY_DIRECTION)
    //   this->m_Mesh->GetPointData()->SetScalars(rgbaArray);
    // else
    //   this->m_Mesh->GetPointData()->SetScalars(scalars);
      
    this->m_Mesh->GetPointData()->SetScalars(scalars);
    }

  m_Glyph->SetInputData(this->m_Mesh);
  m_Glyph->SetScaleFactor(this->m_Scale);

  m_Glyph->ColorGlyphsOn();
  m_Glyph->ClampScalingOn();

  m_Glyph->Update();

  this->m_Mesh = m_Glyph->GetOutput();

}


}


#endif 
