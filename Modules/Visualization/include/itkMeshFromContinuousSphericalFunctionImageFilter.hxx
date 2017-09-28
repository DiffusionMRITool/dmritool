/**
 *       @file  itkMeshFromContinuousSphericalFunctionImageFilter.hxx
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

#ifndef __itkMeshFromContinuousSphericalFunctionImageFilter_hxx
#define __itkMeshFromContinuousSphericalFunctionImageFilter_hxx
#include "itkMeshFromContinuousSphericalFunctionImageFilter.h"

#include "utl.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vtkColorTransferFunction.h"

#include <vtkFloatArray.h>
#include <itkProgressReporter.h>

namespace itk
{

/**
 * Default constructor.
 */
template<class TInputImage, class TOutputMesh>
MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
::MeshFromContinuousSphericalFunctionImageFilter() : Superclass(),
    m_BasisMatrix(new MatrixType())
{
  m_TessellationOrder = 3;
  m_BasicShape = SphereTessellatorType::ICOSAHEDRON;
  m_SphereTessellator = SphereTessellatorType::New();
}

/**
 * Standard PrintSelf method.
 */
template<class TInputImage, class TOutputMesh>
void
MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "BasicShape: " << m_BasicShape << std::endl;
  os << indent << "TessellationOrder: " << m_TessellationOrder << std::endl;
  utl::PrintUtlMatrix(*m_BasisMatrix, "m_BasisMatrix", " ", os<<indent);
}

template <class TInputImage, class TOutputMesh>
typename LightObject::Pointer
MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_SphereTessellator = m_SphereTessellator;
  rval->m_BasicShape = m_BasicShape;
  rval->m_TessellationOrder = m_TessellationOrder;
  rval->m_BasisMatrix = m_BasisMatrix;

  return loPtr;
}

template <class TInputImage, class TOutputMesh>
void 
MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
::BeforeThreadedGenerateData()
{
  utlShowPosition(this->GetDebug());
  InputImageConstPointer inputPtr = this->GetInput();

  this->VerifyInputParameters();

  this->m_SphereTessellator->SetBasicShape(this->m_BasicShape);
  this->m_SphereTessellator->SetOrder(this->m_TessellationOrder);
  this->m_SphereTessellator->TessellateSphere();
  *this->m_Orientations = utl::VnlMatrixToUtlMatrix(this->m_SphereTessellator->GetPointsMatrix());
  // utl::PrintUtlMatrix(*this->m_Orientations, "m_Orientations 0");

  // this->Print(std::cout<<"this");
  this->ComputeBasisMatrix();

  if (this->GetDebug())
    utl::PrintUtlMatrix(*this->m_BasisMatrix, "this->m_BasisMatrix");
  utlGlobalException(this->m_BasisMatrix->Columns()<=0, "no m_BasisMatrix after run ComputeBasisMatrix()");
}

/** Generate data. */
template <class TInputImage, class TOutputMesh>
void 
MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
// ::GenerateData()
::ThreadedGenerateData(const typename TInputImage::RegionType& regionForThread,ThreadIdType threadId ) 
{
  utlShowPosition(this->GetDebug());
  ProgressReporter progress(this, threadId, regionForThread.GetNumberOfPixels());

  // Pointers
  InputImageConstPointer inputPtr = this->GetInput();
  InputImagePixelType inputPixel;
  inputPixel.SetSize(inputPtr->GetNumberOfComponentsPerPixel());

  // Compute Number Of Basis Functions
  unsigned int numberOfBasis = this->m_BasisMatrix->Columns();

  // Input Image Parameters
  InputImageIndexType inputIndex;
  InputImagePointType inputPhysicalPoint;
  
  // iterator for the input image
  ImageRegionConstIteratorWithIndex<InputImageType> inputIt(inputPtr, regionForThread);
  // ImageRegionConstIteratorWithIndex<InputImageType> inputIt(inputPtr, inputPtr->GetLargestPossibleRegion());
    
  // Preparation work for vertices and cells
  unsigned int numberOfPoints = 0;
  unsigned int numberOfCells = 0;
  vnl_matrix<unsigned long> cellMatrix;

  unsigned long offset = 0;
  
  numberOfPoints = this->m_SphereTessellator->GetNumberOfVertices();
  numberOfCells = this->m_SphereTessellator->GetNumberOfFaces();
  cellMatrix = this->m_SphereTessellator->GetCellsMatrix();
      
  // Matrices
  
  VectorType b;
  VectorType x(numberOfBasis);

  // Output Mesh
  vtkSmartPointer<OutputMeshPointsType> outputMeshPoints = vtkSmartPointer<OutputMeshPointsType>::New();
  vtkSmartPointer<OutputMeshCellArrayType> outputMeshCellArray = vtkSmartPointer<OutputMeshCellArrayType>::New();
  vtkSmartPointer<OutputMeshRGBType> outputMeshRGB = vtkSmartPointer<OutputMeshRGBType>::New();
  outputMeshRGB->SetNumberOfComponents(4);
  outputMeshRGB->SetName("RGBA_scalars");

  double outputMeshPoint[3];
  unsigned int numberOfCellPoints = 3;

  vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = 
    vtkSmartPointer<vtkColorTransferFunction>::New();
  //double colorId = 0;
  double rgb[3];
  colorTransferFunction->SetColorSpaceToHSV();
  colorTransferFunction->AddRGBPoint(0.0, 0, 0, 1);
  colorTransferFunction->AddRGBPoint(0.5, 0, 1, 0);
  colorTransferFunction->AddRGBPoint(1.0, 1, 0, 0);
  
  // Populate output mesh
  offset = 0;
  unsigned char RGB[4];
  
  // vtkSmartPointer<vtkFloatArray> newScalars= vtkSmartPointer<vtkFloatArray>::New();
  // newScalars->SetNumberOfComponents(1);

  double bMax= -std::numeric_limits<double>::max();
  if ( this->m_ColorScheme == Superclass::MAGNITUDE )
    {
    inputIt.GoToBegin();
    while( !inputIt.IsAtEnd() )
      {  
      inputIndex = inputIt.GetIndex();
      if (!this->IsPixelIndexVisible(inputIndex))
        {
        ++inputIt;
        continue;
        }

      inputPixel = inputIt.Get(); 

      for (unsigned int k=0;k<numberOfBasis;k++) //isotropic component?
        x(k) = inputPixel[k];

      if (x.GetRootMeanSquares() == 0)
        {
        ++inputIt;
        continue;
        }

      // if (this->GetDebug())
      //   {
      //   std::cout << "index = " << inputIndex << std::endl << std::flush;
      //   utl::PrintUtlVector(x, "x");
      //   }

      if (this->m_Normalization==Superclass::UNIT_INTEGRAL)
        {
        x = this->NormalizeUnitIntegral(x);
        }

      b = (*this->m_BasisMatrix) * x;
      // if (this->GetDebug())
      //   utl::PrintUtlVector(b, "b");

      ScaleSamples(b);

      double maxVal = b.MaxValue();
      if (maxVal>bMax)
        bMax = maxVal;

      ++inputIt;
      }

    }
  
  inputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
  {
    inputIndex = inputIt.GetIndex();
    if (!this->IsPixelIndexVisible(inputIndex))
      {
      ++inputIt;
      progress.CompletedPixel();
      continue;
      }

    inputPtr->TransformIndexToPhysicalPoint(inputIndex, inputPhysicalPoint);
    inputPixel = inputIt.Get(); 

    for (unsigned int k=0;k<numberOfBasis;k++) //isotropic component?
      x(k) = inputPixel[k];

    if (x.GetRootMeanSquares() == 0)
      {
      ++inputIt;
      progress.CompletedPixel();
      continue;
      }

    if (this->GetDebug())
      {
      std::cout << "index = " << inputIndex << std::endl << std::flush;
      utl::PrintUtlVector(x, "x");
      }

    if (this->m_Normalization==Superclass::UNIT_INTEGRAL)
      {
      x = this->NormalizeUnitIntegral(x);
      }
    
    b = (*this->m_BasisMatrix) * x;
    if (this->GetDebug())
      utl::PrintUtlVector(b, "b");
      
    ScaleSamples(b);
    
    for (unsigned int k=0;k<numberOfPoints;k++) 
      {
      if ( this->m_ColorScheme == Superclass::MAGNITUDE )
        {
        colorTransferFunction->GetColor(b(k)/bMax, rgb);      
        for (unsigned int d=0;d<3;d++)
          {
          outputMeshPoint[d] = b(k) * (*this->m_Orientations)(k,d) + inputPhysicalPoint[d];
          RGB[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>(rgb[d]*255.0);
          }
        }
      else if ( this->m_ColorScheme == Superclass::DIRECTION )
        {
        for (unsigned int d=0;d<3;d++)
          {
          outputMeshPoint[d] = b(k) * (*this->m_Orientations)(k,d) + inputPhysicalPoint[d];
          RGB[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>(std::fabs( (*this->m_Orientations)(k,d))*255.0);
          }

        double min_val = std::numeric_limits<double>::max();
        double max_val = -std::numeric_limits<double>::max();

        for (unsigned int d=0;d<3;d++)
          {
          if (RGB[d] > max_val)
            max_val = RGB[d];
          if (RGB[d] < min_val)
            min_val = RGB[d];
          }

        // min_val=0;
        for (unsigned int d=0;d<3;d++)
          {
          RGB[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>( (RGB[d] - min_val)/(max_val - min_val) * 255.0 );
          }

        // double ss;
        // utl::RGBToIndex(std::fabs((*this->m_Orientations)(k,0)),std::fabs((*this->m_Orientations)(k,1)),std::fabs((*this->m_Orientations)(k,2)),ss);
        // newScalars->InsertNextTuple(&ss);
        // newScalars->SetName("rgb");
        }

      RGB[3] = 255.0;
      
      outputMeshPoints->InsertNextPoint( outputMeshPoint );
      outputMeshRGB->InsertNextTupleValue(RGB);
      }
    
    for (vtkIdType c=0;c<numberOfCells;c++) 
      {
      outputMeshCellArray->InsertNextCell(numberOfCellPoints);
      
      for( vtkIdType p = 0; p < numberOfCellPoints; p++ )
        {
        outputMeshCellArray->InsertCellPoint( offset + cellMatrix(c,p) );
        }
      }

    offset += numberOfPoints;
    
    ++inputIt;
    progress.CompletedPixel();

  }

    
  this->m_Mesh->SetPoints( outputMeshPoints );
  this->m_Mesh->SetPolys( outputMeshCellArray );
  this->m_Mesh->GetPointData()->SetScalars( outputMeshRGB );
  // int idx = this->m_Mesh->GetPointData()->SetScalars( newScalars );
  // this->m_Mesh->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
}


template<class TInputImage, class TOutputMesh>
void
MeshFromContinuousSphericalFunctionImageFilter<TInputImage, TOutputMesh>
::ScaleSamples(VectorType& b) const
{
  if (this->m_RemoveNegativeValues)
    {
    for (unsigned int k=0; k<b.Size(); k++ )
      {
      if ( b(k) < 0 )
        b(k) = 0;
      }
    }

  if (std::fabs(this->m_Pow-1.0)>1e-8)
    {
    for ( int k = 0; k < b.Size(); k += 1 ) 
      b(k) = std::pow(b(k), this->m_Pow);
    }


  if ( this->m_Normalization == Superclass::UNIT_MAX )
    {
    b = b/( b.MaxValue() + vnl_math::eps );
    }

  // do not perform min-max normalization for isotropic SHs
  if ( this->m_Normalization==Superclass::MIN_MAX && b.MaxValue()-b.MinValue()>1e-8)
    {
    b = (b - b.MinValue())/( b.MaxValue() - b.MinValue() + vnl_math::eps );
    }

  b %= this->m_Scale;
}

} // end namespace itk

#endif



