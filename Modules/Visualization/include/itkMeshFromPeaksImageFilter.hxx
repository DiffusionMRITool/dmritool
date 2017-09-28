/**
 *       @file  itkMeshFromPeaksImageFilter.hxx
 *      @brief  
 *     Created  "08-26-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkMeshFromPeaksImageFilter_hxx
#define __itkMeshFromPeaksImageFilter_hxx

#include "itkMeshFromPeaksImageFilter.h"
#include <vtkTubeFilter.h>

namespace itk
{
template <class TInputImage, class TOutputMesh>
MeshFromPeaksImageFilter<TInputImage, TOutputMesh>
::MeshFromPeaksImageFilter() : Superclass()
{
  m_PeakType = XYZV;
  m_TubeRadius = 0.1*this->m_Scale;
  m_MaxNumberOfPeaks = -1;
    
  this->m_ColorScheme = Superclass::FIXED;
  m_ColorPeak[0]=255, m_ColorPeak[1]=0, m_ColorPeak[2]=0;
}

template <class TInputImage, class TOutputMesh>
typename LightObject::Pointer
MeshFromPeaksImageFilter<TInputImage, TOutputMesh>
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_PeakType = m_PeakType;
  rval->m_TubeRadius = m_TubeRadius;
  rval->m_MaxNumberOfPeaks = m_MaxNumberOfPeaks;

  for ( int i = 0; i < 3; ++i ) 
    rval->m_ColorPeak[i] = m_ColorPeak[i];

  return loPtr;
}


template <class TInputImage, class TOutputMesh>
void MeshFromPeaksImageFilter<TInputImage, TOutputMesh>
::GenerateData()
{
  VerifyInputParameters();

  // Pointers
  InputImageConstPointer inputPtr = this->GetInput();

  // Input Image Parameters
  InputImageRegionType inputRegion = inputPtr->GetRequestedRegion();
  InputImageIndexType inputIndex;
  InputImagePointType inputPhysicalPoint;
  InputImagePixelType inputPixel;
  
  // iterator for the input image
  ImageRegionConstIterator<InputImageType> inputIt(inputPtr, inputPtr->GetRequestedRegion());
    
  vtkSmartPointer<OutputMeshPointsType> outputMeshPoints = vtkSmartPointer<OutputMeshPointsType>::New();
  vtkSmartPointer<OutputMeshCellArrayType> outputMeshCellArray = vtkSmartPointer<OutputMeshCellArrayType>::New();
  vtkSmartPointer<OutputMeshRGBType> outputMeshRGB = vtkSmartPointer<OutputMeshRGBType>::New();
  outputMeshRGB->SetNumberOfComponents(3);
  outputMeshRGB->SetName("RGB_scalars");
  
  double outputMeshPoint1[3];
  double outputMeshPoint2[3];
  unsigned int numberOfCellPoints = 2;

  unsigned long offset=0;
  unsigned char RGB[3];  
  unsigned int numberOfComponentsPerOrientation = PeakContainerHelper::GetDimensionPerPeak(m_PeakType);
  int numberOfPeaks = PeakContainerHelper::GetNumberOfPeaks(m_PeakType, inputPtr->GetNumberOfComponentsPerPixel());
  int realNumberOfPeaks = m_MaxNumberOfPeaks>0 ? utl::min(numberOfPeaks, m_MaxNumberOfPeaks) : numberOfPeaks;
  double norm = 0;

  std::vector<double> peak(3);
  
  while( !inputIt.IsAtEnd() )
  {
    inputIndex = inputIt.GetIndex();
    inputPtr->TransformIndexToPhysicalPoint(inputIndex, inputPhysicalPoint);
    inputPixel = inputIt.Get();
    
    for (unsigned int k=0;k<realNumberOfPeaks;k++) //isotropic component?
      {
      norm = 0;
      peak = PeakContainerHelper::GetPeak(inputPixel, k, m_PeakType);
      for ( int d = 0; d < 3; ++d ) 
        peak[d] = this->m_Flip[d]? -peak[d] : peak[d];
      
      for (unsigned int d=0; d<3; d++)
        {
        norm += peak[d]*peak[d];
        
        outputMeshPoint1[d] = peak[d] * this->m_Scale + inputPhysicalPoint[d];
        outputMeshPoint2[d] = -peak[d] * this->m_Scale + inputPhysicalPoint[d];
        if (this->m_ColorScheme==Superclass::DIRECTION)
          RGB[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>(std::fabs(peak[d])*255.0);
        else
          RGB[d] = static_cast<VTK_TYPE_NAME_UNSIGNED_CHAR>(m_ColorPeak[d]);
        }

      if ( norm > 0 )
        {
        outputMeshPoints->InsertNextPoint( outputMeshPoint1 );
        outputMeshPoints->InsertNextPoint( outputMeshPoint2 );
        outputMeshRGB->InsertNextTupleValue(RGB);
        outputMeshRGB->InsertNextTupleValue(RGB);

        outputMeshCellArray->InsertNextCell(numberOfCellPoints);
        
        for( vtkIdType p = 0; p < numberOfCellPoints; p++ )
          outputMeshCellArray->InsertCellPoint( offset + p );

        offset += numberOfCellPoints;
        }

      }

    ++inputIt;
  }
  
  this->m_Mesh->SetPoints( outputMeshPoints );
  this->m_Mesh->SetLines( outputMeshCellArray );
  this->m_Mesh->GetPointData()->SetScalars( outputMeshRGB );

  if (m_TubeRadius>0)
    {
    vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkTubeFilter::New();
    tubeFilter->SetInputData(this->m_Mesh);
    tubeFilter->SetRadius(m_TubeRadius);
    tubeFilter->SetNumberOfSides(6);
    tubeFilter->Update();
    this->m_Mesh = tubeFilter->GetOutput();
    }
}


}

#endif 

