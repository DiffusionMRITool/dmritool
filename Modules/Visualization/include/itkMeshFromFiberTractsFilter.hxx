/**
 *       @file  itkMeshFromFiberTractsFilter.hxx
 *      @brief  
 *     Created  "08-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromFiberTractsFilter_hxx
#define __itkMeshFromFiberTractsFilter_hxx

#include "itkMeshFromFiberTractsFilter.h"

#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTubeFilter.h>
#include <vtkCellArray.h>

#include "utlNDArray.h"
#include "utlCore.h"
#include "utlDMRI.h"

namespace itk
{

typename LightObject::Pointer
MeshFromFiberTractsFilter
::InternalClone() const
{
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }

  rval->m_FiberTracts = m_FiberTracts;
  rval->m_ColorScheme = m_ColorScheme;
  rval->m_ShapeMode = m_ShapeMode;
  rval->m_TubeRadius = m_TubeRadius;
  rval->m_Mesh = m_Mesh;
  rval->m_Flip = m_Flip;

  for ( int i = 0; i < 3; ++i ) 
    rval->m_Color[i] = m_Color[i];

  return loPtr;
}

void
MeshFromFiberTractsFilter
::Update()
{
  vtkSmartPointer<vtkPoints>    points    = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  // vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  scalars->SetNumberOfComponents(3);
  scalars->SetName("rgb_color");
  double rgb[3];


  int numFibers = m_FiberTracts->GetNumberOfFibers();

  utl::Vector<double> dirMean(3), dir(3);
  
  long count=0;
  double color = 0;
  for ( int n = 0; n < numFibers; ++n ) 
    {
    auto fiber = m_FiberTracts->GetFiber(n);
    int numPoints = fiber->GetNumberOfPoints();

    utlPrintVar(this->GetDebug(), n, numPoints);


    if (numPoints>0)
      {

      lines->InsertNextCell(numPoints);

      if (m_ColorScheme==COLOR_BY_MEAN_DIRECTION)
        dirMean.Fill(0.0);

      for ( int i = 0; i < numPoints; ++i ) 
        {
        auto vertex = fiber->GetPoint(i);
        utl::FlipVector(vertex, m_Flip, 3);
        points->InsertPoint(count, vertex[0], vertex[1], vertex[2]);
        lines->InsertCellPoint(count);

        if (m_ColorScheme==COLOR_FIXED)
          {
          scalars->InsertTuple(count, m_Color);
          }
        else if (m_ColorScheme==COLOR_BY_POINT_DIRECTION)
          {
          auto diff = fiber->GetDirection(i);
          utl::FlipVector(diff, m_Flip, 3);
          utl::VectorToVector(diff, dir, 3);
          double dirNorm = dir.GetTwoNorm();
          if (dirNorm>0)
            dir /= dirNorm;

          for ( int d = 0; d < 3; ++d ) 
            rgb[d] = std::fabs(dir[d])*255.0;
          scalars->InsertTuple(count, rgb);
          }
        else if (m_ColorScheme==COLOR_BY_MEAN_DIRECTION)
          {
          auto diff = fiber->GetDirection(i);
          utl::FlipVector(diff, m_Flip, 3);
          utl::VectorToVector(diff, dir, 3);
          double dirNorm = dir.GetTwoNorm();
          if (dirNorm>0)
            dir /= dirNorm;
          dirMean += dir;
          }
        else if (m_ColorScheme==COLOR_BY_IMAGE)
          {
          utlGlobalException(true, "TODO");

          }
        else if (m_ColorScheme==COLOR_BY_SCALARS)
          {
          utlGlobalException(true, "TODO");

          }
        else if (m_ColorScheme==COLOR_BY_PROPERTY)
          {
          utlGlobalException(true, "TODO");

          }

        count++;
        }

      if (m_ColorScheme==COLOR_BY_MEAN_DIRECTION)
        {
        double dirNorm = dirMean.GetTwoNorm();
        if (dirNorm>0)
          dirMean /= dirNorm;
        for ( int d = 0; d < 3; ++d ) 
          rgb[d] = std::fabs(dirMean[d])*255.0;
        for ( int i = 0; i < numPoints; ++i ) 
          scalars->InsertTuple(count-i, rgb);
        }
      else if (m_ColorScheme==COLOR_BY_ENDPOINTS_DIRECTION)
        {
        VertexType p0 = fiber->GetPoint(0);
        VertexType p1 = fiber->GetPoint(numPoints-1);
        auto diff = p1-p0;
        utl::VectorToVector(diff, dir, 3);
        double dirNorm = dir.GetTwoNorm();
        if (dirNorm>0)
          dir /= dirNorm;
        for ( int d = 0; d < 3; ++d ) 
          rgb[d] = std::fabs(dir[d])*255.0;
        for ( int i = 0; i < numPoints; ++i ) 
          scalars->InsertTuple(count-i, rgb);
        }

      }
    }


  m_Mesh->SetPoints( points );
  m_Mesh->SetLines( lines );

  if (m_ColorScheme!=COLOR_NONE)
    (m_Mesh->GetPointData())->SetScalars(scalars);

  m_Mesh->Squeeze();

  if (m_TubeRadius>0 && m_ShapeMode==GLYPH_TUBE)
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


