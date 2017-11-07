/**
 *       @file  vtkPolyDataViewer.hxx
 *      @brief  
 *     Created  "06-08-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __vtkPolyDataViewer_hxx
#define __vtkPolyDataViewer_hxx

#include "vtkPolyDataViewer.h"

namespace vtk
{

void vtkPolyDataViewer::Add(const vtkSmartPointer<vtkPolyData>& mesh, const double opacity)
{
  vtkSmartPointer<vtkLODActor> actor = vtkPolyDataToActor(mesh, opacity, {HueRange[0], HueRange[1]}, UseNormal, {ScalarRange[0],ScalarRange[1]}, Lighting);
  this->Add(actor);
}

void vtkPolyDataViewer::Add(const vtkSmartPointer<vtkPolyData>& mesh, const double opacity, const std::vector<double>& hueRange, bool useNormal, const std::vector<double>& scalarRange, bool lighting )
{
  vtkSmartPointer<vtkLODActor> actor;
  if (std::fabs(scalarRange[0]+1.0)<1e-8 && std::fabs(scalarRange[1]+1.0)<1e-8)
    actor = vtkPolyDataToActor(mesh, opacity, {hueRange[0], hueRange[1]}, useNormal, {ScalarRange[0],ScalarRange[1]}, lighting);
  else
    actor = vtkPolyDataToActor(mesh, opacity, {hueRange[0], hueRange[1]}, useNormal, {scalarRange[0],scalarRange[1]}, lighting);
  this->Add(actor);
}

void vtkPolyDataViewer::SetCamara()
{
  vtkSmartPointer<vtkCamera> camera = Renderer->GetActiveCamera();
  camera->Roll(Angle[0]);
  camera->Elevation(Angle[1]);
  camera->Zoom(Zoom);
}

void vtkPolyDataViewer::RenderWindowUpdate()
{
  RenderWindow->AddRenderer( Renderer );
  RenderWindow->SetSize(WindowSize[0], WindowSize[1]);

  Interactor->SetRenderWindow( RenderWindow );
  Interactor->Initialize();
  
  RenderWindow->Render();

  Renderer->ResetCamera();
  this->SetCamara();

  RenderWindow->Render();
}

void vtkPolyDataViewer::Show()
{
  this->RenderWindowUpdate();

  Interactor->Start();
}

void vtkPolyDataViewer::SavePNG(const std::string& pngfile)
{
  this->RenderWindowUpdate();

  vtkSmartPointer<vtkWindowToImageFilter> window2image= vtkWindowToImageFilter::New();
  window2image->SetInput(RenderWindow);
  vtkSmartPointer<vtkPNGWriter> png_writer = vtkPNGWriter::New();
  png_writer->SetInputConnection(window2image->GetOutputPort());
  png_writer->SetFileName(pngfile.c_str());
  png_writer->Write();
}

}


#endif 

