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

void vtkPolyDataViewer::Add(const vtkSmartPointer<vtkPolyData>& mesh)
{
  vtkSmartPointer<vtkLODActor> actor = vtkPolyDataToActor(mesh, {HueRange[0], HueRange[1]}, UseNormal, {ScalarRange[0],ScalarRange[1]});
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

