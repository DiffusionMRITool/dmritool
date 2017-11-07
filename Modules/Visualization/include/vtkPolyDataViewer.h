/**
 *       @file  vtkPolyDataViewer.h
 *      @brief  
 *     Created  "06-08-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __vtkPolyDataViewer_h
#define __vtkPolyDataViewer_h

#include <vtkObject.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkNew.h>
#include <vtkProperty.h>

#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkLODActor.h>
#include <vtkPolyDataNormals.h>
#include <vtkObjectFactory.h>


#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include "utlCoreMacro.h"

namespace vtk
{

/** helper function from vtkPolyData to vtkLODActor  */
inline vtkSmartPointer<vtkLODActor>
vtkPolyDataToActor(const vtkSmartPointer<vtkPolyData>& mesh, const double opacity=1.0, const std::vector<double>& hueRange={0.6667,0.0}, bool useNormal=true, const std::vector<double>& scalarRange={-1.0,-1.0}, bool lighting=true )
  {
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

  if (useNormal)
    {
    vtkSmartPointer<vtkPolyDataNormals> polyDataNormals = vtkPolyDataNormals::New();
    polyDataNormals->SetInputData(mesh);
    mapper->SetInputConnection(polyDataNormals->GetOutputPort());
    }
  else
    mapper->SetInputData(mesh);

  // lut when scalar is used for colors
  if (mesh->GetPointData()->GetScalars() && mesh->GetPointData()->GetScalars()->GetNumberOfComponents()==1)
    {
    vtkSmartPointer<vtkLookupTable> lut = vtkLookupTable::New();
    double* valueRange = mesh->GetScalarRange();
    double rangeUsed[2];
    rangeUsed[0] = std::fabs(scalarRange[0]+1.0) < 1.0e-8 ? valueRange[0] : scalarRange[0]; 
    rangeUsed[1] = std::fabs(scalarRange[1]+1.0) < 1.0e-8 ? valueRange[1] : scalarRange[1]; 
    // utlPrintVar(true, rangeUsed[0], rangeUsed[1]);
    lut->SetTableRange(rangeUsed[0], rangeUsed[1]);
    lut->SetHueRange(hueRange[0],hueRange[1]);
    lut->SetRampToLinear();
    lut->Build();

    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(rangeUsed[0], rangeUsed[1]);
    }

  vtkSmartPointer<vtkLODActor> actor = vtkSmartPointer<vtkLODActor>::New();
  actor->SetMapper(mapper);

  actor->GetProperty()->SetLighting(lighting);
  actor->GetProperty()->SetOpacity(opacity);

  return actor;
  }

/** \class vtkPolyDataViewer
 * \brief  helper class to visualize vtkPolyData.
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
class vtkPolyDataViewer : public vtkObject
  {
public:
  static vtkPolyDataViewer *New();
  vtkTypeMacro(vtkPolyDataViewer, vtkObject);

  vtkSetVector2Macro(ScalarRange, double);

  vtkSetVector2Macro(HueRange, double);

  vtkSetMacro (UseNormal, bool);
  vtkGetMacro (UseNormal, bool);
  
  vtkSetMacro (Zoom, double);
  vtkGetMacro (Zoom, double);
  
  vtkSetMacro (Lighting, bool);
  vtkGetMacro (Lighting, bool);

  vtkSetVector2Macro(Angle, double);
  vtkGetVector2Macro(Angle, double);
  
  vtkSetVector2Macro(WindowSize, int);
  vtkGetVector2Macro(WindowSize, int);

  /** add mesh  */
  void Add(const vtkSmartPointer<vtkPolyData>& mesh, const double opacity=1.0);

  /** add mesh with hueRange, useNormal */
  void Add(const vtkSmartPointer<vtkPolyData>& mesh, const double opacity, const std::vector<double>& hueRange, bool useNormal, const std::vector<double>& scalarRange={-1.0,-1.0}, bool lighting=true );
  
  void Add(const vtkSmartPointer<vtkLODActor>& actor)
    {
    Renderer->AddActor(actor);
    }

  void SetBackground(double r, double g, double b)
    {
    Renderer->SetBackground(r,g,b);
    }
  
  void Show();
  
  void SavePNG(const std::string& pngfile);

protected: 

  vtkPolyDataViewer(){}
  ~vtkPolyDataViewer(){}

  
  void SetCamara();

  void RenderWindowUpdate();

private:

  vtkSmartPointer<vtkRenderer> Renderer = vtkSmartPointer<vtkRenderer>::New();
  
  vtkSmartPointer<vtkRenderWindow> RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> Interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  /** used to normalize the scalars for coloring  */
  double ScalarRange[2]={-1.0,-1.0};

  bool UseNormal=true;

  bool Lighting=true;

  int WindowSize[2]={600,600};

  double Angle[2]={0.0,0.0};

  double HueRange[2]={0.6667,0.0};

  double Zoom=1.0;

  };

vtkStandardNewMacro(vtkPolyDataViewer);

/** helper function to visualize vtkPolyData  */
inline void 
VisualizePolyDataWithScalarRange (const std::vector<vtkPolyData*>& meshes, const std::vector<double>& opacity, const std::vector<double>& scalarRange={-1.0,-1.0}, const std::vector<double>& hueRange={0.6667,0.0}, const std::vector<double>& angle={0.0,0.0}, const std::vector<int>& windowSize={600,600}, const bool useNormal=true, const bool lighting=true, const double zoom=1.0, const std::string& pngfile="", const std::vector<double>& bgColor={0,0,0} )
{
  vtkSmartPointer<vtkPolyDataViewer> viewer = vtkPolyDataViewer::New();

  viewer->SetScalarRange(scalarRange[0], scalarRange[1]);
  viewer->SetHueRange(hueRange[0], hueRange[1]);
  viewer->SetAngle(angle[0], angle[1]);
  viewer->SetWindowSize(windowSize[0], windowSize[1]);
  viewer->SetUseNormal(useNormal);
  viewer->SetLighting(lighting);
  viewer->SetZoom(zoom);
  viewer->SetBackground(bgColor[0], bgColor[1], bgColor[2]);

  utlException(meshes.size()!=opacity.size(), "meshes and opacity should have the same size");

  for ( int i = 0; i < meshes.size(); ++i ) 
    {
    if (meshes[i])
      viewer->Add(meshes[i], opacity[i]);
    }

  if (pngfile=="")
    viewer->Show();
  else
    viewer->SavePNG(pngfile);
}

/** helper function to visualize vtkPolyData  */
inline void 
VisualizePolyDataWithScalarRange ( vtkPolyData* mesh, const std::vector<double>& scalarRange={-1.0,-1.0}, const std::vector<double>& hueRange={0.6667,0.0}, const std::vector<double>& angle={0.0,0.0}, const std::vector<int>& windowSize={600,600}, const bool useNormal=true, const bool lighting=true, const double zoom=1.0, const std::string& pngfile="", const std::vector<double>& bgColor={0,0,0} )
{
  vtkSmartPointer<vtkPolyDataViewer> viewer = vtkPolyDataViewer::New();

  viewer->SetScalarRange(scalarRange[0], scalarRange[1]);
  viewer->SetHueRange(hueRange[0], hueRange[1]);
  viewer->SetAngle(angle[0], angle[1]);
  viewer->SetWindowSize(windowSize[0], windowSize[1]);
  viewer->SetUseNormal(useNormal);
  viewer->SetLighting(lighting);
  viewer->SetZoom(zoom);
  viewer->SetBackground(bgColor[0], bgColor[1], bgColor[2]);

  viewer->Add(mesh);

  if (pngfile=="")
    viewer->Show();
  else
    viewer->SavePNG(pngfile);
}

/** helper function to visualize vtkPolyData  */
inline void 
VisualizePolyData ( vtkPolyData* mesh, const std::vector<double>& angle={0.0,0.0}, const std::vector<int>& windowSize={600,600}, const bool useNormal=true, const bool lighting=true, const double zoom=1.0, const std::string& pngfile="", const std::vector<double>& bgColor={0,0,0} )
{
  VisualizePolyDataWithScalarRange(mesh, {-1.0,-1.0}, {0.6667,0.0}, angle, windowSize, useNormal, lighting, zoom, pngfile, bgColor);
}

}


#if !defined(__vtkPolyDataViewer_hxx)
#include "vtkPolyDataViewer.hxx"
#endif

#endif 
