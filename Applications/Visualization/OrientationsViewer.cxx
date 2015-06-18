/**
 *       @file  OrientationsViewer.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-28-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkLODActor.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDelaunay3D.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkSphereSource.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkGlyph3D.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include "OrientationsViewerCLP.h"
#include "itkSamplingScheme3D.h"
#include "utl.h"

#include "utlVTKMacro.h"


/**
 * \brief  orientations visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main ( int argc, char *argv[] )
{
  
  PARSE_ARGS;

  int numberOfShells = _OrientationFile.size();
  utlGlobalException(numberOfShells==0, "no input orientation file");
  std::cout << numberOfShells << " shells" << std::endl << std::flush;

  if (!_OnlyCombined && _OpacityMesh.size()!=numberOfShells)
    {
    _OpacityMesh.resize(numberOfShells);
    for ( int i = 0; i < numberOfShells; i += 1 ) 
      _OpacityMesh[i]=1.0- 0.7/(double)(numberOfShells-1)*i;
    }
  
  if (!_OnlyCombined && _OpacitySphere.size()!=numberOfShells)
    {
    _OpacitySphere.resize(numberOfShells);
    for ( int i = 0; i < numberOfShells; i += 1 ) 
      _OpacitySphere[i]=_OpacityMesh[i];
    }

  std::vector<std::vector<double> > radiusVec;
  if (_RadiusFileArg.isSet())
    {
    double maxRadius = -1;
    utlGlobalException(_RadiusFile.size()!=numberOfShells, "wrong size of _RadiusFile");
    for ( int i = 0; i < numberOfShells; i += 1 ) 
      {
      std::vector<double> tmp;
      utl::ReadVector(_RadiusFile[i], tmp);
      radiusVec.push_back(tmp);
      double maxTmp = *std::max_element(tmp.begin(), tmp.end());
      if (maxTmp>maxRadius)
        maxRadius = maxTmp;
      }
    for ( int i = 0; i < numberOfShells; i += 1 ) 
      for ( int j = 0; j < radiusVec[i].size(); j += 1 ) 
        radiusVec[i][j] /= maxRadius;
    }
  else
    {
    if (_Radius.size()!=numberOfShells)
      {
      _Radius.resize(numberOfShells);
      for ( int i = 0; i < numberOfShells; i += 1 ) 
        _Radius[i]=i+1;
      }
    }

  if (_Color.size()!=3*numberOfShells)
    {
    _Color.resize(3*numberOfShells);
    if (numberOfShells==1)
      {
      _Color[0]=1.0, _Color[1]=1.0, _Color[2]=1.0;
      }
    else
      {
      if (numberOfShells>1)
        {
        _Color[0]=1.0, _Color[1]=0.0, _Color[2]=0.0;
        if (numberOfShells>=2)
          _Color[3]=0.0, _Color[4]=1.0, _Color[5]=0.0;
        if (numberOfShells==3)
          _Color[6]=0.0, _Color[7]=0.0, _Color[8]=1.0;
        }
      if (numberOfShells>3 && numberOfShells<=6)
        {
        _Color[9]=0.5, _Color[10]=0.0, _Color[11]=0.0;
        if (numberOfShells>=7)
          _Color[12]=0.0, _Color[13]=0.5, _Color[14]=0.0;
        if (numberOfShells==8)
          _Color[15]=0.0, _Color[16]=0.0, _Color[17]=0.5;
        }
      }
    utlGlobalException(numberOfShells>7, "too many shells, please manually set _Color");
    }
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderer->SetBackground(_BackgroundColor[0], _BackgroundColor[1], _BackgroundColor[2]);
    
  typedef itk::SamplingScheme3D<double> SamplingType;

  std::vector<vtkSmartPointer<vtkPoints> > points(numberOfShells+1); 
  std::vector<vtkSmartPointer<vtkPolyData> > polydata(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkGlyph3D> > glyph3D(numberOfShells+1);

  std::vector<vtkSmartPointer<vtkVertexGlyphFilter> > vertexGlyphFilter(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkPolyDataMapper> > pointsMapper(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkLODActor> > pointsActor(numberOfShells+1);
  
  std::vector<vtkSmartPointer<vtkDelaunay3D> > delaunay3D(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkDataSetMapper> > delaunayMapper(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkLODActor> > delaunayActor(numberOfShells+1);
  
  std::vector<vtkSmartPointer<vtkSphereSource> > sphereSource(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkPolyDataMapper> > sphereMapper(numberOfShells+1);
  std::vector<vtkSmartPointer<vtkLODActor> > sphereActor(numberOfShells+1);

  unsigned char colorPoint[3] = {0, 0, 0};
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetName("colors");

  if (_OnlyCombined)
    {
    points.back() = vtkSmartPointer<vtkPoints>::New();
    polydata.back() = vtkSmartPointer<vtkPolyData>::New();
    colors->SetNumberOfComponents(3);
    }

  vtkSmartPointer<vtkSphereSource> sphereSourceTemp = vtkSmartPointer<vtkSphereSource>::New();
  sphereSourceTemp->SetCenter(0.0, 0.0, 0.0);
  sphereSourceTemp->SetRadius(_SpherePointSize);
  sphereSourceTemp->SetThetaResolution(20);
  sphereSourceTemp->SetPhiResolution(20);


  for ( int i = 0; i < _OrientationFile.size(); i += 1 ) 
    {
    SamplingType::Pointer sampling = SamplingType::New();

    if (_NoSymmetricPoints)
      sampling->ReadOrientationFile(_OrientationFile[i], DIRECTION_NODUPLICATE);
    else
      sampling->ReadOrientationFile(_OrientationFile[i], DIRECTION_DUPLICATE);
    // utlPrintVar2(true, _NoSymmetricPoints, sampling->GetNumberOfSamples());
    std::cout << "shell " << i+1 << ", " << sampling->GetNumberOfSamples() << (_NoSymmetricPoints?" points":" antipodal symmetric points") << std::endl << std::flush;

    points[i] = vtkSmartPointer<vtkPoints>::New();

    if (_OnlyCombined)
      {
      colorPoint[0] = (unsigned char)(255*_Color[3*i]);
      colorPoint[1] = (unsigned char)(255*_Color[3*i+1]);
      colorPoint[2] = (unsigned char)(255*_Color[3*i+2]);
      }
    // utlPrintVar3(true, colorPoint[0], colorPoint[1], colorPoint[2]);
    // utlPrintVar3(true, 255*_Color[3*i], 255*_Color[3*i+1], 255*_Color[3*i+2]);

    for ( unsigned int j = 0; j < sampling->GetNumberOfSamples(); j++ )
      {
      double x0 = (*sampling)[j][0];
      double x1 = (*sampling)[j][1];
      double x2 = (*sampling)[j][2];
      // utlPrintVar3(true, x0,x1,x2);
      if (_RadiusFileArg.isSet())
        points[i]->InsertNextPoint( radiusVec[i][j]*x0, radiusVec[i][j]*x1, radiusVec[i][j]*x2 );
      else
        points[i]->InsertNextPoint( _Radius[i]*x0, _Radius[i]*x1, _Radius[i]*x2 );
      if (_OnlyCombined)
        {
        points.back()->InsertNextPoint( x0, x1, x2 );
        colors->InsertNextTupleValue(colorPoint);
        }
      }

    polydata[i] = vtkSmartPointer<vtkPolyData>::New();
    polydata[i]->SetPoints(points[i]);
    // polydata[i]->Print(std::cout<<"polydata[i]\n");

    if (!_NoPoints && !_OnlyCombined)
      {
      if (_PointType=="SPHERE")
        {
        glyph3D[i] = vtkSmartPointer<vtkGlyph3D>::New();
        glyph3D[i]->SetSourceConnection(sphereSourceTemp->GetOutputPort());
        vtkSetInputData(glyph3D[i], polydata[i]);
        glyph3D[i]->Update();
        pointsMapper[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSetInputData(pointsMapper[i], glyph3D[i]->GetOutput());
        }
      else if (_PointType=="VERTEX")
        {
        vertexGlyphFilter[i] = vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vtkSetInputData(vertexGlyphFilter[i], polydata[i]);
        vertexGlyphFilter[i]->Update();
        pointsMapper[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSetInputData(pointsMapper[i], vertexGlyphFilter[i]->GetOutput());
        }


      // unsigned char red[3] = {255, 0, 0};
      // vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
      // colors->SetNumberOfComponents(3);
      // for ( unsigned int i = 0; i < sampling->GetNumberOfSamples(); i++ )
      //   {
      //   colors->InsertNextTupleValue(red);
      //   }
      // polydata[i]->GetPointData()->SetScalars(colors);


      pointsActor[i] = vtkSmartPointer<vtkLODActor>::New();
      pointsActor[i]->SetMapper(pointsMapper[i]);
      pointsActor[i]->GetProperty()->SetColor(_Color[3*i],_Color[3*i+1],_Color[3*i+2]);
      if (_PointType=="VERTEX")
        pointsActor[i]->GetProperty()->SetPointSize(_VertexPointSize);

      renderer->AddActor(pointsActor[i]);
      }

    if (_Mesh && !_OnlyCombined)
      {
      delaunay3D[i] = vtkSmartPointer<vtkDelaunay3D>::New();
      vtkSetInputData(delaunay3D[i], polydata[i]);
      delaunay3D[i]->SetTolerance(1e-8);
      delaunay3D[i]->SetOffset(10);
      delaunay3D[i]->Update();
      // delaunay3D[i]->Print(std::cout<<"delaunay3D[i]=\n");

      delaunayMapper[i] = vtkSmartPointer<vtkDataSetMapper>::New();
      delaunayMapper[i]->SetInputConnection(delaunay3D[i]->GetOutputPort());

      delaunayActor[i] = vtkSmartPointer<vtkLODActor>::New();
      delaunayActor[i]->SetMapper(delaunayMapper[i]);
      delaunayActor[i]->GetProperty()->SetColor(1,1,1);
      delaunayActor[i]->GetProperty()->SetOpacity(_OpacityMesh[i]);

      renderer->AddActor(delaunayActor[i]);
      }

    if (_Sphere && !_OnlyCombined)
      {
      sphereSource[i] = vtkSmartPointer<vtkSphereSource>::New();
      sphereSource[i]->SetCenter(0.0, 0.0, 0.0);
      utlGlobalException(_RadiusFileArg.isSet(), "not support to visualize sphere if _RadiusFile is set");
      sphereSource[i]->SetRadius(_Radius[i]);
      sphereSource[i]->SetThetaResolution(50);
      sphereSource[i]->SetPhiResolution(50);

      sphereMapper[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
      sphereMapper[i]->SetInputConnection(sphereSource[i]->GetOutputPort());

      sphereActor[i] = vtkSmartPointer<vtkLODActor>::New();
      sphereActor[i]->SetMapper(sphereMapper[i]);
      // sphereActor[i]->GetProperty()->SetColor(_Color[3*i],_Color[3*i+1],_Color[3*i+2]);
      sphereActor[i]->GetProperty()->SetColor(1,1,1);
      sphereActor[i]->GetProperty()->SetOpacity(_OpacitySphere[i]);

      renderer->AddActor(sphereActor[i]);
      }
    }

  if (_OnlyCombined)
    {
    polydata.back() = vtkSmartPointer<vtkPolyData>::New();
    polydata.back()->SetPoints(points.back());
    polydata.back()->GetPointData()->AddArray(colors);

    if (!_NoPoints)
      {
      vtkSmartPointer<vtkPolyData> polydata2 = vtkSmartPointer<vtkPolyData>::New();
      if (_PointType=="SPHERE")
        {
        glyph3D.back() = vtkSmartPointer<vtkGlyph3D>::New();
        glyph3D.back()->SetSourceConnection(sphereSourceTemp->GetOutputPort());
        vtkSetInputData(glyph3D.back(), polydata.back());
        glyph3D.back()->Update();
        // polydata2->ShallowCopy(glyph3D.back()->GetOutput());
//        glyph3D.back()->GetOutput()->GetPointData()->SetScalars(colors);
        pointsMapper.back() = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSetInputData(pointsMapper.back(), glyph3D.back()->GetOutput());
        }
      else if (_PointType=="VERTEX")
        {
        vertexGlyphFilter.back() = vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vtkSetInputData(vertexGlyphFilter.back(), polydata.back());
        vertexGlyphFilter.back()->Update();
        // polydata2->ShallowCopy(vertexGlyphFilter.back()->GetOutput());
//        vertexGlyphFilter.back()->GetOutput()->GetPointData()->SetScalars(colors);
        pointsMapper.back() = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSetInputData(pointsMapper.back(), vertexGlyphFilter.back()->GetOutput());
        }

      pointsMapper.back()->ScalarVisibilityOn();
      pointsMapper.back()->SetScalarModeToUsePointFieldData();
      pointsMapper.back()->SetColorModeToDefault();
      pointsMapper.back()->SelectColorArray("colors");

      // polydata2->GetPointData()->SetScalars(colors);

      // pointsMapper.back() = vtkSmartPointer<vtkPolyDataMapper>::New();
      // vtkSetInputData(pointsMapper.back(), polydata2);

      pointsActor.back() = vtkSmartPointer<vtkLODActor>::New();
      pointsActor.back()->SetMapper(pointsMapper.back());
      if (_PointType=="VERTEX")
        pointsActor.back()->GetProperty()->SetPointSize(_VertexPointSize);

      renderer->AddActor(pointsActor.back());
      }

    if (_Mesh)
      {
      delaunay3D.back() = vtkSmartPointer<vtkDelaunay3D>::New();
      vtkSetInputData(delaunay3D.back(), polydata.back());
      delaunay3D.back()->SetTolerance(1e-8);
      delaunay3D.back()->SetOffset(10);
      delaunay3D.back()->Update();
      // delaunay3D.back()->Print(std::cout<<"delaunay3D.back()=\n");

      delaunayMapper.back() = vtkSmartPointer<vtkDataSetMapper>::New();
      delaunayMapper.back()->SetInputConnection(delaunay3D.back()->GetOutputPort());

      delaunayActor.back() = vtkSmartPointer<vtkLODActor>::New();
      delaunayActor.back()->SetMapper(delaunayMapper.back());
      delaunayActor.back()->GetProperty()->SetColor(1,1,1);
      delaunayActor.back()->GetProperty()->SetOpacity(_OpacityMesh[0]);

      renderer->AddActor(delaunayActor.back());
      }

    if (_Sphere)
      {
      sphereSource.back() = vtkSmartPointer<vtkSphereSource>::New();
      sphereSource.back()->SetCenter(0.0, 0.0, 0.0);
      sphereSource.back()->SetRadius(1.0);
      sphereSource.back()->SetThetaResolution(50);
      sphereSource.back()->SetPhiResolution(50);

      sphereMapper.back() = vtkSmartPointer<vtkPolyDataMapper>::New();
      sphereMapper.back()->SetInputConnection(sphereSource.back()->GetOutputPort());

      sphereActor.back() = vtkSmartPointer<vtkLODActor>::New();
      sphereActor.back()->SetMapper(sphereMapper.back());
      // sphereActor.back()->GetProperty()->SetColor(_Color[3*i],_Color[3*i+1],_Color[3*i+2]);
//      sphereActor.back()->GetProperty()->SetColor(1,1,1);
//      sphereActor.back()->GetProperty()->SetOpacity(_OpacitySphere[0]);

      renderer->AddActor(sphereActor.back());
      }
    }
  
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(_Window[0],_Window[1]);
  
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =  vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->SetDesiredUpdateRate(25);

  renderWindow->Render();

  if (_PNGFile!="")
    {
    vtkSmartPointer<vtkWindowToImageFilter> windowToImage = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImage->SetInput(renderWindow);
    windowToImage->SetMagnification(2);
    vtkSmartPointer<vtkPNGWriter> pngWriter = vtkSmartPointer<vtkPNGWriter>::New();
    pngWriter->SetInputConnection(windowToImage->GetOutputPort());
    pngWriter->SetFileName(_PNGFile.c_str());
    pngWriter->Write();
    }
  else
    renderWindowInteractor->Start();



  return EXIT_SUCCESS;
}

