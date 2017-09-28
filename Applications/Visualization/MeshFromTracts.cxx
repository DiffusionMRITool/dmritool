/**
 *       @file  MeshFromTracts.cxx
 *      @brief  
 *     Created  "08-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "MeshFromTractsCLP.h"

#include "vtkPolyDataWriter.h"

#include "itkCommandProgressUpdate.h"

#include "utl.h"
#include "utlVTK.h"

#include "itkFiberTractsReader.h"
#include "itkMeshFromFiberTractsFilter.h"
#include "vtkPolyDataViewer.h"

/**
 * \brief  Gnerate mesh from peaks
 */
  int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;
  
  itk::TimeProbe clock;
  
  utlGlobalException(_Flip.size()!=3, "need 3 parameters in --flip");

  itk::FiberTractsReader::Pointer reader = itk::FiberTractsReader::New();
  reader->SetFileName(_InputFile);
  reader->Update();

  auto fibers = reader->GetFiberTracts();

  // fibers->Print(std::cout<<"fibers=\n");

  typedef itk::MeshFromFiberTractsFilter MeshCreatorType;
  MeshCreatorType::Pointer filter = MeshCreatorType::New();

  if (_ColorScheme == "DIRECTION") { filter->SetColorScheme(MeshCreatorType::COLOR_BY_POINT_DIRECTION); }
  if (_ColorScheme == "MEAN_DIRECTION") { filter->SetColorScheme(MeshCreatorType::COLOR_BY_MEAN_DIRECTION); }
  if (_ColorScheme == "END_POINTS") { filter->SetColorScheme(MeshCreatorType::COLOR_BY_ENDPOINTS_DIRECTION); }
  if (_ColorScheme == "FIXED") { filter->SetColorScheme(MeshCreatorType::COLOR_FIXED); }
  if (_ColorScheme == "IMAGE") { filter->SetColorScheme(MeshCreatorType::COLOR_BY_IMAGE); }
  if (_ColorScheme == "SCALARS") { filter->SetColorScheme(MeshCreatorType::COLOR_BY_SCALARS); }
  if (_ColorScheme == "PROPERTY") { filter->SetColorScheme(MeshCreatorType::COLOR_BY_PROPERTY); }

  filter->SetFiberTracts( fibers );
  filter->SetColor(_ColorFiber[0], _ColorFiber[1], _ColorFiber[2]);
  filter->SetFlip(_Flip[0], _Flip[1], _Flip[2]);
  if (_TubeRadiusArg.isSet())
    {
    filter->SetShapeMode(MeshCreatorType::GLYPH_TUBE);
    filter->SetTubeRadius(_TubeRadius);
    }
  else
    filter->SetShapeMode(MeshCreatorType::GLYPH_LINE);

  filter->SetDebug(_DebugArg.isSet());
  
  std::cout << "Generating mesh ... " << std::flush;
  clock.Start();
  filter->Update();
  clock.Stop();
  std::cout << clock.GetMean() << "s elapsed" << std::endl;



  vtkPolyData* mesh = filter->GetOutput();
  
  if (_OutputFileArg.isSet())
    {
    utl::WriteVtkPolyData(mesh, _OutputFile);
    }
  else
    {
    utlGlobalException(_WindowSize.size()!=2, "wrong window size");
    utlGlobalException(_Angle.size()!=2, "wrong angle size");
    vtk::VisualizePolyDataWithScalarRange(mesh, _ScalarRange, {0.6667,0.0}, _Angle, _WindowSize, !_NoNormalArg.isSet(), _Zoom, _PNGFile);
    }
}
