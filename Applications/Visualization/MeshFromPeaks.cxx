/**
 *       @file  MeshFromPeaks.cxx
 *      @brief  
 *     Created  "08-29-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#include "itkMeshFromPeaksImageFilter.h"
#include "MeshFromPeaksCLP.h"

#include "vtkPolyDataWriter.h"

#include "itkCommandProgressUpdate.h"

#include "utl.h"
#include "utlVTK.h"

#include "vtkPolyDataViewer.h"

/**
 * \brief  Gnerate mesh from peaks
 */
  int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;

  typedef double PixelType;
  typedef itk::VectorImage<PixelType, 3> InputImageType;

  // Time Probe
  itk::TimeProbe clock;

  utlGlobalException(_BoxView.size()!=6, "need 6 parameters in --box");
  utlGlobalException(_SliceView.size()!=3, "need 3 parameters in --box");
  utlGlobalException(_Flip.size()!=3, "need 3 parameters in --flip");

  InputImageType::Pointer inputImage;
  itk::ReadVectorImage(_InputFile, inputImage);

  typedef itk::MeshFromPeaksImageFilter<InputImageType> MeshCreatorType;
  MeshCreatorType::Pointer filter = MeshCreatorType::New();

  filter->SetInput( inputImage );
  filter->SetPeakType(itk::PeakContainerHelper::GetPeakType(_PeakType));
  if (_ColorScheme == "DIRECTION") { filter->SetColorScheme(MeshCreatorType::DIRECTION); }
  if (_ColorScheme == "FIXED") { filter->SetColorScheme(MeshCreatorType::FIXED); }
  utlGlobalException(_ColorPeak.size()!=3, "wrong size of _ColorPeak");
  filter->SetColorPeak(_ColorPeak[0], _ColorPeak[1], _ColorPeak[2]);
  filter->SetMaxNumberOfPeaks(_MaxNumber);

  filter->SetScale( _Scale );
  filter->SetTubeRadius( _TubeRadius );

  filter->SetFlip(_Flip[0], _Flip[1], _Flip[2]);
  filter->SetBoxView(_BoxView[0], _BoxView[1], _BoxView[2], _BoxView[3], _BoxView[4], _BoxView[5]);
  filter->SetSliceView(_SliceView[0], _SliceView[1], _SliceView[2]);

  itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  if (_ShowProgressArg.isSet())
    filter->AddObserver( itk::ProgressEvent(), observer );
  if (_SingleThreadArg.isSet())
    filter->SetNumberOfThreads(1);
  filter->SetDebug(_DebugArg.isSet());
  
  std::cout << "Generating mesh ... " << std::flush;
  clock.Start();
  filter->Update();
  clock.Stop();
  std::cout << clock.GetMean() << "s elapsed" << std::endl;
  
  MeshCreatorType::OutputMeshPolyDataType* mesh = filter->GetOutput();
  
  if (_OutputFileArg.isSet())
    {
    utl::WriteVtkPolyData(mesh, _OutputFile);
    }
  else
    {
    utlGlobalException(_WindowSize.size()!=2, "wrong window size");
    utlGlobalException(_Angle.size()!=2, "wrong angle size");
    vtk::VisualizePolyData(mesh, _Angle, _WindowSize, !_NoNormalArg.isSet(), _Zoom, _PNGFile);
    }

  return EXIT_SUCCESS;

  
  return 0;
}
