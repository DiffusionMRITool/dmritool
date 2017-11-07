/**
 *       @file  MeshFromTensors.cxx
 *      @brief  
 *     Created  "06-19-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#include "itkMeshFromTensorImageFilter.h"
#include "MeshFromTensorsCLP.h"

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

  if (utl::IsEndingWith(_InputFile, ".txt"))
    {
    std::vector<double> tensor;
    utl::ReadVector(_InputFile, tensor, " \t,");
    InputImageType::PixelType pixel(tensor.size());
    utl::VectorToVector(tensor, pixel, tensor.size());
    inputImage = itk::GenerateImageFromSingleVoxel<InputImageType>(pixel);
    }
  else if (utl::IsFileExist(_InputFile))
    itk::ReadVectorImage(_InputFile, inputImage);
  else
    {
    std::vector<double> tensor, tensor6d(6), tensor9d(9);
    utl::SetVector(_InputFile, tensor);
    utlGlobalException(tensor.size()!=6 && tensor.size()!=9, "size should be 6 or 9");
    if (tensor.size()==9)
      {
      tensor9d = tensor;
      utl::ConvertTensor9DTo6D(tensor9d, tensor6d, TENSOR_UPPER_TRIANGULAR);
      }
    else
      {
      tensor6d = tensor;
      utl::ConvertTensor6DTo9D(tensor6d, tensor9d, TENSOR_UPPER_TRIANGULAR);
      }
    utl::PrintVector(tensor9d, "tensor");
    InputImageType::PixelType pixel(tensor6d.size());
    utl::VectorToVector(tensor6d, pixel, tensor6d.size());
    inputImage = itk::GenerateImageFromSingleVoxel<InputImageType>(pixel);
    }

  typedef itk::MeshFromTensorImageFilter<InputImageType> MeshCreatorType;
  MeshCreatorType::Pointer filter = MeshCreatorType::New();

  filter->SetInput( inputImage );
  
  typedef MeshCreatorType::ScalarImageType ScalarImageType;
  ScalarImageType::Pointer scalarImage = ScalarImageType::New();
  if (_ScalarImageFileArg.isSet())
    {
    itk::ReadImage(_ScalarImageFile, scalarImage);
    filter->SetScalarImage(scalarImage);
    }

  if (_GlyphType == "LINE") { filter->SetShapeMode(MeshCreatorType::GLYPH_LINE); }
  if (_GlyphType == "ARROW") { filter->SetShapeMode(MeshCreatorType::GLYPH_ARROW); }
  if (_GlyphType == "DISK") { filter->SetShapeMode(MeshCreatorType::GLYPH_DISK); }
  if (_GlyphType == "CYLINDER") { filter->SetShapeMode(MeshCreatorType::GLYPH_CYLINDER); }
  if (_GlyphType == "CUBE") { filter->SetShapeMode(MeshCreatorType::GLYPH_CUBE); }
  if (_GlyphType == "ELLIPSOID") { filter->SetShapeMode(MeshCreatorType::GLYPH_SPHERE); }
  if (_GlyphType == "SUPERQUADRIC") { filter->SetShapeMode(MeshCreatorType::GLYPH_SUPERQUADRIC); }
  
  if (_ColorScheme == "NONE") { filter->SetTensorColorScheme(MeshCreatorType::COLOR_NONE); }
  if (_ColorScheme == "FA") { filter->SetTensorColorScheme(MeshCreatorType::COLOR_BY_FA); }
  if (_ColorScheme == "MD") { filter->SetTensorColorScheme(MeshCreatorType::COLOR_BY_MD); }
  if (_ColorScheme == "DIRECTION") { filter->SetTensorColorScheme(MeshCreatorType::COLOR_BY_DIRECTION); }
  if (_ColorScheme == "IMAGE") { filter->SetTensorColorScheme(MeshCreatorType::COLOR_BY_IMAGE); }

  filter->SetScale( _Scale );
  filter->SetFlip(_Flip[0], _Flip[1], _Flip[2]);
  filter->SetBoxView(_BoxView[0], _BoxView[1], _BoxView[2], _BoxView[3], _BoxView[4], _BoxView[5]);
  filter->SetSliceView(_SliceView[0], _SliceView[1], _SliceView[2]);

  // itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  // if (_ShowProgressArg.isSet())
  //   filter->AddObserver( itk::ProgressEvent(), observer );
  filter->SetGlyphResolution(_GlyphResolution);
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
    utlGlobalException(_BackgroundColor.size()!=3, "wrong size of background color");
    if (_ColorScheme=="DIRECTION")
      {
      // huerange (0,1) used for direction coloring using RGBToIndex
      vtk::VisualizePolyDataWithScalarRange(mesh, _ScalarRange, {0.0,1.0}, _Angle, _WindowSize, !_NoNormalArg.isSet(), !_NoLightingArg.isSet(), _Zoom, _PNGFile, _BackgroundColor);
      }
    else
      {
      // huerange (0.0667,0) used for scalars 
      vtk::VisualizePolyDataWithScalarRange(mesh, _ScalarRange, {0.6667,0.0}, _Angle, _WindowSize, !_NoNormalArg.isSet(), !_NoLightingArg.isSet(), _Zoom, _PNGFile, _BackgroundColor);
      }
    }

  return EXIT_SUCCESS;
}
