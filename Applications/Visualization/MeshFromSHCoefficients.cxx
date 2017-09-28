/*=========================================================================

 Program:   Mesh From Basis Weights 

 Copyright (c) Pew-Thian Yap, Jian Cheng. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/


#include "itkMeshFromSHCoefficientsImageFilter.h"
#include "MeshFromSHCoefficientsCLP.h"

#include "vtkPolyDataWriter.h"

#include "itkCommandProgressUpdate.h"

#include "utl.h"
#include "utlVTK.h"

#include "vtkPolyDataViewer.h"

int
main(int argc, char *argv[])
{
  PARSE_ARGS;

  // Time Probe
  itk::TimeProbe clock;

  utlGlobalException(_BoxView.size()!=6, "need 6 parameters in --box");
  utlGlobalException(_SliceView.size()!=3, "need 3 parameters in --box");
  utlGlobalException(_Flip.size()!=3, "need 3 parameters in --flip");
  
  // Define Variables
  typedef double PixelType;
  typedef itk::VectorImage<PixelType, 3> InputImageType;

  InputImageType::Pointer inputImage;
  itk::ReadVectorImage(_InputFile, inputImage);

  // Compute MaxOrder if not set
  if ( !_MaxOrderArg.isSet() )
    _MaxOrder = utl::DimToRankSH(inputImage->GetNumberOfComponentsPerPixel());
  
  // Check Parameters
  unsigned numberOfBasisFunctions = (_MaxOrder + 1) * (_MaxOrder + 2) / 2;
  itkAssertOrThrowMacro(numberOfBasisFunctions <=inputImage->GetNumberOfComponentsPerPixel(), "MaxOrder too large. numberOfBasisFunctions="<< numberOfBasisFunctions 
      <<", inputImage->GetNumberOfComponentsPerPixel()="<<inputImage->GetNumberOfComponentsPerPixel());
  
  // Create mesh
  typedef itk::MeshFromSHCoefficientsImageFilter<InputImageType> MeshCreatorType;
  MeshCreatorType::Pointer filter = MeshCreatorType::New();
  filter->SetScale( _Scale );
  
  if (_Normalization=="NONE")  filter->SetNormalization( MeshCreatorType::NONE );
  if (_Normalization=="MIN_MAX")  filter->SetNormalization( MeshCreatorType::MIN_MAX );
  if (_Normalization=="UNIT_MAX")  filter->SetNormalization( MeshCreatorType::UNIT_MAX );
  if (_Normalization=="UNIT_INTEGRAL")  filter->SetNormalization( MeshCreatorType::UNIT_INTEGRAL );
  
  filter->SetPow(_Pow);
  filter->SetRemoveNegativeValues(_RemoveNegativeValuesArg.isSet());
  filter->SetMaxOrder( _MaxOrder );
  filter->SetInput( inputImage );

  filter->SetFlip(_Flip[0], _Flip[1], _Flip[2]);
  filter->SetBoxView(_BoxView[0], _BoxView[1], _BoxView[2], _BoxView[3], _BoxView[4], _BoxView[5]);
  filter->SetSliceView(_SliceView[0], _SliceView[1], _SliceView[2]);
   
  if (_BasicShape == "TETRAHEDRON") { filter->SetBasicShape(MeshCreatorType::SphereTessellatorType::TETRAHEDRON); }
  if (_BasicShape == "OCTAHEDRON")  { filter->SetBasicShape(MeshCreatorType::SphereTessellatorType::OCTAHEDRON); }
  if (_BasicShape == "ICOSAHEDRON") { filter->SetBasicShape(MeshCreatorType::SphereTessellatorType::ICOSAHEDRON); }
  filter->SetTessellationOrder(_TessellationOrder);

  if (_ColorScheme == "DIRECTION") { filter->SetColorScheme(MeshCreatorType::DIRECTION); }
  if (_ColorScheme == "MAGNITUDE") { filter->SetColorScheme(MeshCreatorType::MAGNITUDE); }

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
}

