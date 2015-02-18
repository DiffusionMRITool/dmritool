/*=========================================================================

 Program:   Mesh From Basis Weights Converter

 Copyright (c) Pew-Thian Yap, Jian Cheng. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/


#include "itkTimeProbe.h"
#include "itkMeshFromSHCoefficientsImageFilter.h"
#include "MeshFromSHCoefficientsConverterCLP.h"
#include "vtkPolyDataWriter.h"
#include "utl.h"


int
main(int argc, char *argv[])
{
  PARSE_ARGS;

  // Time Probe
  itk::TimeProbe clock;
  
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
  filter->SetBoxView(_BoxView[0], _BoxView[1], _BoxView[2], _BoxView[3], _BoxView[4], _BoxView[5]);
   
  if (_BasicShape == "TETRAHEDRON") { filter->SetBasicShape(MeshCreatorType::SphereTessellatorType::TETRAHEDRON); }
  if (_BasicShape == "OCTAHEDRON")  { filter->SetBasicShape(MeshCreatorType::SphereTessellatorType::OCTAHEDRON); }
  if (_BasicShape == "ICOSAHEDRON") { filter->SetBasicShape(MeshCreatorType::SphereTessellatorType::ICOSAHEDRON); }
  filter->SetTessellationOrder(_TessellationOrder);

  if (_ColorScheme == "DIRECTION") { filter->SetColorScheme(MeshCreatorType::DIRECTION); }
  if (_ColorScheme == "MAGNITUDE") { filter->SetColorScheme(MeshCreatorType::MAGNITUDE); }

  if (_SingleThreadArg.isSet())
    filter->SetNumberOfThreads(1);
  filter->SetDebug(_DebugArg.isSet());
  
  std::cout << "Generating mesh ... " << std::flush;
  clock.Start();
  filter->Update();
  clock.Stop();
  std::cout << clock.GetMean() << "s elapsed" << std::endl;
  
  MeshCreatorType::OutputMeshPolyDataType* mesh = filter->GetOutput();
  
  // Write Output
  std::cout << "Writing file: " << _OutputFile << std::endl;
  vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  polyDataWriter->SetFileName(_OutputFile.c_str());
  polyDataWriter->SetInputData(mesh);
  polyDataWriter->Write();

  return EXIT_SUCCESS;
}

