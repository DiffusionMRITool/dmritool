/**
 *       @file  MeshFromSphericalFunctionTessellatedSamples.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#include "itkMeshFromSphericalFunctionTessellatedSamplesImageFilter.h"
#include "MeshFromSphericalFunctionTessellatedSamplesCLP.h"
#include "vtkPolyDataWriter.h"
#include "utl.h"


/**
 * \brief  Mesh From spherical samples using a given tessellation
 *  Create mesh (a spherical function) from discrete samples using tessellation.
 *
 *  \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int
main(int argc, char *argv[])
{
  PARSE_ARGS;

  // Define Variables
  typedef float PixelType;
  typedef itk::VectorImage<PixelType, 3> InputImageType;

  InputImageType::Pointer inputImage = InputImageType::New();
  itk::ReadVectorImage(_InputFile, inputImage);

  // Create mesh
  typedef itk::MeshFromSphericalFunctionTessellatedSamplesImageFilter<InputImageType> MeshCreatorType;
  MeshCreatorType::Pointer filter = MeshCreatorType::New();
  filter->SetScale( _Scale );
  
  if (_Normalization=="NONE")  filter->SetNormalization( MeshCreatorType::NONE );
  if (_Normalization=="MIN_MAX")  filter->SetNormalization( MeshCreatorType::MIN_MAX );
  if (_Normalization=="UNIT_MAX")  filter->SetNormalization( MeshCreatorType::UNIT_MAX );
  if (_Normalization=="UNIT_INTEGRAL")  filter->SetNormalization( MeshCreatorType::UNIT_INTEGRAL );
  
  filter->SetPow(_Pow);
  filter->SetRemoveNegativeValues(_RemoveNegativeValuesArg.isSet());
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

  MeshCreatorType::MatrixPointer grad = utl::ReadGrad<double>(_OrientationsFile, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  filter->SetDataOrientations(grad);
  
  std::cout << "Generating mesh ... " << std::flush;
  filter->Update();
  
  MeshCreatorType::OutputMeshPolyDataType* mesh = filter->GetOutput();
  
  // Write Output
  std::cout << "Writing file: " << _OutputFile << std::endl;
  vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  polyDataWriter->SetFileName(_OutputFile.c_str());
  polyDataWriter->SetInputData(mesh);
  polyDataWriter->Write();

  return EXIT_SUCCESS;
}
