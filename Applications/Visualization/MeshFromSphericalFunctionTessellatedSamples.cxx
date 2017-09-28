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
#include "utlVTK.h"
#include "vtkPolyDataViewer.h"

#include "itkSpatiallyDenseSparseVectorImage.h"
#include "itkSpatiallyDenseSparseVectorImageFileReader.h"


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

  utlGlobalException(_BoxView.size()!=6, "need 6 parameters in --box");
  utlGlobalException(_SliceView.size()!=3, "need 3 parameters in --box");
  utlGlobalException(_Flip.size()!=3, "need 3 parameters in --flip");

  // Define Variables
  typedef float PixelType;
  typedef itk::VectorImage<PixelType, 3> InputImageType;

  InputImageType::Pointer inputImage = InputImageType::New();
  if (itk::IsSparseImage(_InputFile))
    {
    typedef itk::SpatiallyDenseSparseVectorImage<PixelType, 3> InputSparseImageType;

    InputSparseImageType::Pointer inputSparseImage = InputSparseImageType::New();
    typedef itk::SpatiallyDenseSparseVectorImageFileReader<InputSparseImageType> ReaderType;

    itk::ReadImage<InputSparseImageType, ReaderType>(_InputFile, inputSparseImage);
    itk::ImageToImage<InputSparseImageType, InputImageType>(inputSparseImage, inputImage);
    }
  else
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
  filter->SetSliceView(_SliceView[0], _SliceView[1], _SliceView[2]);
  filter->SetFlip(_Flip[0], _Flip[1], _Flip[2]);
   
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
