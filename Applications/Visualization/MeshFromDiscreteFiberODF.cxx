/**
 *       @file  MeshFromDiscreteFiberODF.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-20-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#include "vtkPolyDataWriter.h"

#include "MeshFromDiscreteFiberODFCLP.h"
#include "itkMeshFromDiscreteFiberODFImageFilter.h"
#include "itkTensorBasisMatrixGenerator.h"

#include "itkSpatiallyDenseSparseVectorImageFileReader.h"
#include "itkSpatiallyDenseSparseVectorImage.h"

#include "itkCommandProgressUpdate.h"
#include "utl.h"
#include "utlVTK.h"
#include "vtkPolyDataViewer.h"

/**
 * \brief  Mesh From Discrete Fiber ODF
 *  Create mesh (a spherical profile of DWI/EAP/ODF) from discrete fiber ODF for visualization.
 *
 *  \author  Jian Cheng (jian.cheng.1983@gmail.com)
 */
int 
main (int argc, char const* argv[])
{
  PARSE_ARGS;

  // Time Probe
  itk::TimeProbe clock;

  utlGlobalException(_BoxView.size()!=6, "need 6 parameters in --box");
  utlGlobalException(_SliceView.size()!=3, "need 3 parameters in --box");
  utlGlobalException(_Flip.size()!=3, "need 3 parameters in --flip");
  
  // Define Variables
  typedef double PixelType;
  typedef itk::SpatiallyDenseSparseVectorImage<PixelType, 3> InputImageType;
  typedef itk::SpatiallyDenseSparseVectorImageFileReader<InputImageType> ReaderType;

  InputImageType::Pointer inputImage = InputImageType::New();
  itk::ReadImage<InputImageType, ReaderType>(_InputFile, inputImage);
  
  // Create mesh
  typedef itk::MeshFromDiscreteFiberODFImageFilter<InputImageType> MeshCreatorType;
  MeshCreatorType::Pointer filter = MeshCreatorType::New();
  filter->SetScale( _Scale );
  if (_MaskFileArg.isSet())
    filter->SetMaskImage(_MaskFile);
  
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

  itk::CommandProgressUpdate::Pointer observer =itk::CommandProgressUpdate::New();
  if (_ShowProgressArg.isSet())
    filter->AddObserver( itk::ProgressEvent(), observer );
  if (_SingleThreadArg.isSet())
    filter->SetNumberOfThreads(1);
  filter->SetDebug(_DebugArg.isSet());
  
  itk::BasisMatrixGenerator<double>::MatrixPointer gradBasis;
  if (_OrientationsFileArg.isSet())
    gradBasis = utl::ReadGrad<double>(_OrientationsFile, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  else
    gradBasis = utl::ReadGrad<double>(4, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  utlGlobalException((_OutputType=="DWI" || _OutputType=="EAP") && _Radius<0, "need to set _Radius ");
  if (_Radius>0)
    {
    filter->SetRadius(_Radius);
    }

  utlGlobalException(_TensorEigenValues.size()!=2 && _TensorEigenValues.size()!=3, "wrong size of --tensor");
  
  typedef itk::TensorBasisMatrixGenerator<double> TensorBasisGeneratorType;
  TensorBasisGeneratorType::Pointer tensorBasisGenerator = NULL;
  if (_BasisType=="TENSOR")
    {
    utlGlobalException(_TensorEigenValues[0]<_TensorEigenValues[1]-1e-10, "The first eigenvalue should be the largest");
    tensorBasisGenerator = TensorBasisGeneratorType::New();
    tensorBasisGenerator->SetEigenValue1(_TensorEigenValues[0]);
    tensorBasisGenerator->SetEigenValue2(_TensorEigenValues[1]);
    if (_TensorEigenValues.size()==2)
      tensorBasisGenerator->SetEigenValue3(_TensorEigenValues[1]);
    else
      {
      utlGlobalException(_TensorEigenValues[2]<_TensorEigenValues[3]-1e-10, "The second eigenvalue should be more than the third one");
      tensorBasisGenerator->SetEigenValue3(_TensorEigenValues[2]);
      }
    if (_IsotroipicEigenValue>0)
      tensorBasisGenerator->SetEigenValueISO(_IsotroipicEigenValue);

    tensorBasisGenerator->SetBasisOrientations(gradBasis);
    tensorBasisGenerator->SetUseIsotropicTerm( _UseIsotropicComponentArg.isSet() );
    // tensorBasisGenerator->Print(std::cout<<"basis");
    filter->SetBasisMatrixGenerator(tensorBasisGenerator);
    }
  else if (_BasisType=="CYLINDER")
    {
    utlGlobalException(true, "TODO");
    }

  filter->GetBasisMatrixGenerator()->SetODFOrder(_ODFOrder);
  filter->GetBasisMatrixGenerator()->GetSamplingSchemeQSpace()->SetTau(_Tau);
  filter->GetBasisMatrixGenerator()->GetSamplingSchemeRSpace()->SetTau(_Tau);
  if (_OutputType=="DWI")
    filter->GetBasisMatrixGenerator()->SetOutputType(TensorBasisGeneratorType::DWI);
  else if (_OutputType=="EAP") 
    filter->GetBasisMatrixGenerator()->SetOutputType(TensorBasisGeneratorType::EAP);
  else if (_OutputType=="ODF") 
    filter->GetBasisMatrixGenerator()->SetOutputType(TensorBasisGeneratorType::ODF);

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
