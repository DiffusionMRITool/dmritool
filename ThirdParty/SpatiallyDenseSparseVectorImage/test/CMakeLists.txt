itk_module_test()
set(ITKSpatiallyDenseSparseVectorImageTests
  itkVectorToSparseVectorImageTest.cxx
  itkSparseVectorToVectorImageTest.cxx
)

CreateTestDriver(ITKSpatiallyDenseSparseVectorImage "${ITKSparseVectorImage-Test_LIBRARIES}" "${ITKSparseVectorImageTests}")

# Test I/O convertor
itk_add_test( NAME itkVectorToSparseVectorImageTest
  COMMAND ITKSpatiallyDenseSparseVectorImageTestDriver
  --compare DATA{Baseline/itkSpatiallyDenseSparseVectorImageTestBaseline_Sparse_key.nrrd} ${ITK_TEST_OUTPUT_DIR}/itkSpatiallyDenseSparseVectorImageTestOutput_Sparse_key.nrrd
  --compare DATA{Baseline/itkSpatiallyDenseSparseVectorImageTestBaseline_Sparse_value.nrrd} ${ITK_TEST_OUTPUT_DIR}/itkSpatiallyDenseSparseVectorImageTestOutput_Sparse_value.nrrd
  itkVectorToVectorImageTest DATA{Baseline/itkSpatiallyDenseSparseVectorImageTestBaseline_Dense.nii.gz} ${ITK_TEST_OUTPUT_DIR}/itkSpatiallyDenseSparseVectorImageTestOutput_Sparse.spr
  )

itk_add_test( NAME itkSparseVectorToVectorImageTest
  COMMAND ITKSpatiallyDenseSparseVectorImageTestDriver
  --compare DATA{Baseline/itkSpatiallyDenseSparseVectorImageTestBaseline_Dense.nii.gz} ${ITK_TEST_OUTPUT_DIR}/itkSpatiallyDenseSparseVectorImageTestOutput_Dense.nii.gz
  itkSparseVectorToVectorImageTest DATA{Baseline/itkSpatiallyDenseSparseVectorImageTestBaseline_Sparse.spr} ${ITK_TEST_OUTPUT_DIR}/itkSpatiallyDenseSparseVectorImageTestOutput_Dense.nii.gz
  )

