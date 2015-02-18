file( DOWNLOAD http://public.kitware.com/pub/itk/NightlyDoxygen/InsightDoxygenDocTag.gz
  ${ITKTAG_DIR}/InsightDoxygen.tag.gz SHOW_PROGRESS
  )

file( DOWNLOAD http://www.vtk.org/files/nightly/vtkNightlyDoc.tag.gz
  ${VTKTAG_DIR}/vtkNightlyDoc.tag.gz SHOW_PROGRESS
  )
