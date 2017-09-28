/**
 *       @file  itkFiberTractsWriter.hxx
 *      @brief  
 *     Created  "07-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkFiberTractsWriter_hxx
#define __itkFiberTractsWriter_hxx

#include "itkFiberTractsWriter.h"
#include "utlDMRI.h"

namespace itk
{

void FiberTractsWriter::Update()
{
  int format = utl::GetFiberTractsFormatFromFileExtension(m_FileName); 
  if (format==TRACTS_TRK)
    WriteTractsTRK();
  else if (format==TRACTS_TCK)
    WriteTractsTCK();
  else if (format==TRACTS_VTK)
    {
    utlGlobalException(true, "use MeshFromTracts and itk::MeshFromFiberTractsFilter to write vtk file.");
    // WriteTractsVTK();
    }
  else
    {
    utlGlobalException(true, "un-supported tract format");
    }
}

void FiberTractsWriter::WriteTractsTRK()
{
  auto header = m_FiberTracts->GetHeader();
  auto fibers = m_FiberTracts->GetFibers();
  
  FILE* file;
  file = fopen(m_FileName.c_str(), "wb");
  itk::WriteTrackVisHeader(*header, file);

  int n_count = header->n_count;
  int n_p = header->n_properties;
  int n_s = header->n_scalars;

  long offset=1000;
  fseek (file, offset, SEEK_SET);
  int numPoints;
  float val[3+n_s], val2[n_p];
  VertexType vertex;
  for ( int n = 0; n < n_count; ++n ) 
    {
    FiberPointer fiber = fibers->GetElement(n);

    numPoints = fiber->GetNumberOfPoints();
    fwrite((char*)&numPoints, sizeof(int), 1, file);

    auto tract = fiber->GetTract();
    auto properties = fiber->GetProperties();
    auto scalars = fiber->GetScalars();
    auto vertexList = tract->GetVertexList();
    for ( int i = 0; i < numPoints; ++i ) 
      {
      vertex = (*vertexList)[i];
      const STDVectorType& ss = (*scalars)[i];
      for ( int j = 0; j < 3; ++j ) 
          val[j] = vertex[j];
      for ( int j = 0; j < n_s; ++j ) 
          val[j+3] = ss[j];
      fwrite((char*)&val, sizeof(float)*(3+n_s), 1, file);
      }
    for ( int j = 0; j < n_p; ++j ) 
      val2[j+3] = (*properties)[j];
    fwrite((char*)&val2, sizeof(float)*n_p, 1, file);
    }

  fclose (file);
}

void FiberTractsWriter::WriteTractsTCK()
{
  utlGlobalException(true, "TODO: read tracts in tck format");
}

// void FiberTractsWriter::WriteTractsVTK()
// {
// }

}


#endif 



