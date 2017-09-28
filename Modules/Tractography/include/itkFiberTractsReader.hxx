/**
 *       @file  itkFiberTractsReader.hxx
 *      @brief  
 *     Created  "08-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkFiberTractsReader_hxx
#define __itkFiberTractsReader_hxx

#include "itkFiberTractsReader.h"
#include "utlDMRI.h"

namespace itk
{

void FiberTractsReader::Update()
{
  utlGlobalException(!utl::IsFileExist(m_FileName), m_FileName + " does not exist");
  int format = utl::GetFiberTractsFormat(m_FileName); 
  if (format==TRACTS_TRK)
    ReadTractsTRK();
  else if (format==TRACTS_TCK)
    ReadTractsTCK();
  else if (format==TRACTS_VTK)
    ReadTractsVTK();
  else
    {
    utlGlobalException(true, "un-supported tract format");
    }
}

void FiberTractsReader::ReadTractsTRK()
{
  auto header = m_FiberTracts->GetHeader();
  auto fibers = m_FiberTracts->GetFibers();
  
  itk::ReadTrackVisHeader(m_FileName, *header);

  int n_count = m_FiberTracts->GetNumberOfFibers();
  int n_p = header->n_properties;
  int n_s = header->n_scalars;

  FILE* file;
  file = fopen(m_FileName.c_str(), "rb");
  fseek (file, 0, SEEK_END);  
  long fsize = ftell(file);
  long offset=1000;
  int numPoints;
  float val[3+n_s], val2[n_p];
  VertexType vertex;
  while (offset<fsize)
    {
    fseek (file, offset, SEEK_SET);

    fread((char*)&numPoints, sizeof(int), 1, file);

    FiberPointer fiber = FiberType::New();
    auto tract = fiber->GetTract();
    auto properties = fiber->GetProperties();
    auto scalars = fiber->GetScalars();
    for ( int i = 0; i < numPoints; ++i ) 
      {
      fread((char*)&val, sizeof(float)*(3+n_s), 1, file);
      vertex[0]=val[0];
      vertex[1]=val[1];
      vertex[2]=val[2];
      tract->AddVertex(vertex);
      STDVectorType vec;
      for ( int j = 0; j < n_s; ++j ) 
        vec.push_back(val[3+j]);
      scalars->push_back(vec);
      }

    fread((char*)&val2, sizeof(float)*n_p, 1, file);
    for ( int j = 0; j < n_p; ++j ) 
      properties->push_back(val[j]);

    fibers->InsertElement( fibers->Size(), fiber);

    offset += 4+ numPoints*(3+n_s)*4 + n_p*4;
    }

  fclose (file);
}

void FiberTractsReader::ReadTractsTCK()
{
  utlGlobalException(true, "TODO: read tracts in tck format");
}

void FiberTractsReader::ReadTractsVTK()
{
  utlGlobalException(true, "TODO: read tracts in vtk format");
}

}


#endif 


