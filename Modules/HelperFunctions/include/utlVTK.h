/**
 *       @file  utlVTK.h
 *      @brief  
 *     Created  "06-06-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __utlVTK_h
#define __utlVTK_h

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>

#include "utlVTKMacro.h"
#include "utlCore.h"

namespace utl
{

template <typename WriterType>
inline void 
WriteVTK(vtkPolyData* mesh, const std::string& filename)
{
  vtkSmartPointer< WriterType > writer = vtkSmartPointer< WriterType >::New();
  writer->SetInputData(mesh);
  writer->SetFileName(filename.c_str());
  writer->Write();
}



/** Write polydata to file */
inline void 
WriteVtkPolyData ( vtkPolyData* mesh, const std::string& filename)
{
  std::string fileNoExt, ext;
  utl::GetFileExtension(filename, ext, fileNoExt);
  std::cout << "Writing polydata to " << filename << std::endl;
  if (ext=="vtk")
    WriteVTK<vtkPolyDataWriter>(mesh, filename);
  else if (ext=="vtp")
    WriteVTK<vtkXMLPolyDataWriter>(mesh, filename);
  else if (ext=="ply")
    WriteVTK<vtkPLYWriter>(mesh, filename);
  else if (ext=="stl")
    WriteVTK<vtkSTLWriter>(mesh, filename);
  else if (ext=="vtu")
    WriteVTK<vtkXMLUnstructuredGridWriter>(mesh, filename);
  else if (ext=="vti")
    WriteVTK<vtkXMLImageDataWriter>(mesh, filename);
  else if (ext=="vts")
    WriteVTK<vtkXMLStructuredGridWriter>(mesh, filename);
  else if (ext=="vtr")
    WriteVTK<vtkXMLRectilinearGridWriter>(mesh, filename);
  else
    WriteVTK<vtkPolyDataWriter>(mesh, filename);
}

}


#endif 
