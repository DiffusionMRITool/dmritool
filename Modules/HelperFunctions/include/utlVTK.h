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

#include "utlVTKMacro.h"

namespace utl
{


/** Write polydata to file */
inline void 
WriteVtkPolyData ( vtkPolyData* mesh, const std::string& filename)
{
  std::cout << "Writing polydata to " << filename << std::endl;
  vtkSmartPointer<vtkPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  polyDataWriter->SetFileName(filename.c_str());
  polyDataWriter->SetInputData(mesh);
  polyDataWriter->Write();
}

}


#endif 
