/**
 *       @file  SphereTessellator.cxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  10-30-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */  
#include <SphereTessellatorCLP.h> 
#include "utl.h"

#include "itkSphereTessellator.h"


int 
main (int argc, char *argv[]) 
{ 
  PARSE_ARGS;

  typedef double TElement;
  typedef itk::SphereTessellator<TElement> TessellatorType;
  typedef itk::SphereTessellator<TElement>::Pointer TessellatorPointer;
  
  TessellatorPointer tessellator = TessellatorType::New();

  if ( _BasicShape == "TETRAHEDRON" )
    tessellator->SetBasicShape(TessellatorType::TETRAHEDRON);
  if ( _BasicShape == "OCTAHEDRON" )
    tessellator->SetBasicShape(TessellatorType::OCTAHEDRON);
  if ( _BasicShape == "ICOSAHEDRON" )
    tessellator->SetBasicShape(TessellatorType::ICOSAHEDRON);

  tessellator->SetOrder(_Order);
  tessellator->TessellateSphere();

  unsigned int numberOfVertices = tessellator->GetNumberOfVertices();
  unsigned int numberOfFaces = tessellator->GetNumberOfFaces();
  unsigned int numberOfEdges = tessellator->GetNumberOfEdges();
  
  std::ofstream ofs;                          // create ofstream object
  
  ofs.open ( _OutputFile.c_str() );           // open ofstream
  utlAssert (ofs, "\nERROR : failed to open output file " << _OutputFile );

  TessellatorType::PointsMatrixType points = tessellator->GetPointsMatrix();
  TessellatorType::CellsMatrixType cells = tessellator->GetCellsMatrix();

  // if (_OutInOderrArg.isSet())
  //   {
  //   TessellatorType::PointsMatrixType points_back = points;
  //   TessellatorType::CellsMatrixType cells_back = cells;
  //   std::vector<bool> used(points.rows(),false);
  //   int count=0;
  //   for ( int i = 0; i < cells.rows(); i += 1 ) 
  //     {
  //     for ( int j = 0; j < cells.columns(); j += 1 ) 
  //       {
  //       if (used[cells_back(i,j)])
  //         continue;

  //       used[cells_back(i,j)] = true;
  //       points.set_row(count, points_back.get_row(cells_back(i,j)));
  //       for ( int ii = 0; ii < cells.rows(); ii += 1 ) 
  //         for ( int jj = 0; jj < cells.columns(); jj += 1 ) 
  //           {
  //           if (cells_back(ii,jj)==cells_back(i,j))
  //             cells(ii,jj) = count;
  //           }
  //       count++;
  //       }
  //     }
  //   }

  ofs.precision(10);
  if (_NumberOfDirectionsArg.isSet() && !_HemisphereArg.isSet())
    ofs << numberOfVertices << " " << numberOfFaces << " " << numberOfEdges << "\n";
  if (_NumberOfDirectionsArg.isSet() && _HemisphereArg.isSet())
    ofs << numberOfVertices/2 << "\n";

  if (!_HemisphereArg.isSet())
    {
    for ( int i = 0; i < numberOfVertices; i += 1 )
      ofs << points(i,0) << " " << points(i,1) << " " << points(i,2) << "\n";
    }
  else
    {
    for ( int i = 0; i < numberOfVertices; i += 1 )
      {
      if ( points(i,2)>0
        || (std::abs(points(i,2))<1e-20 && points(i,1)<0)
        || (std::abs(points(i,2))<1e-20 && std::abs(points(i,1))<1e-20 && points(i,0)>0))
        ofs << points(i,0) << " " << points(i,1) << " " << points(i,2) << "\n";
      }
    }

  if (_CellArg.isSet() && !_HemisphereArg.isSet())
    {
    for ( int i = 0; i <numberOfFaces; i += 1 )
      ofs << cells(i,0) << " " << cells(i,1) << " " << cells(i,2) << "\n";
    }

  ofs.close();

  return EXIT_SUCCESS; 
} 
