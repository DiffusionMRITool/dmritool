/*=========================================================================

 Program:   Sphere Tessellator
 
 Copyright (c) Pew-Thian Yap. All rights reserved.
 See http://www.unc.edu/~ptyap/ for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/

#ifndef _itkSphereTessellator_hxx
#define _itkSphereTessellator_hxx

#include <cstdio> 
#include <cstdlib> 
#include <cstring> 
#include <cmath>
 
#include "itkSphereTessellator.h"

namespace itk
{

template <typename TElement>
SphereTessellator< TElement >
::SphereTessellator()
{
  m_BasicShape = ICOSAHEDRON;
  m_Order = 1;
  
  vertices = NULL;
  faces = NULL; 
  start = NULL; 
  end = NULL; 
  midpoint = NULL; 
  
  m_NumberOfVertices = 0;
  m_NumberOfFaces = 0;
  m_NumberOfEdges = 0;

  m_Points.clear();
  m_Cells.clear();
}
  
template <typename TElement>
SphereTessellator< TElement >
::~SphereTessellator()
{
  if (vertices) free (vertices); 
  if (faces) free (faces); 
}

/**
 * Tell the container to allocate enough memory to allow at least
 * as many elements as the size given to be stored.
 */
template <typename TElement>
void
SphereTessellator< TElement >
::Reserve(void)
{
  // Do nothing; vnl_matrix cannot reserve elements
}


/**
 * Tell the container to try to minimize its memory usage for storage of
 * the current number of elements.
 */
template <typename TElement>
void
SphereTessellator< TElement >
::Squeeze(void)
{
  // Do nothing; vnl_matrix cannot be squeezed
}


/**
 * Tell the container to release any of its allocated memory.
 */
template <typename TElement>
void
SphereTessellator< TElement >
::Initialize(void)
{
  m_Points.clear();
  m_Cells.clear();
  
  this->Modified();
}

template <typename TElement>
void
SphereTessellator< TElement >
::InitTetrahedron (void) 
{ 
  double sqrt3 = 1 / vcl_sqrt(3.0);
  double tetrahedron_vertices[] = {sqrt3, sqrt3, sqrt3,
          -sqrt3, -sqrt3, sqrt3,
          -sqrt3, sqrt3, -sqrt3,
          sqrt3, -sqrt3, -sqrt3}; 
  int tetrahedron_faces[] = {0, 2, 1, 0, 1, 3, 2, 3, 1, 3, 2, 0};

  n_vertices = 4; 
  n_faces = 4; 
  n_edges = 6; 
  vertices = (double*)malloc(3*n_vertices*sizeof(double)); 
  faces = (int*)malloc(3*n_faces*sizeof(int)); 
  memcpy ((void*)vertices, (void*)tetrahedron_vertices, 3*n_vertices*sizeof(double)); 
  memcpy ((void*)faces, (void*)tetrahedron_faces, 3*n_faces*sizeof(int)); 
} 

template <typename TElement>
void
SphereTessellator< TElement >
::InitOctahedron (void) 
{ 
  double octahedron_vertices[] = {0.0, 0.0, -1.0,
         1.0, 0.0, 0.0,
         0.0, -1.0, 0.0,
         -1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0}; 
  int octahedron_faces[] = {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1, 5, 2, 1, 5, 3, 2, 5, 4, 3, 5, 1, 4}; 

  n_vertices = 6; 
  n_faces = 8;
  n_edges = 12; 
  vertices = (double*)malloc(3*n_vertices*sizeof(double)); 
  faces = (int*)malloc(3*n_faces*sizeof(int)); 
  memcpy ((void*)vertices, (void*)octahedron_vertices, 3*n_vertices*sizeof(double)); 
  memcpy ((void*)faces, (void*)octahedron_faces, 3*n_faces*sizeof(int)); 
} 

template <typename TElement>
void
SphereTessellator< TElement >
::InitIcosahedron (void) 
{ 
  double t = (1+vcl_sqrt(5))/2;
  double tau = t/vcl_sqrt(1+t*t);
  double one = 1/vcl_sqrt(1+t*t);

  double icosahedron_vertices[] = {tau, one, 0.0,
          -tau, one, 0.0,
          -tau, -one, 0.0,
          tau, -one, 0.0,
          one, 0.0 ,  tau,
          one, 0.0 , -tau,
          -one, 0.0 , -tau,
          -one, 0.0 , tau,
          0.0 , tau, one,
          0.0 , -tau, one,
          0.0 , -tau, -one,
          0.0 , tau, -one};
 int icosahedron_faces[] = {4, 8, 7,
          4, 7, 9,
          5, 6, 11,
          5, 10, 6,
          0, 4, 3,
          0, 3, 5,
          2, 7, 1,
          2, 1, 6,
          8, 0, 11,
          8, 11, 1,
          9, 10, 3,
          9, 2, 10,
          8, 4, 0,
          11, 0, 5,
          4, 9, 3,
          5, 3, 10,
          7, 8, 1,
          6, 1, 11,
          7, 2, 9,
          6, 10, 2};
 
  n_vertices = 12; 
  n_faces = 20;
  n_edges = 30;
  vertices = (double*)malloc(3*n_vertices*sizeof(double)); 
  faces = (int*)malloc(3*n_faces*sizeof(int)); 
  memcpy ((void*)vertices, (void*)icosahedron_vertices, 3*n_vertices*sizeof(double)); 
  memcpy ((void*)faces, (void*)icosahedron_faces, 3*n_faces*sizeof(int)); 
} 

template <typename TElement>
int
SphereTessellator< TElement >
::SearchMidpoint (int index_start, int index_end) 
{ 
  unsigned int i;
  for (i=0; i<edge_walk; i++) 
    if ((start[i] == index_start && end[i] == index_end) || 
      (start[i] == index_end && end[i] == index_start)) 
      {
      int res = midpoint[i];

      /* update the arrays */
      start[i]    = start[edge_walk-1];
      end[i]      = end[edge_walk-1];
      midpoint[i] = midpoint[edge_walk-1];
      edge_walk--;

      return res; 
      }

  /* vertex not in the list, so we add it */
  start[edge_walk] = index_start;
  end[edge_walk] = index_end; 
  midpoint[edge_walk] = n_vertices; 

  /* create new vertex */ 
  vertices[3*n_vertices]   = (vertices[3*index_start] + vertices[3*index_end]) / 2.0;
  vertices[3*n_vertices+1] = (vertices[3*index_start+1] + vertices[3*index_end+1]) / 2.0;
  vertices[3*n_vertices+2] = (vertices[3*index_start+2] + vertices[3*index_end+2]) / 2.0;

  /* normalize the new vertex */ 
  double length = vcl_sqrt (vertices[3*n_vertices] * vertices[3*n_vertices] +
    vertices[3*n_vertices+1] * vertices[3*n_vertices+1] +
    vertices[3*n_vertices+2] * vertices[3*n_vertices+2]);
  length = 1/length;
  vertices[3*n_vertices] *= length;
  vertices[3*n_vertices+1] *= length;
  vertices[3*n_vertices+2] *= length;

  n_vertices++;
  edge_walk++;
  return midpoint[edge_walk-1];
} 

template <typename TElement>
void
SphereTessellator< TElement >
::Subdivide (void) 
{ 
  unsigned int n_vertices_new = n_vertices+2*n_edges; 
  unsigned int n_faces_new = 4*n_faces; 
  unsigned int i; 

  edge_walk = 0; 
  n_edges = 2*n_vertices + 3*n_faces; 
  start = (int*)malloc(n_edges*sizeof (int)); 
  end = (int*)malloc(n_edges*sizeof (int)); 
  midpoint = (int*)malloc(n_edges*sizeof (int)); 

  int *faces_old = (int*)malloc (3*n_faces*sizeof(int)); 
  faces_old = (int*)memcpy((void*)faces_old, (void*)faces, 3*n_faces*sizeof(int)); 
  vertices = (double*)realloc ((void*)vertices, 3*n_vertices_new*sizeof(double)); 
  faces = (int*)realloc ((void*)faces, 3*n_faces_new*sizeof(int)); 
  n_faces_new = 0; 

  for (i=0; i<n_faces; i++) 
    { 
      int a = faces_old[3*i]; 
      int b = faces_old[3*i+1]; 
      int c = faces_old[3*i+2]; 

      int ab_midpoint = SearchMidpoint (b, a); 
      int bc_midpoint = SearchMidpoint (c, b); 
      int ca_midpoint = SearchMidpoint (a, c); 

      faces[3*n_faces_new] = a; 
      faces[3*n_faces_new+1] = ab_midpoint; 
      faces[3*n_faces_new+2] = ca_midpoint; 
      n_faces_new++; 
      faces[3*n_faces_new] = ca_midpoint; 
      faces[3*n_faces_new+1] = ab_midpoint; 
      faces[3*n_faces_new+2] = bc_midpoint; 
      n_faces_new++; 
      faces[3*n_faces_new] = ca_midpoint; 
      faces[3*n_faces_new+1] = bc_midpoint; 
      faces[3*n_faces_new+2] = c; 
      n_faces_new++; 
      faces[3*n_faces_new] = ab_midpoint; 
      faces[3*n_faces_new+1] = b; 
      faces[3*n_faces_new+2] = bc_midpoint; 
      n_faces_new++; 
    } 
  n_faces = n_faces_new; 
  free (start); 
  free (end); 
  free (midpoint); 
  free (faces_old); 
} 

/**
 * Do actual work in generating the basis.
 */
template <typename TElement>
void
SphereTessellator< TElement >
::TessellateSphere(void)
{
  switch (m_BasicShape)
    {
  case TETRAHEDRON:
    InitTetrahedron ();
    break;
  case OCTAHEDRON:
    InitOctahedron ();
    break;
  case ICOSAHEDRON:
    InitIcosahedron ();
    break;
    }

  for (unsigned int k=1; k < m_Order; k++) 
    Subdivide ();

  OutputSphere ();
}

template <typename TElement>
void
SphereTessellator< TElement >
::OutputSphere () 
{ 
//  std::cout fprintf (ptr, "OFF\n%d %d %d\n", n_vertices, n_faces, n_edges); 
  m_NumberOfVertices = n_vertices;
  m_NumberOfFaces = n_faces;
  m_NumberOfEdges = n_edges;
  
  m_Points.set_size(n_vertices, 3);
  m_Cells.set_size(n_faces, 3);
  
  for (unsigned int k=0; k<n_vertices; k++)
    {
    m_Points(k,0) = vertices[3*k];
    m_Points(k,1) = vertices[3*k+1];
    m_Points(k,2) = vertices[3*k+2];
    }
  for (unsigned int k=0; k<n_faces; k++)
    {
    m_Cells(k,0) = faces[3*k];
    m_Cells(k,1) = faces[3*k+1];
    m_Cells(k,2) = faces[3*k+2];
    }
} 

template <typename TElement>
void
SphereTessellator< TElement >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "tessellation order = " << m_Order << std::endl;
}

} // end namespace itk

#endif
