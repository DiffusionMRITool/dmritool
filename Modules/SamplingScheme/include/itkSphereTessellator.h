/*=========================================================================

  Program:   Sphere Tessellator
  Module:    $HeadURL $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSphereTessellator_h
#define __itkSphereTessellator_h

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace itk
{

/** \class SphereTessellator
 *  \brief Tesselates a sphere via subdivision of tetrahedron, octahedron, or icosahedron.
 *
 * Generation of a sphere by recursive subdivisions
 * inspired by source code written by Jon Leech
 * (http://www.cs.unc.edu/~jon/sphere.html)
 *
 * \author Pew-Thian Yap, UNC-CH
 * \ingroup DWISamplingScheme
 */

template <typename TElement = double>
class ITK_EXPORT SphereTessellator: public Object
{
public:
  /** Standard class typedefs. */
  typedef SphereTessellator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(SphereTessellator, Object);

  /** Set the orientation tables. */
  const static unsigned int NDimensions = 3;

  /** Save the template parameters. */
  typedef vnl_matrix<TElement> PointsMatrixType;
  typedef vnl_matrix<unsigned long> CellsMatrixType;

  /** Generate basis. */
  void TessellateSphere();

  /** Get/Set basic shape for subdivision */
  typedef enum {TETRAHEDRON, OCTAHEDRON, ICOSAHEDRON} BasicShapeType;
  
  itkGetMacro(BasicShape, BasicShapeType);
  itkSetMacro(BasicShape, BasicShapeType);

  /** Get/Set tesselation order
   oder = numer of subdivision + 1  */
  itkGetMacro(Order, unsigned int);
  itkSetMacro(Order, unsigned int);

  /** Information on tesselation */
  itkGetMacro(NumberOfVertices, unsigned int);
  itkGetMacro(NumberOfFaces, unsigned int);
  itkGetMacro(NumberOfEdges, unsigned int);
  
  PointsMatrixType GetPointsMatrix()
  {
    return m_Points;
  }
  
  PointsMatrixType GetPointsMatrixInHemisphere()
  {
  PointsMatrixType points(m_Points.rows()/2, 3);
  int j=0;
  for ( int i = 0; i < m_NumberOfVertices; i += 1 )
    {
    if ( points(i,2)>0
      || (std::abs(points(i,2))<1e-20 && points(i,1)<0)
      || (std::abs(points(i,2))<1e-20 && std::abs(points(i,1))<1e-20 && points(i,0)>0))
      {
      for ( int kk = 0; kk < 3; ++kk ) 
        points(j,3) = m_Points(i,3);
      j++;
      }
    }
  return points;
  }

  CellsMatrixType GetCellsMatrix()
  {
    return m_Cells;
  }

  /** Release allocated memory. */
  void Initialize(void);
  
  /** This method does nothing. */
  void Reserve(void);

  /** This method does nothing. */
  void Squeeze(void);

protected:
  SphereTessellator();
  virtual ~SphereTessellator();

  /** PrintSelf routine. Normally this is a protected internal method. It is
   * made public here so that Image can call this method.  Users should not
   * call this method but should call Print() instead. */
  void PrintSelf(std::ostream& os, Indent indent) const;

  void InitTetrahedron();
  void InitOctahedron();
  void InitIcosahedron();
  int  SearchMidpoint (int index_start, int index_end);
  void Subdivide() ;
  void OutputSphere();
  
private:
  SphereTessellator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int n_vertices;
  unsigned int n_faces;
  unsigned int n_edges;
  double *vertices;
  int *faces; 

  unsigned int edge_walk; 
  int *start; 
  int *end; 
  int *midpoint; 

  PointsMatrixType m_Points;
  CellsMatrixType m_Cells;
  
  BasicShapeType m_BasicShape;

  unsigned int m_Order;
  
  unsigned int m_NumberOfVertices;
  unsigned int m_NumberOfFaces;
  unsigned int m_NumberOfEdges;
  
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_SphereTessellator(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT SphereTessellator< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef SphereTessellator< ITK_TEMPLATE_2 x > SphereTessellator##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkSphereTessellator+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSphereTessellator_hxx)
# include "itkSphereTessellator.hxx"
#endif

#endif
