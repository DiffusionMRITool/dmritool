/**
 *       @file  itkMeshFromFiberTractsFilter.h
 *      @brief  
 *     Created  "08-23-2017
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkMeshFromFiberTractsFilter_h
#define __itkMeshFromFiberTractsFilter_h

#include <itkLightProcessObject.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "itkFiberTracts.h"
#include "utlCoreMacro.h"

namespace itk
{

/** \class MeshFromFiberTractsFilter
 * \brief Compute mesh from fibers. 
 *
 * \ingroup Visualization
 * \author  Jian Cheng (jian.cheng.1983@gmail.com)
 *
 */
class ITK_EXPORT MeshFromFiberTractsFilter : public LightProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MeshFromFiberTractsFilter       Self;
  typedef LightProcessObject   Superclass;
  typedef SmartPointer<Self>   Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshFromFiberTractsFilter, LightProcessObject);

  typedef FiberTracts<double>                     FiberTractsType;
  typedef typename FiberTractsType::Pointer       FiberTractsPointer;
  typedef typename FiberTractsType::FiberType     FiberType;
  typedef typename FiberTractsType::FiberPointer  FiberPointer;
  typedef typename FiberType::STDVectorType       STDVectorType;
  typedef typename FiberType::VertexType          VertexType;

  typedef vtkPolyData OutputMeshPolyDataType;
  typedef vtkSmartPointer<vtkPolyData> OutputMeshPolyDataPointer;
  
  /** enum of the possible shapes for the glyphs */
  typedef enum 
  {
  /** lines  */
    GLYPH_LINE=0,      //0
    /** tubes   */
    GLYPH_TUBE 
  } GlyphShapeType;
  
  typedef enum
  {
  /** default color  */
    COLOR_NONE=0, 
  /** fixed color  */
    COLOR_FIXED, 
  /** by direction in each point */
    COLOR_BY_POINT_DIRECTION, 
  /** by mean direction in a tract  */
    COLOR_BY_MEAN_DIRECTION,   
  /** by start and end points  */
    COLOR_BY_ENDPOINTS_DIRECTION,   
  /** by an input image  */
    COLOR_BY_IMAGE,
  /** by scalar  */
    COLOR_BY_SCALARS,
  /** by scalar  */
    COLOR_BY_PROPERTY
  } ColorSchemeType;

  itkSetGetMacro(FiberTracts, FiberTractsPointer);
  
  itkSetGetMacro(ColorScheme, ColorSchemeType);
  itkSetGetMacro(ShapeMode, GlyphShapeType);
  
  itkSetGetMacro(TubeRadius, double);
  
  void SetColor(double r, double g, double b)
    {
    m_Color[0]=r, m_Color[1]=g, m_Color[2]=b;
    }
  const double* GetColor() const
    {
    return m_Color;
    }
  
  void SetFlip(const int flipx, const int flipy, const int flipz )
    {
    m_Flip[0]=flipx, m_Flip[1]=flipy, m_Flip[2]=flipz;
    }
  std::vector<int> GetFlip() const
    {
    return m_Flip;
    }

  /** Does the real work. */
  virtual void Update();

  OutputMeshPolyDataPointer GetOutput() const
    {
    return this->m_Mesh;
    }

    
protected:
  MeshFromFiberTractsFilter(){}
  ~MeshFromFiberTractsFilter(){}
  
  typename LightObject::Pointer InternalClone() const;

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    m_FiberTracts->Print(os, indent);
    utl::PrintVector(m_Flip, "m_Flip");
    PrintVar(true, os<<indent, m_TubeRadius, m_ShapeMode, m_ColorScheme);
    PrintVar(true, os<<indent, m_Color[0], m_Color[1], m_Color[2]);
    m_Mesh->Print(std::cout<<"m_Mesh=");
    }
  
  FiberTractsPointer m_FiberTracts = FiberTractsType::New();

  ColorSchemeType m_ColorScheme=COLOR_BY_POINT_DIRECTION;

  GlyphShapeType m_ShapeMode=GLYPH_LINE;
  
  OutputMeshPolyDataPointer m_Mesh = OutputMeshPolyDataType::New();
  
  double m_Color[3]={255,0,0};
  
  double m_TubeRadius=0.2;

  /** flips in x/y/z-axis. 0 means no flip, 1 means flip  */
  std::vector<int> m_Flip = std::vector<int>(3, 0);

private:
  MeshFromFiberTractsFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshFromFiberTractsFilter.hxx"
#endif


#endif 
