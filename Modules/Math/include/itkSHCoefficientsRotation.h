/**
 *       @file  itkSHCoefficientsRotation.h
 *      @brief  
 *     Created  "03-26-2014
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */

#ifndef __itkSHCoefficientsRotation_h
#define __itkSHCoefficientsRotation_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "utlNDArray.h"
#include "utlITKMacro.h"
#include "utl.h"

namespace itk
{


/**
 *   \class   SHCoefficientsRotation
 *   \brief   rotate SH coefficient vector by a given rotation matrix.
 *
 *  Reference: "Efficient and accurate rotation of finite spherical harmonics expansions", 
 *  Journal of Computational Physics 231 (2012) 243–250 
 *
 *   \ingroup Math
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class T = double  >
class ITK_EXPORT SHCoefficientsRotation
  : public Object
{
public:
  /** Standard class typedefs. */
  typedef SHCoefficientsRotation         Self;
  typedef Object  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( SHCoefficientsRotation, Object );
  
  typedef T                                          PreciseType; 

  typedef utl::NDArray<double,1>                     VectorType;
  typedef utl::NDArray<double,2>                     MatrixType;
  typedef utl_shared_ptr<MatrixType>                 MatrixPointer;

  itkSetGetMacro(TessOrder, int);
  itkSetGetMacro(MaxRank, int);

  void Initialize()
    {
    utl::GradientTable<double>::Initialize(m_TessOrder);
    if (m_MaxRank/2>m_SHMatrixInverse->size())
      {
      m_Grad = utl::GradientTable<double>::GetGrad(m_TessOrder, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
      MatrixPointer shMatrix = utl::ComputeSHMatrix(m_MaxRank, *m_Grad, CARTESIAN_TO_SPHERICAL);
      m_SHMatrixInverse = utl_shared_ptr<std::vector<MatrixType> >(new std::vector<MatrixType>(m_MaxRank/2));
      for ( int l = 2; l <= m_MaxRank; l += 2 ) 
        {
        MatrixType shMatrixL;
        shMatrixL = shMatrix->GetNColumns(utl::RankToDimSH(l-2), 2*l+1);
        (*m_SHMatrixInverse)[l/2-1] = utl::PInverseMatrix(shMatrixL, 1e-15);
        }
      }
    }
  
/**
* @brief GetRotatedSHCoefficient 
*        rotate SH coefficient vector by a given rotation matrix
*
* @param shInput input SH coefficients
* @param rotationMatrix rotation matrix
* @author Jian Cheng
*
*  Reference: "Efficient and accurate rotation of finite spherical harmonics expansions", 
*  Journal of Computational Physics 231 (2012) 243–250 
*
*  use itk::SHCoefficientsRotation for multi-thread programming
*/
  VectorType GetRotatedSHCoefficients(const VectorType& shInput, const MatrixType& rotationMatrix) const
    {
    utlSAException(m_SHMatrixInverse->size()==0).msg("need to use Initialize() first");
    int rankReal = utl::DimToRankSH(shInput.Size());
    utlSAException(rankReal>m_MaxRank)(rankReal)(m_MaxRank).msg("rank is larger than m_MaxRank.");
    MatrixType gradt3Rotated =  (*m_Grad) * rotationMatrix;
    MatrixPointer shMatrixRotated = utl::ComputeSHMatrix(rankReal, gradt3Rotated,CARTESIAN_TO_SPHERICAL);
    VectorType shRotated(shInput); 
    MatrixType tmp;

    for ( int l = 2; l <= rankReal; l += 2 ) 
      {
      int colstart = utl::RankToDimSH(l-2);
      tmp = shMatrixRotated->GetNColumns(colstart, 2*l+1);
      VectorType sfValueRotatedL =  tmp* shInput.GetSubVector(utl::GetRange(colstart, colstart+2*l+1));
      VectorType shRotatedL = (*m_SHMatrixInverse)[l/2-1]*sfValueRotatedL;
      shRotated.SetSubVector(utl::GetRange(colstart, colstart+shRotatedL.Size()), shRotatedL);
      }
    return shRotated;
    }

  static Pointer GetInstance()
    {
    if (!m_Instance)
      m_Instance = Self::New();
    return m_Instance;
    }
  

protected:
  SHCoefficientsRotation() : Superclass(),
  m_Grad(new MatrixType()), 
  m_SHMatrixInverse(new std::vector<MatrixType>())
    {
    m_MaxRank = 10;
    m_TessOrder = 3;
    }
  virtual ~SHCoefficientsRotation() {};

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    PrintVar2(true, m_MaxRank, m_TessOrder, os <<indent);
    utl::PrintUtlMatrix(*m_Grad, "m_Grad", " ", os << indent);
    for ( int i = 0; i < m_SHMatrixInverse->size(); ++i ) 
      {
      PrintVar(true, os << indent, i);
      utl::PrintUtlMatrix((*m_SHMatrixInverse)[i], "(*m_SHMatrixInverse)[i]", " ", os << indent);
      }
    }
  
  typename LightObject::Pointer InternalClone() const
    {
    typename LightObject::Pointer loPtr = Superclass::InternalClone();
    typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
    if(rval.IsNull())
      {
      itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
      }
    rval->m_MaxRank = m_MaxRank;
    rval->m_TessOrder = m_TessOrder;
    rval->m_Grad = m_Grad;
    rval->m_SHMatrixInverse = m_SHMatrixInverse;
    return loPtr;
    }

  static Pointer m_Instance;
  
  MatrixPointer m_Grad;
  utl_shared_ptr<std::vector<MatrixType> > m_SHMatrixInverse;

  /** maximal rank */
  int m_MaxRank;
  /** tessellation order */
  int m_TessOrder;

private:
  SHCoefficientsRotation(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented

};

template<>
typename SHCoefficientsRotation<double>::Pointer SHCoefficientsRotation<double>::m_Instance = ITK_NULLPTR;  

inline void 
InitializeSHRotation(const int rank, const int tess)
{
  typedef itk::SHCoefficientsRotation<double> SHRotationFilterType;
  SHRotationFilterType::Pointer shRotate = SHRotationFilterType::GetInstance();  
  shRotate->SetTessOrder(tess); 
  shRotate->SetMaxRank(rank);     
  shRotate->Initialize();    
}

}


#endif 
