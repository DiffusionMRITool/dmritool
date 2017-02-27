/**
 *       @file  itkSHCoefficientsFit.h
 *      @brief  
 *     Created  "10-06-2016
 *
 *     @author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 * =====================================================================================
 */


#ifndef __itkSHCoefficientsFit_h
#define __itkSHCoefficientsFit_h

#include "utlITKMacro.h"
#include "utl.h"


namespace itk
{

namespace Functor
{

/**
 *   \class   SHCoefficientsFit
 *   \brief   fit SH coefficients from spherical function samples.
 *
 *
 *   \ingroup Math
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 */
template < class T = double  >
class ITK_EXPORT SHCoefficientsFit
  : public utl::Functor::VectorFunctorBase<utl::Vector<T> >
{
public:
  /** Standard class typedefs. */
  typedef SHCoefficientsFit         Self;
  typedef utl::Functor::VectorFunctorBase<utl::Vector<T> >  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  
  typedef T                                          PreciseType; 

  typedef utl::NDArray<double,1>                     VectorType;
  typedef utl::NDArray<double,2>                     MatrixType;
  typedef utl_shared_ptr<MatrixType>                 MatrixPointer;
  
  SHCoefficientsFit() : Superclass(),
  m_Orientations(new MatrixType()),
  m_BasisSHMatrix(new MatrixType()),
  m_BasisSHMatrixPInv(new MatrixType())
    {
    }
  virtual ~SHCoefficientsFit() {};

  bool operator==(const Self & other) const  
  { return false; }
  bool operator!=(const Self & other) const 
  { return ! (*this==other); } 

  
  SHCoefficientsFit(const Self& other)
    {
    *this = other;
    }
  
  Self& operator=(const Self& other)
    {
    Superclass::operator=(other);
    if ( this != &other ) 
      {
      m_SHRank = other.m_SHRank;
      m_Lambda = other.m_Lambda;
      m_Power = other.m_Power;
      m_Orientations = other.m_Orientations;
      m_BasisSHMatrix = other.m_BasisSHMatrix;
      m_BasisSHMatrixPInv = other.m_BasisSHMatrixPInv;
      }
    return *this;
    }

  utlSetGetMacro(SHRank, int);
  utlSetGetMacro(Lambda, double);
  utlSetGetMacro(Power, double);
  
  utlSetGetMacro(Orientations, MatrixPointer);

  utlGetMacro(BasisSHMatrix, MatrixPointer);
  utlGetMacro(BasisSHMatrixPInv, MatrixPointer);
  
  void VerifyInputParameters(const int inputVecSize=-1) const
    {
    utlVLogPosition(LOG_DEBUG);
    Superclass::VerifyInputParameters();
    utlGlobalException(m_SHRank<0, "need to set m_SHRank");
    utlSAGlobalException(m_Orientations->Rows()==0).msg("need to set m_Orientations");

    if (inputVecSize>0)
      {
      utlSAGlobalException(m_Orientations->Rows()!=inputVecSize)
        (m_Orientations->Rows())(inputVecSize).msg("need to set orientations correctly");
      }
    }

  void Initialize()
    {
    utlVLogPosition(LOG_DEBUG);
    VerifyInputParameters();

    m_BasisSHMatrix = utl::ComputeSHMatrix(m_SHRank, *m_Orientations, CARTESIAN_TO_SPHERICAL);
    int J = m_BasisSHMatrix->Cols();
    int N = m_BasisSHMatrix->Rows();

    MatrixType shT = m_BasisSHMatrix->GetTranspose();

    if (m_Lambda>0)
      {
      MatrixType regMat(J, J, 0.0);
      for ( int i = 0; i < J; ++i ) 
        {
        std::vector<int> lm = utl::GetIndexSHlm(i);
        regMat(i,i) = m_Lambda*lm[0]*lm[0]*(lm[0]+1.)*(lm[0]+1.);
        }
      if (utl::IsLogDebug(this->m_LogLevel))
        utl::PrintUtlMatrix(regMat, "regMat");
      *m_BasisSHMatrixPInv = utl::InverseSymmericMatrix( utl::Eval<double,2>(shT * (*m_BasisSHMatrix) + regMat) )* shT;
      }
    else
      *m_BasisSHMatrixPInv = utl::InverseSymmericMatrix( shT * (*m_BasisSHMatrix) )* shT;

    if (utl::IsLogDebug(this->m_LogLevel))
      {
      utl::PrintUtlMatrix(*m_BasisSHMatrix, "m_BasisSHMatrix");
      utl::PrintUtlMatrix(*m_BasisSHMatrixPInv, "m_BasisSHMatrixPInv");
      }
    }
  
  VectorType operator()(const VectorType& sf) const
    {
    if (std::fabs(m_Power-1.0)>1e-8)
      {
      VectorType sf_power = utl::Pow(sf, m_Power);
      return (*m_BasisSHMatrixPInv) * sf_power;
      }
    else
      return (*m_BasisSHMatrixPInv) * sf;
    }

  int GetOutputDimension(const int ) const
    {
    utlSAException(m_SHRank<0)(m_SHRank).msg("need to set m_SHRank");
    return utl::RankToDimSH(m_SHRank);
    }
  
  void Print(std::ostream & os=std::cout) const
    {
    Superclass::Print(os);
    utlLogOSVar(os , m_Power, m_SHRank, m_Lambda);
    utl::PrintUtlMatrix(*m_Orientations, "m_Orientations", " ", os);
    utl::PrintUtlMatrix(*m_BasisSHMatrix, "m_BasisSHMatrix", " ", os);
    utl::PrintUtlMatrix(*m_BasisSHMatrixPInv, "m_BasisSHMatrixPInv", " ", os);
    }

protected:

  

  /** rank for SH basis  */
  int m_SHRank=-1;

  /** L2 regularization parameter.  */
  double m_Lambda=0;

  double m_Power=1.0;
  
  /** Orientation matrix in Cartesian format  */
  MatrixPointer  m_Orientations;
  
  /** SH basis matrix  */
  MatrixPointer  m_BasisSHMatrix;
  /** pseudoinverse of m_BasisSHMatrix, considering L2 regularization using m_Lambda  */
  MatrixPointer  m_BasisSHMatrixPInv;


private:

};

}

}


#endif 
