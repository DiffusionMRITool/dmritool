/**
 *       @file  itkSpamsWeightedLassoSolver.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "02-21-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkSpamsWeightedLassoSolver_h
#define __itkSpamsWeightedLassoSolver_h

#include "itkSolverBase.h"

#include <linalg.h>
#include <decomp.h>

namespace itk
{

/**
 *   \class   SpamsWeightedLassoSolver
 *   \brief   solve weighted LASSO using spams
 *
 *     For all columns b of B, and w of W, it computes one column x of X
 *     which is the solution of
 *
 *   0) when mode=0
 *   0) when mode=1
 *   2) when mode=2 (default)
 *   \f[
 *      \min_{x} \|Ax-b\|^2 + lambda* w^T abs(x)
 *   \f]
 *   where A is a matrix, w is a vector, b, x are vectors, and lambda is a scalar value.
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup OptimizationSolver
 */
template <class TPrecision>
class ITK_EXPORT SpamsWeightedLassoSolver 
  : public SolverBase<TPrecision>
{
public:
  /** Standard class typedefs. */

  typedef SpamsWeightedLassoSolver  Self;
  typedef SolverBase<TPrecision>                      Superclass;
  typedef SmartPointer<Self>                          Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SpamsWeightedLassoSolver, SolverBase);

  typedef typename Superclass::ValueType       ValueType;
  typedef typename Superclass::MatrixType      MatrixType;
  typedef typename Superclass::VectorType      VectorType;
  typedef typename Superclass::MatrixPointer   MatrixPointer;
  typedef typename Superclass::VectorPointer   VectorPointer;
  
  typedef typename spams::Matrix<double>       SpamsMatrixType;
  typedef typename spams::SpMatrix<double>     SpamsSpMatrixType;
  typedef typename spams::Vector<double>       SpamsVectorType;
  typedef typename utl_shared_ptr<spams::Matrix<double> >    SpamsMatrixPointer;
  typedef typename utl_shared_ptr<spams::SpMatrix<double> >  SpamsSpMatrixPointer;
  typedef typename utl_shared_ptr<spams::Vector<double> >    SpamsVectorPointer;

  typedef enum 
    {
    L1CONS=0, 
    L2CONS,
    PENALTY
    } ConstraintType;
  
  itkSetMacro(ConstraintType, ConstraintType);
  itkGetMacro(ConstraintType, ConstraintType);
  
  itkSetMacro(NumberOfThreads, int);
  itkGetMacro(NumberOfThreads, int);
  
  itkSetMacro(Lambda, double);
  itkGetMacro(Lambda, double);
  
  void SetA(const MatrixPointer& mat);
  itkGetMacro(A, MatrixPointer);
  void SetW(const MatrixPointer& W);
  void Setw(const VectorPointer& w);
  itkGetMacro(W, MatrixPointer);
  void SetB(const MatrixPointer& B);
  void Setb(const VectorPointer& b);
  itkGetMacro(B, MatrixPointer);
  
  itkGetMacro(X, MatrixPointer);
  
  itkSetMacro(Positive,bool);
  itkGetMacro(Positive,bool);
  itkBooleanMacro(Positive);

  
  int GetXDimension() const
    {
    int N = m_A->Columns();
    utlException(N==0, "wrong size! m_A.Columns()="<<m_A->Columns());
    return N;
    }
  int GetXNumber() const
    {
    int M = m_B->Columns();
    utlException(M==0, "wrong size! m_B.Columns()="<<m_B->Columns());
    return M;
    }

  void Clear() 
    { 
    Superclass::Clear(); 
    ClearA();
    ClearW();
    ClearB();
    }
  void ClearA() 
    { 
    m_A=MatrixPointer(new MatrixType()); m_As=SpamsMatrixPointer(new SpamsMatrixType());
    }
  void ClearW() 
    { 
    m_W=MatrixPointer(new MatrixType()); m_Ws=SpamsMatrixPointer(new SpamsMatrixType());
    }
  void ClearB() 
    { 
    m_B=MatrixPointer(new MatrixType()); m_Bs=SpamsMatrixPointer(new SpamsMatrixType());
    }
  

  void VerifyInputs() const;
  
  void Solve(); 
  
  ValueType EvaluateCostFunction(const VectorType& x, const int col) const;
  ValueType EvaluateCostFunction(const MatrixType& x=MatrixType()) const;
  

protected:
  SpamsWeightedLassoSolver();
  virtual ~SpamsWeightedLassoSolver() {}
  
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual typename LightObject::Pointer InternalClone() const;
  
  /** MxN matrix  */
  MatrixPointer m_A;
  MatrixPointer m_B;
  MatrixPointer m_W;
  MatrixPointer m_X;
  
  SpamsMatrixPointer m_As;
  SpamsMatrixPointer m_Bs;
  SpamsMatrixPointer m_Ws;
  // SpamsMatrixType m_Xs;
  SpamsSpMatrixPointer m_Xs;

  ConstraintType m_ConstraintType;

  bool m_Positive;

  double m_Lambda;
  int m_NumberOfThreads;

private:
  SpamsWeightedLassoSolver(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_SpamsWeightedLassoSolver(_, EXPORT, TypeX, TypeY)         \
  namespace itk                                                 \
  {                                                             \
  _( 1 ( class EXPORT SpamsWeightedLassoSolver< ITK_TEMPLATE_1 TypeX > ) )     \
  namespace Templates                                           \
  {                                                             \
  typedef SpamsWeightedLassoSolver< ITK_TEMPLATE_1 TypeX > SpamsWeightedLassoSolver##TypeY; \
  }                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkSpamsWeightedLassoSolver+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkSpamsWeightedLassoSolver_hxx)
#include "itkSpamsWeightedLassoSolver.hxx"
#endif

#endif 
