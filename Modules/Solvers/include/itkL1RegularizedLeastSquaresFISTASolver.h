/**
 *       @file  itkL1RegularizedLeastSquaresFISTASolver.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "05-03-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkL1RegularizedLeastSquaresFistaSolver_h
#define __itkL1RegularizedLeastSquaresFistaSolver_h

#include "itkIterativeSolverBase.h"
#include "itkL2RegularizedLeastSquaresSolver.h"

namespace itk
{

/**
 *   \class   L1RegularizedLeastSquaresFISTASolver
 *   \brief   solve least square problem with L1 regularization using FISTA
 *
 *   The least square with L2 regularization is 
 *   \f[
 *      \min_{x} \|Ax-b\|^2 + |w x|
 *   \f]
 *   where A is a matrix, w is a diagonal matrix (vnl_vector)
 *   b, x are vectors. 
 *
 *  reference: Fast Iterative Shrinkage-Thresholding Algorithm (FISTA), SIAM J. Imaging Sciences 2009
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup OptimizationSolver
 */
template <class TPrecision>
class ITK_EXPORT L1RegularizedLeastSquaresFISTASolver : public IterativeSolverBase<TPrecision>
{
public:
  /** Standard class typedefs. */

  typedef L1RegularizedLeastSquaresFISTASolver  Self;
  typedef IterativeSolverBase<TPrecision>                      Superclass;
  typedef SmartPointer<Self>                          Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(L1RegularizedLeastSquaresFISTASolver, IterativeSolverBase);

  typedef typename Superclass::ValueType ValueType;
  typedef typename Superclass::MatrixType      MatrixType;
  typedef typename Superclass::VectorType      VectorType;
  typedef typename Superclass::MatrixPointer   MatrixPointer;
  typedef typename Superclass::VectorPointer   VectorPointer;

  typedef typename Superclass::ValueContainerType ValueContainerType;
  typedef typename Superclass::UpdateInfomationType UpdateInfomationType;
  typedef L2RegularizedLeastSquaresSolver<TPrecision>  L2SolverType;
  
  void SetA(const MatrixPointer& mat);
  itkGetMacro(A, MatrixPointer);
  void Setw(const VectorPointer& w);
  void SetwForInitialization(const VectorPointer& w);
  itkGetMacro(w, VectorPointer);
  void Setb(const VectorPointer& b);
  itkGetMacro(b, VectorPointer);
  
  itkSetMacro(UseL2SolverForInitialization,bool);
  itkGetMacro(UseL2SolverForInitialization,bool);
  itkBooleanMacro(UseL2SolverForInitialization);
  
  int GetXDimension() const
    {
    int N = m_A->Columns();
    utlException(N==0, "wrong size! m_A->Columns()="<<m_A->Columns());
    return N;
    }

  void Clear() 
    { 
    Superclass::Clear(); 
    m_A=MatrixPointer(new MatrixType()); 
    m_AtA=MatrixPointer(new MatrixType()); 
    m_Atb=VectorPointer(new VectorType()); 
    m_b=VectorPointer(new VectorType()); 
    m_w=VectorPointer(new VectorType()); 
    m_L2Solver->Clear();
    }
  void ClearA() 
    { 
    m_A=MatrixPointer(new MatrixType()); 
    m_AtA=MatrixPointer(new MatrixType()); 
    m_Atb=VectorPointer(new VectorType()); 
    m_L2Solver->ClearA();
    }
  void Clearw() 
    { 
    m_w=VectorPointer(new VectorType()); 
    m_L2Solver->ClearLamabda();
    }
  void Clearb() 
    { 
    m_b=VectorPointer(new VectorType()); 
    m_Atb=VectorPointer(new VectorType()); 
    m_L2Solver->Clearb();
    }
  
  /** Update history information and monitor stop conditions  */
  void HistoryUpdateAndConvergenceCheck();

  void VerifyInputs() const;
  
  void Solve(const VectorType& xInitial=VectorType()); 
  void Iterate(); 
  void Initialize(const VectorType& xInitial=VectorType()); 
  
  ValueType EvaluateCostFunction(const VectorType& x=VectorType()) const;
  

protected:
  L1RegularizedLeastSquaresFISTASolver();
  virtual ~L1RegularizedLeastSquaresFISTASolver() {}
  
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual typename LightObject::Pointer InternalClone() const;
  
  /** MxN matrix  */
  MatrixPointer m_A;
  /** Mx1 vector  */
  VectorPointer m_b;
  
  /** N dimensional vector for L1 regularization. */
  VectorPointer m_w;

  bool m_UseL2SolverForInitialization;

private:
  L1RegularizedLeastSquaresFISTASolver(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MatrixPointer m_At;
  MatrixPointer m_AtA;
  VectorPointer m_Atb;
  double m_Step;
  
  // private members
  VectorPointer m_xOld;

  typename L2SolverType::Pointer m_L2Solver;
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_L1RegularizedLeastSquaresFISTASolver(_, EXPORT, TypeX, TypeY)         \
  namespace itk                                                 \
  {                                                             \
  _( 1 ( class EXPORT L1RegularizedLeastSquaresFISTASolver< ITK_TEMPLATE_1 TypeX > ) )     \
  namespace Templates                                           \
  {                                                             \
  typedef L1RegularizedLeastSquaresFISTASolver< ITK_TEMPLATE_1 TypeX > L1RegularizedLeastSquaresFISTASolver##TypeY; \
  }                                                             \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkL1RegularizedLeastSquaresFISTASolver+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkL1RegularizedLeastSquaresFISTASolver_hxx)
#include "itkL1RegularizedLeastSquaresFISTASolver.hxx"
#endif

#endif 

