/**
 *       @file  itkL2RegularizedLeastSquaresSolver.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-25-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkL2RegularizedLeastSquaresSolver_h
#define __itkL2RegularizedLeastSquaresSolver_h

#include "vnl/vnl_matrix.h"
#include "itkSolverBase.h"


namespace itk
{

/**
 *   \class   L2RegularizedLeastSquaresSolver
 *   \brief   solve least square problem with L2 regularization
 *
 *   The least square with L2 regularization is 
 *   \f[
 *      \min_{x} \|Ax-b\|^2 + x^T \Lambda x
 *   \f]
 *   where A, \Lambda are matrices, 
 *   b, x are vectors. 
 *
 *   \author  Jian Cheng (JC), jian.cheng.1983@gmail.com
 *   \ingroup OptimizationSolver
 */
template <class TPrecision>
class ITK_EXPORT L2RegularizedLeastSquaresSolver : public SolverBase<TPrecision>
{
public:
  /** Standard class typedefs. */
  typedef L2RegularizedLeastSquaresSolver         Self;
  typedef SolverBase<TPrecision>                  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( L2RegularizedLeastSquaresSolver, SolverBase );

  typedef typename Superclass::ValueType       ValueType;
  typedef typename Superclass::MatrixType      MatrixType;
  typedef typename Superclass::VectorType      VectorType;
  typedef typename Superclass::MatrixPointer   MatrixPointer;
  typedef typename Superclass::VectorPointer   VectorPointer;

  typedef typename Superclass::ValueContainerType ValueContainerType;

  /** m_LS is released when setting m_A and m_lambda  */
  void SetA(const MatrixPointer& mat);
  itkGetMacro(A, MatrixPointer);
  void SetLambda(const MatrixPointer& mat);
  itkGetMacro(Lambda, MatrixPointer);
  itkSetMacro(b, VectorPointer);
  itkGetMacro(b, VectorPointer);
  
  itkGetMacro(LS, MatrixPointer);
  itkGetMacro(ConditionNumber, ValueType);

  int GetXDimension() const ITK_OVERRIDE
    {
    int N = m_A->Columns();
    utlException(N==0, "wrong size! m_A->Columns()="<<m_A->Columns());
    return N;
    }

  /** m_LS is released when releasing m_A and m_lambda  */
  void Clear() ITK_OVERRIDE;
  void ClearA() 
    { 
    m_A=MatrixPointer(new MatrixType()); 
    m_LS=MatrixPointer(new MatrixType());
    m_ConditionNumber=-1;
    }
  void ClearLambda() 
    { 
    m_Lambda=MatrixPointer(new MatrixType()); 
    m_LS=MatrixPointer(new MatrixType()); 
    m_ConditionNumber=-1;
    m_IsLambdaSymmetric = true;
    }
  void Clearb() 
    { 
    m_b=VectorPointer(new VectorType()); 
    }

  void VerifyInputs() const ITK_OVERRIDE;
  
  void Solve(const VectorType& xInitial=VectorType()) ITK_OVERRIDE; 
  void Initialize(const VectorType& xInitial=VectorType()) ITK_OVERRIDE; 
  
  ValueType EvaluateCostFunction(const VectorType& x=VectorType()) const ITK_OVERRIDE;

protected:
  L2RegularizedLeastSquaresSolver();
  virtual ~L2RegularizedLeastSquaresSolver() {};

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
  
  virtual typename LightObject::Pointer InternalClone() const ITK_OVERRIDE;
  
  /** MxN matrix  */
  MatrixPointer m_A;
  /** Mx1 vector  */
  VectorPointer m_b;
  
  /** NxN matrix  */
  MatrixPointer m_Lambda;

  bool m_IsLambdaSymmetric;
  

private:
  L2RegularizedLeastSquaresSolver(const Self&);  //purposely not implemented
  void operator=(const Self&);  //purposely not implemented
  
  /** MxN matrix  */
  MatrixPointer m_LS;
  ValueType m_ConditionNumber;

};

}

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkL2RegularizedLeastSquaresSolver+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkL2RegularizedLeastSquaresSolver_hxx)
#include "itkL2RegularizedLeastSquaresSolver.hxx"
#endif


#endif 
