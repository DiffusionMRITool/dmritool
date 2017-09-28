/**
 *       @file  itkDiffusionTensor.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "11-01-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */


#ifndef __itkDiffusionTensor_h
#define __itkDiffusionTensor_h

#include <itkDiffusionTensor3D.h>
#include "utlNDArray.h"

namespace itk
{

/**
 *   \class  DiffusionTensor
 *   \brief  tensor with some useful functions
 *   
 *   The elements is stored as a 6 dimentional array.
 *   \f[
 *    \begin{array}{ccc}
 *        0 & 1 & 2   \\
 *        X & 3 & 4   \\
 *        X & X & 5   
 *    \end{array} 
 *   \f]
 *
 *
 *   \ingroup DiffusionModels
 *   \author Jian Cheng
 */
template< class TPrecision >
class DiffusionTensor : public DiffusionTensor3D<TPrecision>
{

public:
  /** Standard class typedefs. */
  typedef DiffusionTensor      Self;
  typedef DiffusionTensor3D<TPrecision> Superclass;

  const static unsigned int NDimension = 3;
  itkStaticConstMacro (NDegreesOfFreedom, unsigned int, NDimension*(NDimension+1)/2);

  typedef typename Superclass::ValueType     ValueType;
  typedef typename Superclass::ComponentType ComponentType;
  
  typedef typename Superclass::EigenValuesArrayType   EigenValuesArrayType;
  typedef typename Superclass::EigenVectorsMatrixType EigenVectorsMatrixType;
  
  typedef typename Superclass::AccumulateValueType AccumulateValueType;
  typedef typename Superclass::RealValueType       RealValueType;
  
  typedef typename Superclass::SymmetricEigenAnalysisType  SymmetricEigenAnalysisType;
    
  typedef Vector<TPrecision, NDimension>            VectorType;

  /** Defaut constructor */
  inline DiffusionTensor(): Superclass(){};

  /** copy constructors. */
  inline DiffusionTensor (const Self& r): Superclass (r){}

  template< class TTensorValueType > 
    inline DiffusionTensor(const DiffusionTensor< TTensorValueType>& v): Superclass(v) {}
  
  inline bool operator==(const Self & mat) const
    {
    for ( int i = 0; i < NDegreesOfFreedom; ++i ) 
      {
      if ((*this)[i]!=mat[i])
        return false;
      }
    return true;
    }
  inline bool operator!=(const Self& mat) const
    {
    return !operator==(mat);
    }
  
  inline bool IsSame(const Self & mat, const double eps=1e-8) const
    {
    for ( int i = 0; i < NDegreesOfFreedom; ++i ) 
      {
      if ( std::fabs((*this)[i]-mat[i]) > eps )
        return false;
      }
    return true;
    }
  
  template< class TMatrixType > 
  inline Self& operator=(const TMatrixType & mat)
    {
    for ( unsigned int i = 0; i < NDimension; i++ )
      for(unsigned int j=i; j<NDimension; j++)
        (*this)[NDimension*i-(i+1)*i/2+j] = mat(i,j);
    return *this;
    }

  inline Self operator-(const Self & mat) const
    {
    Self result;
    for ( unsigned int i = 0; i < NDimension; i++ )
      result[i] = ( *this )[i] - mat[i];
    return result;
    }
  inline Self operator+(const Self & mat) const
    {
    Self result;
    for ( unsigned int i = 0; i < NDimension; i++ )
      result[i] = ( *this )[i] + mat[i];
    return result;
    }

  inline Self operator*(const Self& A) const
    {
    Self res;
    res.SetVnlMatrix ( this->GetVnlMatrix() * A.GetVnlMatrix() );
    return res;
    }
  inline vnl_matrix<TPrecision> operator*(const vnl_matrix<TPrecision>& A) const
    {
    return this->GetVnlMatrix() * A;
    }
    
  /** Tensor multiplication by a vector - returns a vector*/
  inline VectorType operator*(const VectorType& v) const
    {
    VectorType result;
    for( unsigned int i=0; i<NDimension; i++ )
      {
      result[i]=0;
      for( unsigned int j=0; j<NDimension; j++)
        result[i] += (*this)(i,j)*v[j];
      }
    return result;
    }

  /** access matrix values, which is similar with the ones in \ref DiffusionTensor3D, 
   * but with boundary check. */
  ValueType & operator()(unsigned int row, unsigned int col);

  const ValueType & operator()(unsigned int row, unsigned int col) const;

  void Flip(const int flipx, const int flipy, const int flipz)
    {
    double sign_xy = flipx^flipy ? -1 : 1;
    double sign_xz = flipx^flipz ? -1 : 1;
    double sign_yz = flipy^flipz ? -1 : 1;
    (*this)[1] = sign_xy*(*this)[1];
    (*this)[2] = sign_xz*(*this)[2];
    (*this)[4] = sign_yz*(*this)[4];
    }

  /** 
   * \brief Return an array containing EigenValues in ascending order, and a matrix containing the corresponding Eigenvectors. 
   *
   * \param eigenVectors eigenvectors stored in vnl_matrix<TPrecision>, where each column is an eigenvector.
   * \note eigenvectors in \ref itk::SymmetricSecondRankTensor::ComputeEigenAnalysis stored as rows of \ref itk::Matrix
   * */
  void GetEigenValuesVectors(vnl_diag_matrix<TPrecision> & eigenValues, vnl_matrix<TPrecision> & eigenVectors) const;
  
  /** 
   *  \brief analytic way to calculate eigenValues (in ascending order) and eigenVectors.  
   *  In mathematica, Eigenvectors[( { {a, b, c}, {b, d, e}, {c, e, f}  } )] 
   *
   * \param eigenVectors eigenvectors stored in vnl_matrix<TPrecision>, where each column is an eigenvector.
   *
   * \note: Eigenvalues can be obtained analytically. 
   * But there is a numerical issue when dividing a small number for calculating eigenvectors. 
   * In this case, use a numerical way for eigenvectors.
   * */
  template<class TArrayType, class TMatrixType >
  void GetEigenValuesVectorsAnalytic(TArrayType& eigenValues, TMatrixType& eigenVectors) const;

  template <class TArrayType >
  void SetEigenValues(const TArrayType& array);
  template <class TMatrixType >
  void SetEigenVectors(const TMatrixType& vectors);

  /** Set a vnl_matrix_ref referencing the same memory block. */
  inline void SetVnlMatrix( const vnl_matrix<TPrecision> & mat);

  /** Get a vnl_matrix with a copy of the internal memory  block. */
  inline vnl_matrix<TPrecision> GetVnlMatrix () const;

  inline double GetDeterminant() const;
  inline double GetNorm() const;

  /** Transform the tensor with a matrix */
  inline Self GetRotate (const typename Self::MatrixType& mat) const;
  inline Self& Rotate (const typename Self::MatrixType& mat);

  /** Transform the tensor with a matrix */
  inline vnl_matrix<TPrecision> GetRotate (const vnl_matrix<TPrecision> & vec) const;
  
  double GetQuadraticForm( const double x, const double y, const double z) const;
  
  /** get spherical samples based on given gradients (Cartesian form).
   * \note No boundary checking. The size of the input vec should be the same as the gradients.rows(), and vec should have the operator[].
   * */
  template <class TVectorType >
  void GetSphericalSamples( TVectorType& vec, const utl::NDArray<TPrecision,2>& gradients) const;

  /** get the DWI samples based on given gradients (Cartesian form) and b values.
   * \note No boundary checking. The size of the input vec should be the same as the gradients.rows(), and vec should have the operator[].
   * */
  template <class TVectorType >
  void GetDWISamples( TVectorType& vec, const utl::NDArray<TPrecision,2>& gradients, const std::vector<TPrecision>& bValues) const;
  
  /** get the EAP samples based on given gradients (Cartesian form) and r values.
   * \note No boundary checking. The size of the input vec should be the same as the gradients.rows(), and vec should have the operator[].
   * */
  template <class TVectorType >
  void GetEAPSamples( TVectorType& vec, const utl::NDArray<TPrecision,2>& gradients, const std::vector<TPrecision>& rValues, const TPrecision& tau) const;

  /** get the ODF samples based on given gradients (Cartesian form) and odf order. 
   * When gradients are evenly distributed in \f$S^2\f$, the odf samples can be numerically normalized if isNormalize=true. 
   * \note No boundary checking. The size of the input vec should be the same as the gradients.rows(), and vec should have the operator[].
   * */
  template <class TVectorType >
  void GetODFSamples( TVectorType& vec, const utl::NDArray<TPrecision,2>& gradients, const int& odfOrder=2, const bool& isNormalize=false) const;

  /** Return-To-Origin probability  */
  RealValueType GetReturnToOrigin (const TPrecision tau) const;
  /** Mean Squared Displacement  */
  RealValueType GetMeanSquaredDisplacement (const TPrecision tau) const;
  /** GA  */
  RealValueType GetGA() const;
  /** FA  */
  RealValueType GetFA() const;
  /** Mode  */
  RealValueType GetMODE() const;
  /** RA  */
  RealValueType GetRA() const;
  /** Mean Diffusivity  */
  RealValueType GetADC() const;
  RealValueType GetMD() const 
    { return GetADC(); }

  bool IsFinite() const;
  bool HasNans() const;
  bool IsZero(const double eps=1e-10) const;

  bool IsPositive() const;
  bool IsDiagonal(const double eps=1e-10) const;
    
  void Positivize() const;
  
  /** Return the matrix logarithm */
  inline Self Log() const;

  /** Return the matrix exponential */
  inline Self Exp () const;

  /** Power of the tensor */
  inline Self Pow (const double& n) const;
  inline Self Inv () const;
  inline Self Sqrt () const;
  inline Self InvSqrt () const;

  inline TPrecision EuclideanDistance(const DiffusionTensor<TPrecision>& X) const;
  inline TPrecision KLDistance(const DiffusionTensor<TPrecision>& X) const;
  inline TPrecision GeodesicDistance(const DiffusionTensor<TPrecision>& X) const;
  inline TPrecision LogEucDistance(const DiffusionTensor<TPrecision>& X) const;

  void Print(std::ostream& os, Indent indent=0) const
    {
    os << indent << "[";
    for ( int i = 0; i < NDimension; i += 1 ) 
      {
      for ( int j = 0; j < NDimension -1; j += 1 ) 
        os<< (*this)(i,j) << " ";  
      os << (*this)(i,NDimension-1) << "; ";
      }
    os <<"];" << std::endl;
    }

};


}

#define ITK_TEMPLATE_DiffusionTensor(_, EXPORT, TypeX, TypeY)     \
  namespace itk                                                             \
  {                                                                         \
  _( 2 ( class EXPORT DiffusionTensor< ITK_TEMPLATE_2 TypeX > ) ) \
  namespace Templates                                                       \
  {                                                                         \
  typedef DiffusionTensor< ITK_TEMPLATE_2 TypeX >                 \
  DiffusionTensor##TypeY;                                       \
  }                                                                         \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkDiffusionTensor+-.h"
#endif

#if !defined(ITK_MANUAL_INSTANTIATION) && !defined(__itkDiffusionTensor_hxx)
#include "itkDiffusionTensor.hxx"
#endif


#endif 


