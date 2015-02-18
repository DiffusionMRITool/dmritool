/**
 *       @file  itkTensor.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  11-01-2012
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2012, Jian Cheng
 *
 * =====================================================================================
 */  

#ifndef __itkDiffusionTensor_hxx
#define __itkDiffusionTensor_hxx

#include "itkDiffusionTensor.h"
#include "utl.h"

#include "vnl/vnl_matrix.h"
#include <vnl/vnl_math.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_trace.h"
#include "vnl/algo/vnl_determinant.h"


namespace itk
{

/**
 * Matrix notation access to elements
 */
template< class TPrecision >
const typename DiffusionTensor< TPrecision >::ValueType &
DiffusionTensor< TPrecision >
::operator()(unsigned int row, unsigned int col) const
{
  unsigned int k;

  if ( row < col )
    {
    k = row * NDimension + col - row * ( row + 1 ) / 2;
    }
  else
    {
    k = col * NDimension + row - col * ( col + 1 ) / 2;
    }

  utlAssert(k < NDegreesOfFreedom, "wrong dimension");

  return ( *this )[k];
}

/**
 * Matrix notation access to elements
 */
template< class TPrecision >
typename DiffusionTensor< TPrecision >::ValueType &
DiffusionTensor< TPrecision >
::operator()(unsigned int row, unsigned int col)
{
  unsigned int k;

  if ( row < col )
    {
    k = row * NDimension + col - row * ( row + 1 ) / 2;
    }
  else
    {
    k = col * NDimension + row - col * ( col + 1 ) / 2;
    }
  utlAssert(k < NDegreesOfFreedom, "wrong dimension");

  return ( *this )[k];
}

template<class TPrecision >
void
DiffusionTensor<TPrecision>
::GetEigenValuesVectors(vnl_diag_matrix<TPrecision> & eigenValues, vnl_matrix<TPrecision> & eigenVectors) const
{
  EigenValuesArrayType values;
  EigenVectorsMatrixType vectors;
  this->ComputeEigenAnalysis(values, vectors);

  vnl_diag_matrix<TPrecision> vnl_eigenValues(NDimension);
  for ( int i = 0; i < NDimension; i += 1 ) 
    eigenValues[i] = values[i];
  eigenVectors = vectors.GetVnlMatrix().transpose();
}

template<class TPrecision>
template<class TArrayType >
void
DiffusionTensor<TPrecision>
::SetEigenValues(const TArrayType& array)
{
  if (this->IsDiagonal())
    {
    for ( int i = 0; i < NDimension; i += 1 ) 
      (*this)(i,i) = array[i];
    }
  else
    {
    vnl_diag_matrix<TPrecision> eigenValues(NDimension);
    vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
    this->GetEigenValuesVectors(eigenValues, eigenVectors);
    for ( int i = 0; i < NDimension; i += 1 ) 
      eigenValues[i] = array[i];
    *this = eigenVectors * eigenValues.asMatrix() * eigenVectors.transpose();
    }
}

template<class TPrecision>
template<class TMatrixType >
void
DiffusionTensor<TPrecision>
::SetEigenVectors(const TMatrixType& vectors)
{
  vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
  for ( int i = 0; i < NDimension; i += 1 ) 
    for ( int j = 0; j < NDimension; j += 1 ) 
      eigenVectors(i,j) = (TPrecision)vectors(i,j);
  if (this->IsDiagonal())
    *this = eigenVectors * (*this) * eigenVectors.transpose();
  else
    {
    vnl_diag_matrix<TPrecision> eigenValues(NDimension);
    this->GetEigenValuesVectors(eigenValues, eigenVectors);
    *this = eigenVectors * eigenValues.asMatrix() * eigenVectors.transpose();
    }
}


template<class TPrecision >
vnl_matrix<TPrecision>
DiffusionTensor<TPrecision>
::GetVnlMatrix(void) const
{
  TPrecision* block = new TPrecision[9];
  for(unsigned int i=0;i<NDimension;i++)
    for(unsigned int j=0;j<NDimension;j++)
      block[NDimension*i+j] = (*this)(i,j);

  vnl_matrix< TPrecision > result = vnl_matrix< TPrecision >(block, NDimension, NDimension);
  delete [] block;
  return result;
}

template<class TPrecision>
void
DiffusionTensor<TPrecision>
::SetVnlMatrix( const vnl_matrix<TPrecision> & mat)
{
  for(unsigned int i=0; i<NDimension; i++)
    {
    for(unsigned int j=i; j<NDimension; j++)
      (*this)[NDimension*i-(i+1)*i/2+j] = mat(i,j);
    }
}

template<class TPrecision>
DiffusionTensor<TPrecision>&
DiffusionTensor<TPrecision>
::Rotate (const typename Self::MatrixType& R)
{
  this->SetVnlMatrix ( this->GetRotate (R.GetVnlMatrix() ) );
  return *this;
}

template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::GetRotate (const typename Self::MatrixType& R) const
{
  Self res;
  res.SetVnlMatrix ( this->GetRotate (R.GetVnlMatrix() ) );
  return res;
}


template<class TPrecision>
vnl_matrix<TPrecision>
DiffusionTensor<TPrecision>
::GetRotate (const vnl_matrix<TPrecision> & R) const
{
  vnl_matrix<TPrecision> res;
  res = ( R*(this->GetVnlMatrix())*R.transpose() );
  return res;
}

template < class TPrecision>
template < class TVectorType >
void
DiffusionTensor<TPrecision>
::GetDWISamples ( TVectorType& dwisignal, const utl::NDArray<TPrecision,2>& gradients, const std::vector<TPrecision>& bValues) const
{
  utlAssert(gradients.Rows()==bValues.size(), "wrong size! gradients.Rows()="<< gradients.Rows() << ", bValues.size()="<<bValues.size());

  for ( int i = 0; i < gradients.Rows(); i += 1 ) 
    {
    TPrecision x=gradients(i,0), y=gradients(i,1), z=gradients(i,2);
    TPrecision quadratic = (*this)(0,0)*x*x + 2*(*this)(0,1)*x*y + 2*(*this)(0,2)*x*z +
        (*this)(1,1)*y*y + 2*(*this)(1,2)*y*z + (*this)(2,2)*z*z; 
    dwisignal[i] = std::exp(-1.0*bValues[i]*quadratic);
    }
}   // -----  end of method itkDiffusionTensor<TPrecision>::GetDWISignal  -----

template < class TPrecision >
template < class TVectorType >
void
DiffusionTensor<TPrecision>
::GetODFSamples ( TVectorType& odf, const utl::NDArray<TPrecision,2>& gradients, const int& odfOrder, const bool& isNormalize) const
{
  DiffusionTensor<TPrecision> D_inv = this->Inv();
  TPrecision D_det = this->GetDeterminant();
  for ( int i = 0; i < gradients.Rows(); i += 1 ) 
    {
    TPrecision x=gradients(i,0), y=gradients(i,1), z=gradients(i,2);
    TPrecision quadratic_inv = D_inv(0,0)*x*x + 2*D_inv(0,1)*x*y + 2*D_inv(0,2)*x*z +
        D_inv(1,1)*y*y + 2*D_inv(1,2)*y*z + D_inv(2,2)*z*z; 
    
    odf[i] = 1.0 / ( std::pow(D_det, 0.5) * std::pow(quadratic_inv,0.5*(odfOrder+1)) ) ;
    if (odfOrder==2)
      odf[i] *= 1.0/(4*M_PI);
    }

  // when odfOrder==2, it is naturally normalized.
  if (isNormalize && odfOrder!=2)
    {
    TPrecision normFactor = 4*M_PI*utl::GetSumOfVector<TVectorType>(odf,gradients.Rows()) / gradients.Rows();
    if (normFactor!=0)
      {
      for ( int i = 0; i < gradients.Rows(); i += 1 ) 
        odf[i] /= normFactor;
      }
    }
}   // -----  end of method itkDiffusionTensor<TPrecision>::GetDWISignal  -----

template < class TPrecision >
template < class TVectorType >
void
DiffusionTensor<TPrecision>
::GetEAPSamples ( TVectorType& eap, const utl::NDArray<TPrecision,2>& gradients, const std::vector<TPrecision>& rValues, const TPrecision& tau) const
{
  utlAssert(gradients.Rows()==rValues.size(), "wrong size! gradients.Rows()="<< gradients.Rows() << ", rValues.size()="<<rValues.size());

  DiffusionTensor<TPrecision> D_inv = this->Inv();
  TPrecision D_det = this->GetDeterminant();
  for ( int i = 0; i < gradients.Rows(); i += 1 ) 
    {
    TPrecision x=gradients(i,0), y=gradients(i,1), z=gradients(i,2);
    TPrecision quadratic_inv = D_inv(0,0)*x*x + 2*D_inv(0,1)*x*y + 2*D_inv(0,2)*x*z +
        D_inv(1,1)*y*y + 2*D_inv(1,2)*y*z + D_inv(2,2)*z*z; 
    eap[i] = 1.0 / (std::pow(4.0*M_PI*tau,(double)1.5)*std::sqrt(D_det) ) * std::exp(-1.0/(4.0*tau)*rValues[i]*rValues[i]*quadratic_inv);
    }
}   // -----  end of method itkDiffusionTensor<TPrecision>::GetDWISignal  -----

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetReturnToOrigin (const TPrecision tau) const
{
  RealValueType D_det = this->GetDeterminant();
  return 1.0 / (std::pow(4*M_PI*tau,(double)1.5)*std::sqrt(D_det) );
}   // -----  end of method itkDiffusionTensor<TPrecision>::GetDWISignal  -----

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetMeanSquaredDisplacement (const TPrecision tau) const
{
  return this->GetTrace()*2.0*tau;
}   // -----  end of method itkDiffusionTensor<TPrecision>::GetDWISignal  -----

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetGA () const
{
  vnl_diag_matrix<TPrecision> eigenValues(NDimension);
  vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
  this->GetEigenValuesVectors(eigenValues, eigenVectors);

  vnl_diag_matrix<TPrecision> logeigenValues(NDimension);
  for ( int i = 0; i < NDimension; i += 1 ) 
    logeigenValues[i] = std::log(eigenValues[i]);

  TPrecision meanEigenValues=0.0, ga=0.0;
  for ( int i = 0; i < NDimension; i += 1 ) 
    meanEigenValues += logeigenValues[i];
  meanEigenValues /= NDimension;
  for ( int i = 0; i < NDimension; i += 1 ) 
    ga += (logeigenValues[i]-meanEigenValues[i])*(logeigenValues[i]-meanEigenValues[i]);
  return std::sqrt(ga);
}

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetFA () const
{ 
  return this->GetFractionalAnisotropy();  
}

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetADC () const
{ 
  return this->GetTrace()/RealValueType(NDimension);
}

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetRA () const
{ 
  return this->GetRelativeAnisotropy(); 
}

template<class TPrecision>
bool
DiffusionTensor<TPrecision>
::IsPositive () const
{
  utlAssert (this->IsFinite(), "Tensor is not finite");
  EigenValuesArrayType eigenValues;
  this->ComputeEigenValues(eigenValues);
  if( eigenValues[0] <= NumericTraits<double>::Zero )
    return false;
  return true;
}

template<class TPrecision>
void
DiffusionTensor<TPrecision>
::Positivize () const
{
  if (IsPositive())
    return;
  
  vnl_diag_matrix<TPrecision> eigenValues(NDimension);
  vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
  this->GetEigenValuesVectors(eigenValues, eigenVectors);

  for ( int i = 0; i < NDimension; i += 1 ) 
    {
    if (eigenValues[i]<0)
      eigenValues[i]=0;
    }
  *this = eigenVectors* eigenValues.asMatrix()* eigenVectors.transpose();
}

template<class TPrecision>
bool
DiffusionTensor<TPrecision>
::IsZero () const
{
  for ( int i = 0; i < NDegreesOfFreedom; i += 1 ) 
    {
    if (std::abs((*this)[i])>1e-10)
      return false;
    }
  return true;
}

template<class TPrecision>
bool
DiffusionTensor<TPrecision>
::IsDiagonal () const
{
  for ( int i = 0; i < NDimension; i += 1 ) 
    for ( int j = i+1; j < NDimension; j += 1 ) 
      {
      if (std::abs((*this)(i,j))>1e-10)
        return false;
      }
  return true;
}

template<class TPrecision>
bool
DiffusionTensor<TPrecision>
::IsFinite () const
{
  return this->GetVnlMatrix().is_finite();
}

template<class TPrecision>
bool
DiffusionTensor<TPrecision>
::HasNans () const
{
  return this->GetVnlMatrix().has_nans();
}

template<class TPrecision>
TPrecision
DiffusionTensor<TPrecision>
::GetDeterminant () const
{
  vnl_matrix< TPrecision > M = this->GetVnlMatrix();
  return vnl_determinant( M );
}

template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::Log () const
{
  utlException (!this->IsFinite(), "Tensor is not finite");
  Self result;
  vnl_diag_matrix<TPrecision> eigenValues(NDimension);
  vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
  this->GetEigenValuesVectors(eigenValues, eigenVectors);
  for(unsigned int i=0;i<NDimension;i++)
    {
    utlException (eigenValues[i] < 0.0, "Negative eigenvalue encountered.");
    eigenValues[i] = vcl_log (eigenValues[i]);
    }
  
  vnl_matrix<TPrecision> vnl_result = eigenVectors * eigenValues.asMatrix() * eigenVectors.transpose();
  result.SetVnlMatrix ( vnl_result );
  return result;
}

template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::Exp () const
{
  utlAssert(this->IsFinite(),"Tensor is not finite");
  Self result;
  vnl_diag_matrix<TPrecision> eigenValues(NDimension);
  vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
  this->GetEigenValuesVectors(eigenValues, eigenVectors);
  for(unsigned int i=0;i<NDimension;i++)
    eigenValues[i] = vcl_exp (eigenValues[i]);

  vnl_matrix<TPrecision> vnl_result = eigenVectors * eigenValues.asMatrix() * eigenVectors.transpose();
  result.SetVnlMatrix ( vnl_result );
  return result;
}

template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::Pow (const double& n) const
{
  utlAssert (this->IsFinite(), "Tensor is not finite");
  Self result;

  vnl_diag_matrix<TPrecision> eigenValues(NDimension);
  vnl_matrix<TPrecision> eigenVectors(NDimension, NDimension);
  this->GetEigenValuesVectors(eigenValues, eigenVectors);
  // std::cout << *this << std::endl << std::flush;
  // std::cout << "eigenValues = " << eigenValues << std::endl << std::flush;
  // std::cout << "eigenVectors = " << eigenVectors << std::endl << std::flush;

  // typedef vnl_symmetric_eigensystem< TPrecision >  SymEigenSystemType;
  // SymEigenSystemType eig (this->GetVnlMatrix());
  // std::cout << "eig.D = " << eig.D << std::endl << std::flush;
  // std::cout << "eig.V = " << eig.V << std::endl << std::flush;

  for(unsigned int i=0;i<NDimension;i++)
    {
    utlException(n<0 && std::abs(eigenValues[i])<=1e-10, "pow="<<n<<", eigenValues["<<i<<"]="<<eigenValues[i] << " close to zero");
    eigenValues[i] = static_cast<TPrecision> (vcl_pow (static_cast<double>(eigenValues[i]), n));
    }

  vnl_matrix<TPrecision> vnl_result = eigenVectors * eigenValues.asMatrix() * eigenVectors.transpose();
  result.SetVnlMatrix ( vnl_result );
  return result;
}

template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::Inv (void) const
{
  return this->Pow (-1.0);
}


template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::Sqrt (void) const
{
  return this->Pow (0.5);
}


template<class TPrecision>
DiffusionTensor<TPrecision>
DiffusionTensor<TPrecision>
::InvSqrt () const
{
  return this->Pow (-0.5);
}

template<typename TPrecision>
TPrecision 
DiffusionTensor<TPrecision>
::EuclideanDistance(const DiffusionTensor<TPrecision>& X) const 
{
  return (*this - X).GetInnerScalarProduct();
}

template<typename TPrecision>
TPrecision 
DiffusionTensor<TPrecision>
::KLDistance(const DiffusionTensor<TPrecision>& X) const 
{
  const DiffusionTensor<TPrecision> ithis = this->Inv();
  const DiffusionTensor<TPrecision> iX = X.Inv();
  const TPrecision tr = (ithis*X + iX*(*this)).GetTrace();
  return 0.5*std::sqrt(tr - 2*NDimension);
}

template<typename TPrecision>
TPrecision 
DiffusionTensor<TPrecision>
::GeodesicDistance(const DiffusionTensor<TPrecision>& X) const 
{
  const DiffusionTensor<TPrecision> isqrtT = this->InvSqrt();
  const DiffusionTensor<TPrecision> A = (isqrtT * X * isqrtT).Log();
  return std::sqrt(0.5*(A*A).GetTrace());
}

template<typename TPrecision>
TPrecision 
DiffusionTensor<TPrecision>
::LogEucDistance(const DiffusionTensor<TPrecision>& X) const 
{
  DiffusionTensor<TPrecision> logThis = this->Log();
  DiffusionTensor<TPrecision> logX = X.Log();
  return (logThis-logX).GetInnerScalarProduct();
}


}

#endif 
