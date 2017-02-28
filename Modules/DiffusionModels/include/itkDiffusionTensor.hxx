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

  eigenValues.set_size(NDimension);
  for ( int i = 0; i < NDimension; i += 1 ) 
    eigenValues[i] = values[i];
  eigenVectors = vectors.GetVnlMatrix().transpose();
}

template<class TPrecision>
template<class TArrayType, class TMatrixType >
void
DiffusionTensor<TPrecision>
::GetEigenValuesVectorsAnalytic(TArrayType& eigenValues, TMatrixType& eigenVectors) const
{
  if (IsDiagonal(1e-10))
    {
    for ( int i = 0; i < 3; ++i ) 
      {
      eigenValues[i] = (*this)(i,i);
      eigenVectors(i,i) = 1.0;
      for ( int j = 0; j < 3; ++j ) 
        {
        if (i!=j)
          eigenVectors(i,j) = 0.0;
        }
      }

    // sort three values in ascending order
    if (eigenValues[0]>eigenValues[1])
      {
      std::swap(eigenValues[0], eigenValues[1]);
      for ( int i = 0; i < 3; ++i ) 
        std::swap(eigenVectors(i,0), eigenVectors(i,1));
      }
    if (eigenValues[1]>eigenValues[2])
      {
      std::swap(eigenValues[1], eigenValues[2]);
      for ( int i = 0; i < 3; ++i ) 
        std::swap(eigenVectors(i,1), eigenVectors(i,2));
      }
    if (eigenValues[0]>eigenValues[1])
      {
      std::swap(eigenValues[0], eigenValues[1]);
      for ( int i = 0; i < 3; ++i ) 
        std::swap(eigenVectors(i,0), eigenVectors(i,1));
      }
    return;
    }

  std::vector<double> detCoef(4);
  std::vector<std::complex<double> > lambdaVec;

  const TPrecision* p = this->GetDataPointer();
  detCoef[0] = p[2]*p[2]*p[3] - 2.0*p[1]*p[2]*p[4] + p[0]*p[4]*p[4] + p[1]*p[1]*p[5] - p[0]*p[3]*p[5];
  detCoef[1] = -p[1]*p[1] - p[2]*p[2] + p[0]*p[3] - p[4]*p[4] + p[0]*p[5] + p[3]*p[5]; 
  detCoef[2] = -p[0] - p[3] - p[5];
  detCoef[3] = 1;
  lambdaVec = utl::PolynomialRoot(detCoef);

  for ( int i = 0; i < 3; ++i ) 
    eigenValues[i] = std::real(lambdaVec[i]);

  // sort three values in ascending order
  if (eigenValues[0]>eigenValues[1])
    std::swap(eigenValues[0], eigenValues[1]);
  if (eigenValues[1]>eigenValues[2])
    std::swap(eigenValues[1], eigenValues[2]);
  if (eigenValues[0]>eigenValues[1])
    std::swap(eigenValues[0], eigenValues[1]);
  
  bool numericalZero =false;
  if (std::fabs(p[2])<1e-5)
    numericalZero = true;

  if (!numericalZero)
    {
    double norm, v1,v2,v3;
    // double p2=p[2]!=0?p[2]:1e-10, p1=p[1]!=0?p[1]:1e-10;
    double p2=p[2], p1=p[1];

    double cebf = -p2*p[4]+p1*p[5];
    double cdbe = -p2*p[3]+p1*p[4];
    for ( int i = 0; i < 3; ++i ) 
      {
      v1 = -(p[5]-eigenValues[i])*(cdbe+p2*eigenValues[i]) + p[4]*(cebf-p1*eigenValues[i]);
      v2 = -p2*(cebf-p1*eigenValues[i]);
      v3 = p2*(cdbe+p2*eigenValues[i]);
      norm = std::sqrt(v1*v1 + v2*v2 + v3*v3);
      if (std::fabs(p2)<1e-5 || std::fabs(cdbe+p2*eigenValues[i])<1e-5 )
        {
        // NOTE: when p[2] or cdbe+p2*eigenValues[i] equal 0, we have divide by zero problem.
        numericalZero=true;
        break;
        }
      eigenVectors(0,i) = v1/norm;
      eigenVectors(1,i) = v2/norm;
      eigenVectors(2,i) = v3/norm;
      }
    }

  if (numericalZero)
    {
    vnl_diag_matrix<TPrecision> val;
    vnl_matrix<double> vec;
    GetEigenValuesVectors(val, vec);
    for ( int i = 0; i < 3; ++i ) 
      for ( int j = 0; j < 3; ++j ) 
        eigenVectors(i,j) = vec(i,j);
    }
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
    this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
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
    this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
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
double
DiffusionTensor<TPrecision>
::GetQuadraticForm ( const double x, const double y, const double z) const
{
  return (*this)[0]*x*x + 2*(*this)[1]*x*y + 2*(*this)[2]*x*z +
    (*this)[3]*y*y + 2*(*this)[4]*y*z + (*this)[5]*z*z; 
}

template < class TPrecision>
template < class TVectorType >
void
DiffusionTensor<TPrecision>
::GetSphericalSamples ( TVectorType& samples, const utl::NDArray<TPrecision,2>& gradients) const
{
  for ( int i = 0; i < gradients.Rows(); i += 1 ) 
    {
    TPrecision x=gradients(i,0), y=gradients(i,1), z=gradients(i,2);
    double quadratic = GetQuadraticForm(x,y,z);
    samples[i] = GetQuadraticForm(x,y,z);
    }
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
    double quadratic = GetQuadraticForm(x,y,z);
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
  this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);

  vnl_diag_matrix<TPrecision> logeigenValues(NDimension);
  for ( int i = 0; i < NDimension; i += 1 ) 
    logeigenValues[i] = std::log(eigenValues[i]);

  TPrecision meanEigenValues=0.0, ga=0.0;
  for ( int i = 0; i < NDimension; i += 1 ) 
    meanEigenValues += logeigenValues[i];
  meanEigenValues /= NDimension;
  for ( int i = 0; i < NDimension; i += 1 ) 
    {
    double tmp = logeigenValues[i]-meanEigenValues;
    ga += tmp*tmp;
    }
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

template < class TPrecision >
typename DiffusionTensor<TPrecision>::RealValueType
DiffusionTensor<TPrecision>
::GetMODE () const
{ 
  double md = GetMD();
  
  Self dTensor(*this);
  dTensor[0] -= md;
  dTensor[3] -= md;
  dTensor[5] -= md;

  double dtNorm = dTensor.GetNorm();
  dTensor /= dtNorm;
  return 3.0*std::sqrt(6.0)*dTensor.GetDeterminant();
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
  this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);

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
::IsZero (const double eps) const
{
  for ( int i = 0; i < NDegreesOfFreedom; i += 1 ) 
    {
    if (std::fabs((*this)[i])>eps)
      return false;
    }
  return true;
}

template<class TPrecision>
bool
DiffusionTensor<TPrecision>
::IsDiagonal (const double eps) const
{
  for ( int i = 0; i < NDimension; i += 1 ) 
    for ( int j = i+1; j < NDimension; j += 1 ) 
      {
      if (std::fabs((*this)(i,j))>eps)
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
double
DiffusionTensor<TPrecision>
::GetDeterminant () const
{
  const TPrecision* p = this->GetDataPointer();
  return -p[2]*p[2]*p[3] + 2.0*p[1]*p[2]*p[4] - p[0]*p[4]*p[4] - p[1]*p[1]*p[5] + p[0]*p[3]*p[5];
}

template<class TPrecision>
double
DiffusionTensor<TPrecision>
::GetNorm () const
{
  const TPrecision* p = this->GetDataPointer();
  return std::sqrt(p[0]*p[0]+2.0*p[1]*p[1]+2.0*p[2]*p[2] + p[3]*p[3]+ 2.0*p[4]*p[4] + p[5]*p[5]);
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
  this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
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
  this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
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
  this->GetEigenValuesVectorsAnalytic(eigenValues, eigenVectors);
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
  double det = GetDeterminant();
  utlGlobalException(det==0, "det=0, cannot be inverted");
  double detInv=1.0/det;
  const TPrecision* p = this->GetDataPointer();
  DiffusionTensor<TPrecision> result;
  result[0] = (-p[4]*p[4] + p[3]*p[5]) *detInv;
  result[1] = ( p[2]*p[4] - p[1]*p[5]) *detInv;
  result[2] = (-p[2]*p[3] + p[1]*p[4]) *detInv;
  result[3] = (-p[2]*p[2] + p[0]*p[5]) *detInv;
  result[4] = ( p[1]*p[2] - p[0]*p[4]) *detInv;
  result[5] = (-p[1]*p[1] + p[0]*p[3]) *detInv;
  return result;
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
