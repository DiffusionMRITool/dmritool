/**
 *       @file  utlVNL.h
 *      @brief  help functions for VNL
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-16-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlVNL_h
#define __utlVNL_h


#include "utlCore.h"
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_svd.h>

namespace utl 
{

/** @addtogroup utlHelperFunctions
@{ */

template <class T>
vnl_matrix<T> 
SphericalToCartesian ( const vnl_matrix<T>& in )
{
  utlAssert(in.columns()==3 || in.rows()==3, "wrong dimension");
  vnl_matrix<T> out(in);
  if (in.columns()==3)
    {
    for ( int i = 0; i < in.rows(); i += 1 ) 
      utl::spherical2Cartesian(in(i,0), in(i,1), in(i,2), out(i,0), out(i,1), out(i,2));
    }
  else
    {
    for ( int i = 0; i < in.columns(); i += 1 ) 
      utl::spherical2Cartesian(in(0,i), in(1,i), in(2,i), out(0,i), out(1,i), out(2,i));
    }
  return out;
}

template <class T>
vnl_matrix<T> 
CartesianToSpherical ( const vnl_matrix<T>& in )
{
  utlAssert(in.columns()==3 || in.rows()==3, "wrong dimension");
  vnl_matrix<T> out(in);
  if (in.columns()==3)
    {
    for ( int i = 0; i < in.rows(); i += 1 ) 
      utl::cartesian2Spherical(in(i,0), in(i,1), in(i,2), out(i,0), out(i,1), out(i,2));
    }
  else
    {
    for ( int i = 0; i < in.columns(); i += 1 ) 
      utl::cartesian2Spherical(in(0,i), in(1,i), in(2,i), out(0,i), out(1,i), out(2,i));
    }
  return out;
}

/** min-max normalization for each column of the given vnl_matrix  */
template <class T>
vnl_matrix<T> 
NormalizeMinMax ( const vnl_matrix<T>& matrix )
{
  std::vector<T> vec(matrix.rows());
  vnl_matrix<T> result(matrix);
  for ( int i = 0; i < matrix.columns(); i += 1 ) 
    {
    for ( int j = 0; j < matrix.rows(); j += 1 ) 
      vec[j] = matrix(j,i);
    vec = utl::NormalizeMinMax(vec);
    for ( int j = 0; j < matrix.rows(); j += 1 ) 
      result(j,i) = vec[j];
    }
  return result;
}

template <class T>
vnl_matrix<T>
ConnectVnlMatrix ( const vnl_matrix<T>& m1, const vnl_matrix<T>& m2, const bool isConnectRow )
{
  if (isConnectRow)
    {
    utlException(m1.columns()!=m2.columns(), "wrong column size! m1.columns()="<<m1.columns()<<", m2.columns()"<<m2.columns());
    vnl_matrix<T> result(m1.rows()+m2.rows(), m1.columns());
    int m1Rows = m1.rows();
    for ( int j = 0; j < m1.columns(); j += 1 ) 
      {
      for ( int i = 0; i < m1.rows(); i += 1 ) 
        result(i,j) = m1(i,j);
      for ( int i = 0; i < m2.rows(); i += 1 ) 
        result(i+m1Rows,j) = m2(i,j);
      }
    return result;
    }
  else
    {
    utlException(m1.rows()!=m2.rows(), "wrong row size! m1.rows()="<<m1.rows()<<", m2.rows()"<<m2.rows());
    vnl_matrix<T> result(m1.rows(), m1.columns()+m2.columns());
    int m1Columns = m1.columns();
    for ( int i = 0; i < m1.rows(); i += 1 ) 
      {
      for ( int j = 0; j < m1.columns(); j += 1 ) 
        result(i,j) = m1(i,j);
      for ( int j = 0; j < m2.columns(); j += 1 ) 
        result(i,j+m1Columns) = m2(i,j);
      }
    return result;
    }
}

template <class T>
vnl_vector<T>
ConnectVnlVector ( const vnl_vector<T>& m1, const vnl_vector<T>& m2 )
{
  vnl_vector<T> result(m1.size()+m2.size());
  int m1Size = m1.size();
  for ( int i = 0; i < m1Size; i += 1 ) 
    result[i] = m1[i];
  for ( int i = 0; i < m2.size(); i += 1 ) 
    result[i+m1Size] = m2[i];
  return result;
}

template <class T>
vnl_vector<T>
GetAbsoluteVnlVector ( const vnl_vector<T>& vec )
{
  vnl_vector<T> result(vec);
  for ( int i = 0; i < vec.size(); i += 1 ) 
    result[i] = std::fabs(vec[i]);
  return result;
}

template <class T>
T
GetMedianVnlVector(const vnl_vector<T>& values)
{
  vnl_vector<T> vec(values);
  const unsigned int N = vec.size();
  std::nth_element(vec.begin(),vec.begin()+N/2,vec.end());
  return vec[N/2];
}

template <class T>
std::vector<T> 
GetVnlVectorStats ( const vnl_vector<T>& values )
{
  return utl::GetContainerStats(values.begin(), values.end());
}

template <class T>
std::vector<T> 
GetVnlMatrixStats ( const vnl_matrix<T>& values )
{
  return utl::GetContainerStats(values.begin(), values.end());
}

/** get the minimal angle from a given evenly distributed gradients in Cartesian format */
template <class T>
inline double
GetMinAngle( const vnl_matrix<T>& grad ) // cartesian
{
  utlException(grad.rows()<2 || grad.columns()!=3, "grad has wrong size!");
  double angle = 90.0;  

  /* finding minimum angle between samplings */
  for(unsigned int i = 1; i < grad.rows(); i++) 
    {
    double dot = grad(0,0)*grad(i,0) + grad(0,1)*grad(i,1) + grad(0,2)*grad(i,2);

    if(dot >= -1.0 - M_EPS && dot <= -1.0 + M_EPS)
      dot = 180;
    else if(dot <= 1.0 + M_EPS && dot >= 1.0 - M_EPS)
      dot = 0;
    else 
      dot = 180*std::acos(dot)/M_PI;

    if(dot < angle)
      angle = dot;
    }
  return angle;
}

/** Get the pseudoinverse of a symmetric vnl_matrix. The small eigenvalues are set to zero  */
template <class T>
vnl_matrix<T> 
GetVnlSymmericMatrixPInverse ( const vnl_matrix<T>& mat, const double eps=1e-8)
{
  utlException(mat.rows()!=mat.columns(), "wrong size! mat.rows()="<<mat.rows()<<", mat.columns()="<<mat.columns());
  typedef vnl_symmetric_eigensystem< T >  SymEigenSystemType;
  SymEigenSystemType eig (mat);
  unsigned n = eig.D.rows();
  vnl_diag_matrix<T> invD(eig.D);
  for (unsigned i=0; i<n; ++i)
    {
    if ( eig.D(i,i)>eps || eig.D(i,i)<-eps ) 
      invD(i,i) = 1.0/eig.D(i,i);
    else
      invD(i, i) = 0.0;
    }
  return eig.V * invD * eig.V.transpose();
}

/** Get the pseudoinverse of a vnl_matrix. The small eigenvalues are set to zero  */
template <class T>
vnl_matrix<T> 
GetVnlMatrixPInverse ( const vnl_matrix<T>& mat, const double eps=1e-8)
{
  vnl_svd<T> svd(mat);
  svd.zero_out_absolute(eps);
  return svd.pinverse();
}

template <class T>
void
SaveVnlMatrix ( const vnl_matrix<T>& matrix, const std::string file )
{
  int NumberRows = matrix.rows(), NumberColumns=matrix.columns();
  utl::SaveMatrix<vnl_matrix<T> >(matrix, NumberRows, NumberColumns, file);
  return;
}

template <class T>
void 
PrintVnlMatrixStats (  const vnl_matrix<T>& matrix, const std::string str="", const char* separate=" ", std::ostream& os=std::cout )
{
  int NumberRows = matrix.rows(), NumberColumns=matrix.columns();
  utl::PrintMatrixStats<vnl_matrix<T> >(matrix, NumberRows, NumberColumns, str, separate, os);
}

template <class T>
void 
PrintVnlMatrix (  const vnl_matrix<T>& matrix, const std::string str="", const char* separate=" ", std::ostream& os=std::cout )
{
  int NumberRows = matrix.rows(), NumberColumns=matrix.columns();
  utl::PrintMatrix<vnl_matrix<T> >(matrix, NumberRows, NumberColumns, str, separate, os);
}

template <class T>
void 
PrintVnlVector ( const vnl_vector<T>& vec, const std::string str="", const char* separate=" ", std::ostream& os=std::cout)
{
  int NSize = vec.size();
  utl::PrintContainer(vec.begin(), vec.end(), str, separate, os);
}

template <class T>
vnl_vector<T>
GetValuesVnlVector ( const vnl_vector<T>& vec, const int colstart, const int n )
{
  vnl_vector<T> result(n);
  for ( int i = 0; i < n; i += 1 ) 
    result[i] = vec[i+colstart];
  return result;
}

template <class T>
vnl_vector<T>
GetValuesVnlVector ( const vnl_vector<T>& vec, const std::vector<int>& index )
{
  vnl_vector<T> result(index.size());
  for ( int i = 0; i < index.size(); i += 1 ) 
    result[i] = vec[index[i]];
  return result;
}

template <class T>
void
SetValuesVnlVector ( const vnl_vector<T>& subvec, vnl_vector<T>& vec, const int colstart )
{
  utlException(subvec.size()+colstart>vec.size(), "wrong size! subvec.size()="<<subvec.size() << ", vec.size()="<<vec.size() );
  for ( int i = 0; i < subvec.size(); i += 1 ) 
    vec[i+colstart] = subvec[i];
}

template <class T>
void
SetValuesVnlVector ( const vnl_vector<T>& subvec, vnl_vector<T>& vec, const std::vector<int>& index )
{
  utlException(subvec.size()!=index.size(), "wrong size!");
  for ( int i = 0; i < index.size(); i += 1 ) 
    vec[index[i]] = subvec[i];
}

template <class T>
vnl_matrix<T> 
GetRowsVnlMatrix ( const vnl_matrix<T>& mat, const std::vector<int>& index )
{
  utlException(index.size()>mat.rows(), "wrong size!");
  vnl_matrix<T> result(index.size(), mat.columns());
  for ( int i = 0; i < index.size(); i += 1 ) 
    for ( int j = 0; j < mat.columns(); j += 1 ) 
      result(i, j) = mat(index[i], j);
  return result;
}

template <class T>
vnl_matrix<T> 
GetColumnsVnlMatrix ( const vnl_matrix<T>& mat, const std::vector<int>& index )
{
  utlException(index.size()>mat.columns(), "wrong size!");
  vnl_matrix<T> result(mat.rows(), index.size());
  for ( int i = 0; i < mat.rows(); i += 1 ) 
    for ( int j = 0; j < index.size(); j += 1 ) 
      result(i, j) = mat(i, index[j]);
  return result;
}

template <class T>
void
SetRowsVnlMatrix ( const vnl_matrix<T>& submat, vnl_matrix<T>& mat, const std::vector<int>& index )
{
  utlException(submat.columns()!=mat.columns(), "wrong size!");
  utlException(submat.rows()!=index.size(), "wrong size!");
  for ( int i = 0; i < index.size(); i += 1 ) 
    for ( int j = 0; j < mat.columns(); j += 1 ) 
      mat(index[i], j) = submat(i,j);
}

template <class T>
void
SetColumnsVnlMatrix ( const vnl_matrix<T>& submat, vnl_matrix<T>& mat, const std::vector<int>& index )
{
  utlException(submat.rows()!=mat.rows(), "wrong size!");
  utlException(submat.columns()!=index.size(), "wrong size!");
  for ( int i = 0; i < mat.rows(); i += 1 ) 
    for ( int j = 0; j < index.size(); j += 1 ) 
      mat(i, index[j]) = submat(i,j);
}

template <class T>
std::vector<T> 
VnlVectorToStdVector ( const vnl_vector<T>& vec )
{
  std::vector<T> v(vec.size());
  for ( int i = 0; i < vec.size(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
vnl_vector<T> 
StdVectorToVnlVector ( const std::vector<T>& vec )
{
  vnl_vector<T> v(vec.size());
  for ( int i = 0; i < vec.size(); i += 1 ) 
    v[i] = vec[i];
  return v;
}

template <class T>
vnl_matrix<T> 
GetDiagonalMatrix ( const vnl_vector<T>& vec )
{
  vnl_matrix<T> mat(vec.size(), vec.size());
  mat.fill(0.0);
  mat.set_diagonal(vec);
  return mat;
}

template <class T>
bool 
IsMatrixSymmetric ( const vnl_matrix<T>& mat, const double eps=1e-10 )
{
  if (mat.rows()!=mat.columns())
    return false;
  for ( int i = 1; i < mat.rows(); i += 1 ) 
    {
    for ( int j = 0; j < i; j += 1 ) 
      {
      if ( std::fabs(mat(i,j)- mat(j,i))>eps )
        return false;
      }
    }
  return true;
}

    /** @} */

}


#endif 
