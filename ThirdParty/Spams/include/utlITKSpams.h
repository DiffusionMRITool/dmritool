/**
 *       @file  utlITKSpams.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "01-02-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlITKSpams_h
#define __utlITKSpams_h



#include "linalg.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace spams
{

template <class T>
void
VnlMatrixToMatrix ( const vnl_matrix<T>& matVnl, Matrix<T>& matSpams )
{
  matSpams.resize(matVnl.rows(), matVnl.columns());
  // for ( int i = 0; i < matVnl.rows(); i += 1 ) 
  //   for ( int j = 0; j < matVnl.columns(); j += 1 ) 
  //     matSpams(i,j) = matVnl(i,j);
  T const* matVnl_data = matVnl.data_block();
  T* matSpams_data = matSpams.rawX();
  if (matVnl.rows() >= matVnl.columns())
    {
    for ( int j = 0; j < matVnl.columns(); j += 1 ) 
      cblas_copy<T>(matVnl.rows(), (T*)matVnl_data+j, matVnl.columns(), matSpams_data+j*matVnl.rows(),1);
    }
  else
    {
    for ( int i = 0; i < matVnl.rows(); i += 1 ) 
      cblas_copy<T>(matVnl.columns(), (T*)matVnl_data+i*matVnl.columns(), 1, matSpams_data+i,matVnl.rows());
    }
}

template <class T>
void
MatrixToVnlMatrix ( const Matrix<T>& matSpams, vnl_matrix<T>& matVnl )
{
  matVnl.set_size(matSpams.m(), matSpams.n());
  // for ( int i = 0; i < mat.m(); i += 1 ) 
  //   for ( int j = 0; j < mat.n(); j += 1 ) 
  //     matVnl(i,j) = matSpams(i,j);
  T* matVnl_data = matVnl.data_block();
  const T * matSpams_data = matSpams.X();
  if (matVnl.rows() >= matVnl.columns())
    {
    for ( int j = 0; j < matVnl.columns(); j += 1 ) 
      cblas_copy<T>(matVnl.rows(), (T*)matSpams_data+j*matVnl.rows(),1, matVnl_data+j, matVnl.columns());
    }
  else
    {
    for ( int i = 0; i < matVnl.rows(); i += 1 ) 
      cblas_copy<T>(matVnl.columns(), (T*)matSpams_data+i,matVnl.rows(), matVnl_data+i*matVnl.columns(), 1);
    }
}

template <class T>
void
VnlVectorToVector ( const vnl_vector<T>& v, Vector<T>& vec )
{
  vec.resize(v.size());
  // for ( int i = 0; i < v.size(); i += 1 ) 
  //   vec[i] = v[i];
  cblas_copy<T>(v.size(), v.data_block(),1, vec.rawX(), 1);
}

template <class T>
void
VectorToVnlVector ( const Vector<T>& v, vnl_vector<T>& vec )
{
  vec.set_size(v.n());
  // for ( int i = 0; i < v.n(); i += 1 ) 
  //   vec[i] = v[i];
  cblas_copy<T>(vec.size(), v.rawX(),1, vec.data_block(), 1);
}

template <class T>
void 
SpMatrixToVnlMatrix ( const SpMatrix<T>& mat, vnl_matrix<T>& result )
{
  result.set_size(mat.m(), mat.n());
  result.fill(0.0);
  int* pB = mat.pB();
  int* pE = mat.pE();
  int* r = mat.r();
  for ( int i = 0; i < mat.n(); i += 1 ) 
    for (int j = pB[i]; j<pE[i]; ++j)
      result(r[j], i) = mat.v(j);
}


}

#endif 

