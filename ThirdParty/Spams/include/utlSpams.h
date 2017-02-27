/**
 *       @file  utlSpams.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-29-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlSpams_h
#define __utlSpams_h

#include "linalg.h"
#include "utlNDArray.h"

namespace spams
{

template <class T>
void
UtlMatrixToMatrix ( const utl::NDArray<T,2>& matUtl, Matrix<T>& matSpams )
{
  matSpams.resize(matUtl.Rows(), matUtl.Columns());
  T const* matUtl_data = matUtl.GetData();
  T* matSpams_data = matSpams.rawX();
  if (matUtl.Rows() >= matUtl.Columns())
    {
    for ( int j = 0; j < matUtl.Columns(); j += 1 ) 
      cblas_copy<T>(matUtl.Rows(), (T*)matUtl_data+j, matUtl.Columns(), matSpams_data+j*matUtl.Rows(),1);
    }
  else
    {
    for ( int i = 0; i < matUtl.Rows(); i += 1 ) 
      cblas_copy<T>(matUtl.Columns(), (T*)matUtl_data+i*matUtl.Columns(), 1, matSpams_data+i,matUtl.Rows());
    }
}

template <class T>
void
MatrixToUtlMatrix ( const Matrix<T>& matSpams, utl::NDArray<T,2>& matUtl )
{
  matUtl.ReSize(matSpams.m(), matSpams.n());
  T* matUtl_data = matUtl.GetData();
  const T * matSpams_data = matSpams.X();
  if (matUtl.Rows() >= matUtl.Columns())
    {
    for ( int j = 0; j < matUtl.Columns(); j += 1 ) 
      cblas_copy<T>(matUtl.Rows(), (T*)matSpams_data+j*matUtl.Rows(),1, matUtl_data+j, matUtl.Columns());
    }
  else
    {
    for ( int i = 0; i < matUtl.Rows(); i += 1 ) 
      cblas_copy<T>(matUtl.Columns(), (T*)matSpams_data+i,matUtl.Rows(), matUtl_data+i*matUtl.Columns(), 1);
    }
}

template <class T>
void
UtlVectorToVector ( const utl::NDArray<T,1>& v, Vector<T>& vec )
{
  vec.resize(v.size());
  cblas_copy<T>(v.Size(), v.GetData(),1, vec.rawX(), 1);
}

template <class T>
void
VectorToUtlVector ( const Vector<T>& v, utl::NDArray<T,1>& vec )
{
  vec.ReSize(v.n());
  cblas_copy<T>(vec.Size(), v.rawX(),1, vec.GetData(), 1);
}

template <class T>
void 
SpMatrixToUtlMatrix ( const SpMatrix<T>& mat, utl::NDArray<T,2>& result )
{
  result.ReSize(mat.m(), mat.n());
  result.Fill(0.0);
  int* pB = mat.pB();
  int* pE = mat.pE();
  int* r = mat.r();
  for ( int i = 0; i < mat.n(); i += 1 ) 
    for (int j = pB[i]; j<pE[i]; ++j)
      result(r[j], i) = mat.v(j);
}


}


#endif 
