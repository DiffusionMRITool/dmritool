/**
 *       @file  utlVNLIO.h
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "08-27-2014
 *    Revision  1.0
 *    Compiler  gcc/g++
 *   Copyright  Copyright (c) 2014, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __utlVNLIO_h
#define __utlVNLIO_h

#include "utlCore.h"
#include "utlVector.h"
#include "utlNDArray.h"
#include "utlVNL.h"
#include <vnl/vnl_matrix_fixed.h>

namespace utl
{

/** @addtogroup utlHelperFunctions
@{ */

template <class T>
void 
VnlMatrixToUtlMatrix( const vnl_matrix<T>& mat, utl::NDArray<T,2>& matUtl )
{
  matUtl.ReSize(mat.rows(),mat.cols());
  matUtl.CopyData((T* const)mat.data_block(), mat.rows(), mat.cols() );
}

template <class T, unsigned rows, unsigned cols>
void 
VnlMatrixToUtlMatrix( const vnl_matrix_fixed<T, rows, cols>& mat, utl::NDArray<T,2>& matUtl )
{
  matUtl.ReSize(mat.rows(),mat.cols());
  matUtl.CopyData((T* const)mat.data_block(), mat.rows(), mat.cols() );
}

template <class T>
utl::NDArray<T,2> 
VnlMatrixToUtlMatrix( const vnl_matrix<T>& mat )
{
  utl::NDArray<T,2> result;
  VnlMatrixToUtlMatrix(mat, result);
  return result;
}

template <class T>
void 
UtlMatrixToVnlMatrix( const NDArray<T,2>& mat, vnl_matrix<T>& matVnl )
{
  matVnl.set_size(mat.Rows(), mat.Cols());
  matVnl.copy_in(mat.GetData());
}

template <class T>
vnl_matrix<T> 
UtlMatrixToVnlMatrix( const NDArray<T,2>& mat )
{
  vnl_matrix<T> result(mat.GetData(), mat.Rows(), mat.Cols());
  return result;
}

template <class T>
utl::NDArray<T,1> 
VnlVectorToUtlVector( const vnl_vector<T>& vec )
{
  utl::NDArray<T,1> result(vec.data_block(), vec.size());
  return result;
}

template <class T>
vnl_vector<T> 
UtlVectorToVnlVector( const NDArray<T,1>& vec )
{
  vnl_vector<T> result(vec.GetData(), vec.Size());
  return result;
}

    /** @} */
}

#endif 
